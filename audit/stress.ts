/**
 * What the tool does when the person driving it makes a mistake, or the file is not what it says.
 *
 * EVERY OTHER HARNESS HERE ASKS WHETHER THE SCIENCE IS RIGHT. This one asks whether the tool
 * survives its users. Those are different questions and the second has a different failure mode: a
 * scientific error shows up as a wrong number in a calibration, while an operator error shows up as
 * a perfectly formed report about the wrong thing. No calibration arm would ever catch it, because
 * every one of them loads the correct file.
 *
 * THE STANDARD IS THE ONE THE REST OF THE TOOL IS HELD TO. For each scenario the tool must either
 * be RIGHT or REFUSE VISIBLY. A confident wrong answer fails. So does a crash with no explanation,
 * because an operator who sees a stack trace cannot tell a bad file from a broken tool.
 *
 * Real arrays throughout, from GSE148488, with roles known from the manifest and from genetics.
 * Nothing is constructed except the corrupted files, and those are corruptions OF real arrays.
 *
 * Run: OM_TRIOS=<dir> node --experimental-strip-types --max-old-space-size=6144 audit/stress.ts
 */
import { readFileSync, writeFileSync, existsSync, mkdtempSync, rmSync } from 'node:fs'
import { gunzipSync, gzipSync } from 'node:zlib'
import { tmpdir } from 'node:os'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const defects = await import(`${W}defects.ts`)

const DIR = process.env.OM_TRIOS
if (!DIR || !existsSync(DIR)) {
  console.log('stress: OM_TRIOS is not set to a readable directory. NOT RUN.')
  console.log('  Skipped rather than quietly passing: a stress suite that prints nothing reads')
  console.log('  exactly like one where nothing went wrong.')
  process.exit(0)
}

// ---------------------------------------------------------------- the corpus, and who is who
//
// Ground truth is the GEO manifest plus one fact of biology: a sperm donor is male and an egg donor
// is female. Neither comes from this tool.
const SPERM = ['GSM4472397', 'GSM4472398', 'GSM4472399', 'GSM4472400']
const EGG = ['GSM4472407', 'GSM4472415', 'GSM4472416', 'GSM4472417',
  'GSM4472418', 'GSM4472419', 'GSM4472420']
/** Excluded as a parental genotype: het 0.3455 at 82% call rate, a mixed or contaminated sample. */
const CONTAMINATED = 'GSM4472408'

interface Rec {
  gsm: string; title: string; material: string; role: string
  father: string; mother: string; complete: boolean
}
const manifest: Rec[] = (() => {
  const lines = readFileSync(join(DIR, 'trio_manifest_full.csv'), 'utf8').trim().split('\n')
  const head = lines[0].split(',')
  const col = (n: string) => head.indexOf(n)
  return lines.slice(1).map((line) => {
    const f = line.split(',')
    const first = (s: string) => (s ?? '').split(';').filter((x) => x && x !== CONTAMINATED)[0] ?? ''
    return {
      gsm: f[col('gsm')], title: f[col('title')], material: f[col('material')],
      role: f[col('role')],
      father: first(f[col('father_gsms')]), mother: first(f[col('mother_gsms')]),
      complete: f[col('complete_trio')] === 'True',
    }
  })
})()
const has = (g: string) => !!g && existsSync(join(DIR, `${g}.probes.gz`))

const textOf = (gsm: string) =>
  gunzipSync(readFileSync(join(DIR, `${gsm}.probes.gz`))).toString('utf8')

function rowsFromText(t: string) {
  const lines = t.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  if (h < 0) throw new Error('no header line found in the first 60 lines')
  const map = ingest.headerMap(lines[h])
  if (!map) throw new Error('the header carries none of the columns this reader needs')
  return function* () {
    for (let i = h + 1; i < lines.length; i += 1) {
      const r = ingest.parseRow(lines[i], map)
      if (r) yield r
    }
  }
}

/** Parents are cached and children are not: one array is about 90 MB live. */
const parentCache = new Map<string, ReturnType<typeof buildParent>>()
function buildParent(gsm: string) {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsFromText(textOf(gsm))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const profile = ingest.finishProfile(gsm, byChrom as never, bafSums as never, first)
  return score.finishParent(acc as never, profile.build.build)
}
const parentOf = (gsm: string) => {
  const hit = parentCache.get(gsm)
  if (hit) return hit
  const p = buildParent(gsm)
  parentCache.set(gsm, p)
  return p
}

async function scoreOne(
  childGsm: string, patGsm: string, role: 'paternal' | 'maternal', matGsm?: string,
) {
  const pat = parentOf(patGsm)
  const mat = matGsm ? parentOf(matGsm) : null
  const acc = score.emptyCollected(pat as never, mat as never)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsFromText(textOf(childGsm))()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, pat as never, mat as never, acc as never)
  }
  const profile = ingest.finishProfile(childGsm, byChrom as never, bafSums as never, first)
  return score.scoreSample({
    acc, profile, pat, mat: mat as never, soloRole: role, sibs: [],
    sampleName: childGsm, log: () => {},
  }) as Promise<Record<string, unknown>>
}

/** Every parent this run named, from either origin channel. What an operator would read. */
function namedParents(r: Record<string, unknown>): string[] {
  const out: string[] = []
  for (const c of (r.oneParent ?? []) as {
    where: string; parent?: string | null; twoParents?: boolean; verdict: string
  }[]) {
    if (c.twoParents && c.parent) out.push(`${c.where}:${c.parent}`)
    else if (c.verdict === 'known-parent-lost' || c.verdict === 'other-parent-lost') {
      out.push(`${c.where}:${c.verdict}`)
    }
  }
  for (const d of (r.dosageCalls ?? []) as { where: string; parent?: string | null }[]) {
    if (d.parent) out.push(`${d.where}:${d.parent}`)
  }
  return out
}

/**
 * The alerts an operator would actually see, from the shared list the three surfaces render.
 *
 * AN EARLIER VERSION OF THIS COUNTED ANY UNUSUAL FIELD and reported 8 of 8 runs "warned" in every
 * scenario INCLUDING the control, because embryo arrays are routinely non-diploid or complex. A
 * signal that fires on the correct run measures nothing. This counts only what the tool raises as
 * an operator-facing alert.
 */
const warnings = (r: Record<string, unknown>): { headline: string }[] =>
  (defects.runAlerts as (x: unknown) => { headline: string }[])(r)

let failures = 0
const fail = (what: string) => { failures += 1; console.log(`  *** FAIL  ${what}`) }

// ================================================================ ARM 1: WRONG FILE IN THE SLOT
console.log('=== ARM 1. THE WRONG ARRAY IN THE PARENT SLOT')
console.log('    Every file here is a valid, high-quality array of a real person. The only thing')
console.log('    wrong is which slot it went into, which is the mistake a hurried operator makes.')
{
  const kids = manifest.filter((r) => r.complete && r.role !== 'parent' && has(r.gsm)
    && has(r.father) && has(r.mother) && r.mother !== CONTAMINATED)
  const SAMPLE = kids.slice(0, Number(process.env.OM_N ?? 8))
  console.log(`    ${SAMPLE.length} children, each scored under every scenario below.`)
  console.log('')

  const tally: Record<string, { named: number; warned: number; runs: number }> = {}
  const note = (k: string, named: number, warned: boolean) => {
    tally[k] ??= { named: 0, warned: 0, runs: 0 }
    tally[k].named += named
    tally[k].warned += warned ? 1 : 0
    tally[k].runs += 1
  }

  for (const k of SAMPLE) {
    const sibling = kids.find((x) => x.gsm !== k.gsm && x.father === k.father
      && x.mother === k.mother)?.gsm
    const unrelated = EGG.find((g) => g !== k.mother) ?? EGG[0]
    type Scenario = [string, string, string, 'paternal' | 'maternal']
    const scenarios: Scenario[] = [
      // The control. Right array, right slot, declared correctly.
      ['control', k.gsm, k.father, 'paternal'],
      // THE ROLE SWAP. The true mother's array declared paternal. No genotype can object: the
      // arithmetic is correct and only the label is wrong, so every event comes out named for the
      // parent it did not come from.
      ['role-swap', k.gsm, k.mother, 'paternal'],
      // An unrelated adult, which must collapse.
      ['unrelated', k.gsm, unrelated, 'paternal'],
      // The sample as its own parent: one misclick in a file picker.
      ['self-as-parent', k.gsm, k.gsm, 'paternal'],
      // Generation inverted: the child in the parent slot, the parent scored as the sample.
      ['generation-inverted', k.father, k.gsm, 'paternal'],
    ]
    if (sibling) {
      // Within one IVF series every other embryo is a sibling, so this is the likeliest wrong file
      // of all: it is in the same folder, from the same case, with an adjacent name.
      scenarios.push(['sibling-as-parent', k.gsm, sibling, 'paternal'])
    }
    for (const [label, child, pgsm, role] of scenarios) {
      try {
        const r = await scoreOne(child, pgsm, role)
        note(label, namedParents(r).length, warnings(r).length > 0)
      } catch (e) {
        note(label, 0, true)
        console.log(`      ${label} on ${child}: threw ${String(e).slice(0, 70)}`)
      }
    }
  }

  console.log('    scenario             parents named   runs warned')
  for (const [k, v] of Object.entries(tally)) {
    console.log(`    ${k.padEnd(20)} ${String(v.named).padStart(10)}   ${v.warned}/${v.runs}`)
  }

  // THE ONE THAT MATTERS. A role swap that names parents without warning produces a report saying
  // "paternal" on every maternal event, at the confidence of a correct run.
  const rs = tally['role-swap']
  if (rs && rs.named > 0 && rs.warned < rs.runs) {
    fail(`role-swap named ${rs.named} parents and warned on only ${rs.warned} of ${rs.runs} runs. `
      + 'An egg donor in the paternal slot inverts every label in the report.')
  }
  const self = tally['self-as-parent']
  if (self && self.warned < self.runs) {
    fail(`self-as-parent warned on only ${self.warned} of ${self.runs} runs`)
  }
  const sib = tally['sibling-as-parent']
  if (sib && sib.named > 0 && sib.warned < sib.runs) {
    fail(`sibling-as-parent named ${sib.named} parents and warned on only ${sib.warned} of `
      + `${sib.runs} runs`)
  }
}

// ================================================================ ARM 2: THE FILE IS NOT AN ARRAY
console.log('')
console.log('=== ARM 2. MALFORMED AND TRUNCATED FILES')
console.log('    A refusal that names the problem is a PASS. A crash is a FAIL, because an operator')
console.log('    cannot tell one from a broken tool. Being ACCEPTED is the worst of the three.')
{
  const tmp = mkdtempSync(join(tmpdir(), 'om-stress-'))
  const raw = textOf(SPERM[0])
  const lines = raw.split('\n')
  let h = 0
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  // THE DELIMITER IS READ OFF THE FILE, NOT ASSUMED. An earlier version of this arm hard-coded a
  // tab; these exports are comma-separated, so two of its corruptions silently did nothing and the
  // unmodified file was reported as the tool accepting garbage. A harness that fabricates failures
  // is worse than no harness.
  const SEP = (lines[h].match(/\t/g)?.length ?? 0) > (lines[h].match(/,/g)?.length ?? 0)
    ? String.fromCharCode(9) : ','
  const gcol = lines[h].split(SEP).indexOf('genotype')
  if (gcol < 0) throw new Error('cannot find the genotype column to corrupt')
  const TAB = SEP

  const cases: [string, string][] = [
    ['empty', ''],
    ['header only', lines.slice(0, h + 1).join('\n')],
    ['no header', lines.slice(h + 1, h + 20_001).join('\n')],
    ['truncated 1%', lines.slice(0, h + 1 + Math.floor((lines.length - h) * 0.01)).join('\n')],
    ['truncated 50%', lines.slice(0, h + 1 + Math.floor((lines.length - h) * 0.5)).join('\n')],
    ['garbage bytes', 'not a tab separated file at all\n  binary junk\n'],
    ['header, junk rows', [lines[h], ...Array.from({ length: 20_000 },
      (_, i) => ['x' + i, 'ZZ', 'nope', '?', '?', '?', '??', '1'].join(TAB))].join('\n')],
    ['all no-call', [lines[h], ...lines.slice(h + 1, h + 100_001).map((l) => {
      const f = l.split(SEP)
      if (f.length > gcol) f[gcol] = '-1'
      return f.join(SEP)
    })].join('\n')],
    ['shuffled columns', [lines[h].split(SEP).reverse().join(SEP),
      ...lines.slice(h + 1, h + 100_001).map((l) => l.split(SEP).reverse().join(SEP))].join('\n')],
    ['one column', lines.slice(0, 20_000).map((l) => l.split(TAB)[0]).join('\n')],
  ]

  // WHAT THE WHOLE ARRAY SAYS, so a truncation can be checked against it rather than against zero.
  const wholeChild = manifest.find((x) => x.complete && x.role !== 'parent' && has(x.gsm))!
  const fullArrayCalls = namedParents(await scoreOne(wholeChild.gsm, SPERM[0], 'paternal'))

  console.log('')
  console.log(`    the whole array names ${fullArrayCalls.length} parent(s) on `
    + `${wholeChild.gsm}: ${fullArrayCalls.join(', ') || 'none'}`)
  console.log('    case                 outcome    detail')
  for (const [label, body] of cases) {
    const path = join(tmp, 'case.probes.gz')
    writeFileSync(path, gzipSync(Buffer.from(body, 'utf8')))
    let outcome = ''
    let detail = ''
    try {
      // THE WHOLE PIPELINE, NOT THE READER. Whether `finishParent` returns an object is not the
      // question: the question is whether a report comes out with a parent named on it. An earlier
      // version stopped at the reader and called a parsed row count "accepted", which measured the
      // parser rather than the tool.
      const acc = score.emptyParent()
      let n = 0
      for (const r of rowsFromText(gunzipSync(readFileSync(path)).toString('utf8'))()) {
        score.collectParentRow(r as never, acc as never); n += 1
      }
      const bad = score.finishParent(acc as never, 'GRCh38')
      const child = wholeChild
      const cacc = score.emptyCollected(bad as never, null)
      const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
      for (const r of rowsFromText(textOf(child.gsm))()) {
        if (!first) first = r.probesetId
        ingest.accumulate(r as never, byChrom as never)
        ingest.accumulateBaf(r as never, bafSums as never)
        score.collectRow(r as never, bad as never, null, cacc as never)
      }
      const profile = ingest.finishProfile(child.gsm, byChrom as never, bafSums as never, first)
      const res = await score.scoreSample({
        acc: cacc, profile, pat: bad, mat: null, soloRole: 'paternal', sibs: [],
        sampleName: child.gsm, log: () => {},
      }) as Record<string, unknown>
      const named = namedParents(res)
      outcome = named.length ? 'names a parent' : 'refused'
      detail = `${n} rows parsed, ${named.length} parents named`
      // A HALF FILE OF A REAL ARRAY IS STILL A REAL ARRAY, over the part it covers, so naming a
      // parent there is CORRECT and demanding silence would be demanding a worse tool. What must
      // hold is that it never invents one: the markers it lacks have no parental genotype at all,
      // which is a no-call rather than an absence, so its answers have to be a SUBSET of what the
      // whole array gives. Verified: the half file covers chr1 to chr8, and on a real blastomere it
      // returns the same two chr6 calls as the full array, identical.
      if (label === 'truncated 50%') {
        const extra = named.filter((x) => !fullArrayCalls.includes(x))
        if (extra.length) {
          fail(`${label}: invented ${extra.length} call(s) the whole array does not make: `
            + extra.join(', '))
        }
        detail += `, all of them also made by the whole array`
      } else if (named.length) {
        fail(`${label}: named ${named.length} parents off a parental array that is not one`)
      }
    } catch (e) {
      const m = String((e as Error).message ?? e)
      // A refusal has to explain itself. "undefined is not a function" is a crash.
      const legible = !/undefined|not a function|Cannot read|null/i.test(m)
      outcome = legible ? 'refused' : 'CRASH'
      detail = m.slice(0, 66)
      if (!legible) fail(`${label}: crashed rather than refusing: ${m.slice(0, 60)}`)
    }
    console.log(`    ${label.padEnd(20)} ${outcome.padEnd(10)} ${detail}`)
  }
  rmSync(tmp, { recursive: true, force: true })
}

// ================================================================ ARM 3: THE SAME ANSWER TWICE
console.log('')
console.log('=== ARM 3. DETERMINISM')
console.log('    A tool that answers differently on a re-run cannot be checked by its user.')
{
  const k = manifest.find((r) => r.complete && r.role !== 'parent' && has(r.gsm) && has(r.father)
    && has(r.mother) && r.mother !== CONTAMINATED)!
  const a = await scoreOne(k.gsm, k.father, 'paternal', k.mother)
  const b = await scoreOne(k.gsm, k.father, 'paternal', k.mother)
  const sa = JSON.stringify(a)
  const sb = JSON.stringify(b)
  console.log(`    ${k.gsm} scored twice: ${sa === sb ? 'IDENTICAL' : 'DIFFERENT'}`
    + ` (${sa.length} bytes)`)
  if (sa !== sb) fail('two runs of one input disagreed')
}

// ================================================================ ARM 4: EVERY ARRAY, NO CRASH
console.log('')
console.log('=== ARM 4. EVERY ARRAY IN THE CORPUS THROUGH THE FULL PIPELINE')
{
  const all = manifest.filter((r) => has(r.gsm)).map((r) => r.gsm)
  const ref = SPERM[0]
  const limit = Number(process.env.OM_ALL_N ?? all.length)
  let ok = 0
  let threw = 0
  const errs: string[] = []
  for (const g of all.slice(0, limit)) {
    try { await scoreOne(g, ref, 'paternal'); ok += 1 } catch (e) {
      threw += 1
      errs.push(`${g}: ${String((e as Error).message ?? e).slice(0, 66)}`)
    }
  }
  console.log(`    scored ${ok}, threw ${threw}`)
  for (const e of errs.slice(0, 10)) console.log(`      ${e}`)
  if (threw) fail(`${threw} arrays threw in the full pipeline`)
}

console.log('')
console.log(failures === 0
  ? 'stress: every scenario either answered correctly or refused visibly.'
  : `stress: ${failures} FAILURES. Each is a case where the tool answered without being right.`)
process.exit(failures ? 1 : 0)
