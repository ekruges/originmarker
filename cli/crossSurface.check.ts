// THE TWO SURFACES MUST GIVE ONE ANSWER. Run: node cli/crossSurface.check.ts
//
// The browser and the command line each answer "which parent's copy is this". They are different
// programs reading the same file, and twice they answered differently about the same losses on the
// same array: once because the command line assembled its tally without the B-allele frequency,
// which moved the sample across the diploid heterozygosity boundary, and once because a refusal
// verdict was renamed and only one surface's gate still matched it.
//
// Neither was reachable by any check. Every other check in this repo tests a unit; the disagreement
// lived in how each surface ASSEMBLED those units, which is exactly the seam a unit test cannot
// see. So this one runs a whole array through both surfaces and compares the answers.
//
// It runs on a synthetic array built here, so it runs in CI with no data to supply. Set OM_ARRAYS
// to a directory holding `parent.probes` and `sample.probes` and it additionally pins the answer on
// real material; when that is unset it says so rather than passing quietly.
import assert from 'node:assert/strict'
import { execFileSync } from 'node:child_process'
import { mkdtempSync, writeFileSync, readdirSync, readFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const defects = await import(`${W}defects.ts`)

// ---------------------------------------------------------------- a synthetic array
//
// Deliberately simple and deliberately DECIDABLE. The sample carries one parental genome and is
// missing two whole chromosomes, which is the shape that exercises the path the two surfaces
// disagreed on: no dosage channel can reach a whole-chromosome loss on this material, so the
// answer has to come from the genome-level zygosity call, inherited by every event.

const CHROMS = Array.from({ length: 22 }, (_, i) => String(i + 1))
// A chromosome needs 200 informative markers to be reported at all, and on a lost one only the
// surviving quarter of calls are informative, so the per-chromosome count has to clear that floor
// AFTER the collapse rather than before it.
const PER_CHROM = 2400
const LOST = new Set(['13', '19'])

// THE FIXTURE IS BUILT TO THE NUMBERS A REAL RUN PRODUCED, and that is the point of these three.
//
// Zygosity is read from the fraction of B-allele frequencies in the heterozygous band, and falls
// back to genotype heterozygosity against half the parent's own rate only where no frequencies
// exist. Those two measures do not have to agree, and on the array this fixture is modelled on
// they landed on OPPOSITE SIDES of the 0.08 boundary: 0.058 by band, 0.090 by genotype, against a
// parent measuring 0.170 and therefore a fallback threshold of 0.085.
//
// A fixture built on round numbers agrees with itself either way and cannot tell the two measures
// apart. This one straddles the boundary exactly as the real array does, so discarding the
// frequencies flips the genome-level call, which is what it did in production.
const PARENT_HET = 0.170
const SAMPLE_HET = 0.090
/** Share of heterozygous CALLS whose intensity is not balanced: 0.090 by genotype, 0.058 by band. */
const BAND_EXTREME = 0.355

/** A deterministic stream, so a failure is the code changing and never the draw changing. */
function rng(seed: number): () => number {
  let x = seed >>> 0
  return () => {
    x ^= x << 13; x >>>= 0
    x ^= x >> 17
    x ^= x << 5; x >>>= 0
    return x / 0x100000000
  }
}

const HEAD = 'probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset'

/**
 * One parent and one gynogenetic sample.
 *
 * The parent is an ordinary diploid. The sample carries none of it: at every marker where the
 * parent is homozygous, the sample reads the allele the parent does not have, which is what a
 * genome with no paternal contribution looks like to this tool. Two chromosomes are absent
 * outright, with their calls dropped and their intensity down, which is how a whole-chromosome
 * loss presents.
 */
function synthesise(): { parent: string; sample: string } {
  const r = rng(20260820)
  const parent = [HEAD]
  const sample = [HEAD]
  for (const chrom of CHROMS) {
    for (let i = 0; i < PER_CHROM; i += 1) {
      const id = `AX-${chrom}-${i}`
      const pos = 100_000 + i * 90_000
      // 17% heterozygous, which is what the loaded parent of the real run measures.
      const pg = r() < PARENT_HET ? 1 : (r() < 0.5 ? 0 : 2)
      parent.push([id, chrom, pos, (r() * 0.1 - 0.05).toFixed(4),
        pg === 1 ? (0.5 + (r() - 0.5) * 0.06).toFixed(4) : (pg === 0 ? 0.02 : 0.98).toFixed(4),
        '2.0', String(pg), '1'].join('\t'))

      const lost = LOST.has(chrom)
      // Homozygous end to end, and at the OPPOSITE allele wherever the parent is homozygous.
      // Some markers read heterozygous anyway: that is amplification error and every real array
      // carries it. The rate is set where a real one sat, for the reason under BAND_EXTREME.
      const noise = r() < SAMPLE_HET
      const sg = noise ? 1 : pg === 0 ? 2 : pg === 2 ? 0 : (r() < 0.5 ? 0 : 2)
      // A heterozygous CALL whose intensity is nowhere near balanced. This is the difference
      // between the two measures of zygosity, and it is the whole reason the fixture is built to
      // these numbers rather than round ones.
      const bandExtreme = noise && r() < BAND_EXTREME
      const called = lost ? r() < 0.25 : true
      sample.push([id, chrom, pos,
        (lost ? -0.85 : 0) + (r() * 0.1 - 0.05),
        sg === 1 && !bandExtreme ? 0.5 + (r() - 0.5) * 0.06
          : sg === 1 ? (r() < 0.5 ? 0.03 + r() * 0.02 : 0.95 + r() * 0.02)
            : sg === 0 ? 0.02 + r() * 0.02 : 0.96 + r() * 0.02,
        lost ? '1.0' : '2.0',
        called ? String(sg) : '-1', '1'].join('\t'))
    }
  }
  return { parent: parent.join('\n'), sample: sample.join('\n') }
}

/**
 * One parent and one BIPARENTAL sample missing the loaded parent's copy of two chromosomes.
 *
 * The other fixture is uniparental, where the parent is not in question and every event inherits
 * one genome-level answer. This is the case a lab actually brings: a cell with both parents in it,
 * where naming the parent takes Mendelian evidence. At a marker where the loaded parent is
 * homozygous AA, losing that parent's copy leaves the sample reading BB, an allele the parent does
 * not have. Dropout removes alleles and never invents one, so BB there means the parent's copy is
 * gone. That is what `callOneParentOrigin` reads, and it is the only fixture here that exercises
 * the verdict-to-parent mapping: the uniparental path names its parent without consulting it.
 */
function synthesiseBiparental(): { parent: string; sample: string } {
  const r = rng(20260821)
  const parent = [HEAD]
  const sample = [HEAD]
  for (const chrom of CHROMS) {
    for (let i = 0; i < PER_CHROM; i += 1) {
      const id = `AX-${chrom}-${i}`
      const pos = 100_000 + i * 90_000
      const pg = r() < PARENT_HET ? 1 : (r() < 0.5 ? 0 : 2)
      parent.push([id, chrom, pos, (r() * 0.1 - 0.05).toFixed(4),
        pg === 1 ? (0.5 + (r() - 0.5) * 0.06).toFixed(4) : (pg === 0 ? 0.02 : 0.98).toFixed(4),
        '2.0', String(pg), '1'].join('\t'))

      const lost = LOST.has(chrom)
      // What the loaded parent transmitted, and what the other parent transmitted alongside it.
      const fromLoaded = pg === 1 ? (r() < 0.5 ? 0 : 1) : (pg === 0 ? 0 : 1)
      const fromOther = r() < 0.5 ? 0 : 1
      // On a lost chromosome the loaded parent's allele is not there at all, so the sample is
      // homozygous for whatever the other parent gave: at a parent-AA marker that reads BB.
      const alleles = lost ? [fromOther, fromOther] : [fromLoaded, fromOther]
      const sg = alleles[0] + alleles[1]
      const called = lost ? r() < 0.55 : r() < 0.98
      sample.push([id, chrom, pos,
        (lost ? -0.60 : 0) + (r() * 0.1 - 0.05),
        sg === 1 ? 0.5 + (r() - 0.5) * 0.08 : sg === 0 ? 0.02 + r() * 0.02 : 0.96 + r() * 0.02,
        lost ? '1.0' : '2.0',
        called ? String(sg) : '-1', '1'].join('\t'))
    }
  }
  return { parent: parent.join('\n'), sample: sample.join('\n') }
}

// ---------------------------------------------------------------- surface A: the browser's path
//
// This is what SyngamyPage does per sample, with the File reader replaced by a string split. Every
// decision below it is scoreSample's, which is the point: there is nothing here to keep in step.

function rowsOf(text: string) {
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  const map = ingest.headerMap(lines[h])
  assert.ok(map, 'the synthetic header must parse through the shipped ingest')
  const out = []
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = ingest.parseRow(lines[i], map)
    if (r) out.push(r)
  }
  return out
}

async function browserAnswer(parentText: string, sampleText: string, role = 'paternal') {
  const pacc = score.emptyParent()
  let byChrom = new Map(); let bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rowsOf(parentText)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom)
    ingest.accumulateBaf(r, bafSums)
    score.collectParentRow(r, pacc)
  }
  const pProfile = ingest.finishProfile('parent', byChrom, bafSums, first)
  const pat = score.finishParent(pacc, pProfile.build.build)

  const acc = score.emptyCollected(pat, null)
  byChrom = new Map(); bafSums = ingest.emptyBafSums(); first = ''
  for (const r of rowsOf(sampleText)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom)
    ingest.accumulateBaf(r, bafSums)
    score.collectRow(r, pat, null, acc)
  }
  const profile = ingest.finishProfile('sample', byChrom, bafSums, first)
  const result = await score.scoreSample({
    acc, profile, pat, mat: null, soloRole: role, sibs: [], sampleName: 'sample',
    log: () => {},
  })
  return result
}

/** The rows a reader compares, from either surface, in one shape. */
const comparable = (result: {
  originClass?: string; zygosity?: string
  dosageCalls?: unknown[]; oneParent?: unknown[]
}, role: 'paternal' | 'maternal') => {
  const rows = new Map<string, { verdict: string; origin: string | null; band?: string }>()
  for (const c of (result.dosageCalls ?? [])) {
    const d = c as { where: string; verdict: string; band?: string; parent?: string }
    rows.set(d.where, {
      verdict: d.verdict,
      origin: d.parent ?? defects.parentNamed(d.verdict, role),
      band: d.band,
    })
  }
  for (const c of (result.oneParent ?? [])) {
    const g = c as { where: string; verdict: string; band?: string }
    const named = defects.parentNamed(g.verdict, role)
    if (!named && rows.has(g.where)) continue
    rows.set(g.where, { verdict: g.verdict, origin: named, band: g.band })
  }
  return {
    originClass: result.originClass,
    zygosity: result.zygosity,
    rows: [...rows.entries()].sort(([a], [b]) => a.localeCompare(b)),
  }
}

// ---------------------------------------------------------------- surface B: the command line

function cliAnswer(dir: string, role = 'paternal') {
  const outText = execFileSync(process.execPath, [
    '--experimental-strip-types',
    new URL('./om.ts', import.meta.url).pathname,
    'origin', join(dir, 'parent.probes'), join(dir, 'sample.probes'),
    '--role', role, '--json',
  ], { encoding: 'utf8', maxBuffer: 64 * 1024 * 1024 })
  return JSON.parse(outText) as {
    originClass?: string; zygosity?: string
    events: { locus: string; verdict: string; origin: string | null; band?: string }[]
  }
}

const cliComparable = (a: ReturnType<typeof cliAnswer>) => ({
  originClass: a.originClass,
  zygosity: a.zygosity,
  rows: a.events.map((e) => [e.locus, { verdict: e.verdict, origin: e.origin, band: e.band }])
    .sort((x, y) => String(x[0]).localeCompare(String(y[0]))),
})

// ---------------------------------------------------------------- 1. the synthetic array
{
  const { parent, sample } = synthesise()
  const dir = mkdtempSync(join(tmpdir(), 'om-cross-'))
  writeFileSync(join(dir, 'parent.probes'), parent)
  writeFileSync(join(dir, 'sample.probes'), sample)

  const web = comparable(await browserAnswer(parent, sample) as never, 'paternal')
  const cli = cliComparable(cliAnswer(dir))

  // FIRST, IS THE ANSWER RIGHT. Agreement alone is not correctness: both surfaces now share every
  // decision, so a defect in the shared code moves them together and this comparison would still
  // pass. What is asserted here is the truth of the fixture, which is known because it was built:
  // one parental genome, and exactly the two chromosomes that were removed.
  assert.equal(web.originClass, 'gynogenetic',
    `the sample carries one parental genome and none of the loaded parent's, so it must read `
    + `gynogenetic. Got ${web.originClass}.`)
  assert.equal(web.zygosity, 'uniparental_homozygous',
    `a genome with one parental contribution is homozygous end to end. Got ${web.zygosity}.`)
  assert.deepEqual(web.rows.map(([where]) => where).sort(), ['chr13', 'chr19'],
    `exactly the two removed chromosomes must be reported, no more and no fewer. Got `
    + `${JSON.stringify(web.rows.map(([w]) => w))}. Extra rows are noise being read as biology; `
    + 'missing rows are events being lost.')
  for (const [where, row] of web.rows) {
    assert.equal(row.origin, 'maternal',
      `${where} was removed from a genome with no paternal copy in it, so the lost copy was `
      + `maternal by construction. Got ${row.origin}.`)
  }

  assert.equal(cli.originClass, web.originClass, 'the genome-level call must match')
  assert.equal(cli.zygosity, web.zygosity, 'the zygosity must match')
  assert.deepEqual(cli.rows, web.rows,
    'THE TWO SURFACES DISAGREE ABOUT THE SAME FILE. This is the failure this check exists for: '
    + 'one of them is reading rows differently, mapping verdicts differently, or has grown a '
    + 'second implementation of a channel.')

  const named = web.rows.filter(([, v]) => v.origin).length
  console.log(`  synthetic array: ${web.rows.length} event(s), ${named} with a named parent, `
    + `identical on both surfaces (${web.originClass}, ${web.zygosity})`)
}

// ---------------------------------------------------------------- 1b. a biparental array
//
// The case a lab brings. Both parents are in the cell, so the parent is genuinely in question and
// has to be MEASURED rather than inherited from a genome-level call.
{
  const { parent, sample } = synthesiseBiparental()
  const dir = mkdtempSync(join(tmpdir(), 'om-cross-bi-'))
  writeFileSync(join(dir, 'parent.probes'), parent)
  writeFileSync(join(dir, 'sample.probes'), sample)

  const web = comparable(await browserAnswer(parent, sample) as never, 'paternal')
  const cli = cliComparable(cliAnswer(dir))

  assert.equal(web.zygosity, 'diploid',
    `a cell carrying both parents is heterozygous where they differ, so it must read diploid. `
    + `Got ${web.zygosity}. If this fails the fixture stopped being biparental, and the mapping `
    + 'below is no longer being exercised.')
  assert.notEqual(web.originClass, 'gynogenetic',
    'the loaded parent IS present here, so the genome-level absence call must not fire')

  const named = web.rows.filter(([, v]) => v.origin)
  assert.ok(named.length >= 1,
    `the loaded parent's copy was removed from ${[...LOST].join(' and ')}, which is Mendelian `
    + 'evidence needing no detection floor, so at least one event must name a parent. Got '
    + `${web.rows.length} event(s), ${named.length} named.`)
  for (const [where, row] of named) {
    assert.equal(row.origin, 'paternal',
      `${where} is missing the LOADED parent's copy, and the loaded parent here is paternal, so `
      + `the lost copy is paternal. Got ${row.origin}. An inverted verdict-to-parent mapping `
      + 'reads exactly like this.')
  }

  assert.equal(cli.zygosity, web.zygosity, 'the zygosity must match')
  assert.deepEqual(cli.rows, web.rows,
    'THE TWO SURFACES DISAGREE ABOUT THE SAME BIPARENTAL FILE.')
  console.log(`  biparental array: ${web.rows.length} event(s), ${named.length} naming a `
    + `measured parent, identical on both surfaces (${web.zygosity})`)
}

// ---------------------------------------------------------------- 2. no third implementation
//
// A shared module only helps while it is the only caller. These are the channels that decide a
// parent; anything outside scoreSample.ts calling one directly is a second answer being born.
{
  // Each channel, and the one module allowed to contain it: its own, which declares it, and
  // scoreSample.ts, which calls it.
  const CHANNELS: [string, string][] = [
    ['callDosageOrigin', 'dosageOrigin.ts'],
    ['uniparentalOrigin', 'uniparentalOrigin.ts'],
    ['callOneParentOrigin', 'oneParentOrigin.ts'],
  ]
  const src = new URL('../web/src/', import.meta.url).pathname
  const offenders: string[] = []
  for (const f of readdirSync(src)) {
    if (!/\.tsx?$/.test(f) || f.endsWith('.check.ts')) continue
    if (f === 'scoreSample.ts') continue
    const text = readFileSync(join(src, f), 'utf8')
    for (const [c, home] of CHANNELS) {
      if (f === home) continue
      if (new RegExp(`\\b${c}\\s*\\(`).test(text)) offenders.push(`${f} calls ${c}`)
    }
  }
  assert.deepEqual(offenders, [],
    'a parental-origin channel is being called outside scoreSample.ts, which is how the two '
    + 'surfaces came to disagree in the first place')
  console.log(`  ${CHANNELS.length} origin channels, called from scoreSample.ts and nowhere else`)
}

// ---------------------------------------------------------------- 3. real material, when supplied
{
  const dir = process.env.OM_ARRAYS
  if (!dir) {
    console.log('  real arrays: NOT CHECKED. Set OM_ARRAYS to a directory holding parent.probes '
      + 'and sample.probes to pin the answer on real material.')
  } else {
    const web = comparable(await browserAnswer(
      readFileSync(join(dir, 'parent.probes'), 'utf8'),
      readFileSync(join(dir, 'sample.probes'), 'utf8'),
    ) as never, 'paternal')
    const cli = cliComparable(cliAnswer(dir))
    assert.equal(cli.originClass, web.originClass)
    assert.equal(cli.zygosity, web.zygosity)
    assert.deepEqual(cli.rows, web.rows, 'the two surfaces disagree on real material')
    console.log(`  real arrays: ${web.rows.length} event(s), identical on both surfaces `
      + `(${web.originClass}, ${web.zygosity})`)
  }
}

console.log('crossSurface.check.ts: one answer, whichever surface asks')
