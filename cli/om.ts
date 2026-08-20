#!/usr/bin/env -S node --experimental-strip-types
/**
 * OriginMarker on the command line.
 *
 * The browser tool is deliberately a mindless drag and drop: it infers what it can and refuses
 * what it cannot, and it does not offer you a knob whose correct setting you would have to know.
 * This is the other end of that. Every constant the modules use is exposed here as a flag, the
 * intermediate quantities are printable, and the batch commands take a directory rather than a
 * gesture. It is for the person who already knows what they are changing and why.
 *
 * IT IS THE SAME CODE. Every command imports the module the web tool runs, so nothing here can
 * drift from what the app does or skip a guard the app applies. The one thing this adds is the
 * ability to move a threshold, and every command that has been moved off its shipped default
 * says so in its own output, because a number produced under a changed constant is not
 * comparable with one produced under the shipped configuration.
 *
 *   om stage      <array>                  what material an array is, and the dropout it implies
 *   om origin     <parent> <sample>        which parent's copy is missing, per region
 *   om link       <parent> <sample>...     parent, duplicate, unrelated, or refused
 *   om cohort     <dir> --ref <array>      origin for every confirmed child in a directory
 *   om census     <dir>                    haploid products per group, for reconstruction
 *   om reconstruct <product>...            a parent's genotypes from that parent's haploid cells
 *   om enrich     <regions.tsv> --track    positional enrichment, marker-matched null
 *   om constants                           every tunable, its shipped value and its provenance
 *
 * Run any command with --help for its own flags, or --json for machine-readable output.
 */
import { readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs'
import { join, basename } from 'node:path'
import { gunzipSync } from 'node:zlib'

const W = new URL('../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const parentage = await import(`${W}parentage.ts`)
const stageMod = await import(`${W}stage.ts`)
const oneParent = await import(`${W}oneParentOrigin.ts`)
const uniparental = await import(`${W}uniparentalOrigin.ts`)
const obligate = await import(`${W}obligateHet.ts`)
const segments = await import(`${W}segments.ts`)
const inferredRef = await import(`${W}inferredReference.ts`)
const features = await import(`${W}features.ts`)
const dosage = await import(`${W}dosageOrigin.ts`)
const post = await import(`${W}originPosterior.ts`)
const tax = await import(`${W}abnormalities.ts`)
const score = await import(`${W}scoreSample.ts`)
const defects = await import(`${W}defects.ts`)

type AB = 'AA' | 'AB' | 'BB' | 'NC'

// ------------------------------------------------------------------ argument parsing

interface Args {
  cmd: string
  positional: string[]
  flags: Map<string, string>
  bools: Set<string>
}

function parseArgs(argv: string[]): Args {
  const flags = new Map<string, string>()
  const bools = new Set<string>()
  const positional: string[] = []
  const cmd = argv[0] && !argv[0].startsWith('-') ? argv[0] : 'help'
  for (let i = cmd === 'help' ? 0 : 1; i < argv.length; i += 1) {
    const a = argv[i]
    if (!a.startsWith('--')) { positional.push(a); continue }
    const eq = a.indexOf('=')
    if (eq > 0) { flags.set(a.slice(2, eq), a.slice(eq + 1)); continue }
    const next = argv[i + 1]
    if (next && !next.startsWith('--')) { flags.set(a.slice(2), next); i += 1 } else bools.add(a.slice(2))
  }
  return { cmd, positional, flags, bools }
}

const args = parseArgs(process.argv.slice(2))
const JSON_OUT = args.bools.has('json')

/** A numeric flag, with the shipped default. Every override is recorded so output can say so. */
const overrides: { name: string, from: number, to: number }[] = []
function num(name: string, shipped: number): number {
  const raw = args.flags.get(name)
  if (raw === undefined) return shipped
  const v = Number(raw)
  if (!Number.isFinite(v)) die(`--${name} needs a number, got ${raw}`)
  if (v !== shipped) overrides.push({ name, from: shipped, to: v })
  return v
}

function die(msg: string): never {
  process.stderr.write(`om: ${msg}\n`)
  process.exit(2)
}

function out(obj: unknown, text: () => void): void {
  if (JSON_OUT) {
    process.stdout.write(`${JSON.stringify(
      overrides.length ? { ...(obj as object), tuning: overrides } : obj, null, 1,
    )}\n`)
  } else {
    text()
    if (overrides.length) {
      process.stdout.write('\nNOT THE SHIPPED CONFIGURATION. '
        + `${overrides.map((o) => `${o.name} ${o.from} -> ${o.to}`).join(', ')}. `
        + 'Numbers produced under a changed constant are not comparable with numbers produced '
        + 'under the defaults, and the validation figures in the audit do not apply to them.\n')
    }
  }
}

// ------------------------------------------------------------------ loading

interface Loaded {
  id: string
  gt: Map<string, AB>
  pos: Map<string, { chrom: string, pos: number }>
  /** B-allele frequency per marker. Read whether or not a genotype was emitted, which is the
   *  whole reason the dosage channel survives a collapsed call rate. */
  baf: Map<string, number>
  cnByChrom: Map<string, { chrom: string, pos: number, called: boolean, log2R: number | null }[]>
  probeAt: Map<string, string>
  rows: number
  profile: ReturnType<typeof ingest.finishProfile>
}

/**
 * Read one array through the shipped ingest.
 *
 * `stride` subsamples for screening only. It is exposed because a linkage rate over 100,000
 * markers has a standard error under 0.001 against boundaries of 0.02 and 0.105, so a full read
 * buys nothing there and costs eight times as much. Anything reported, rather than screened, is
 * read at stride 1.
 */
function load(path: string, stride = 1, full = true): Loaded {
  const raw = readFileSync(path)
  const text = (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  if (h < 0) die(`${path}: no header found in the first 60 lines`)
  const map = ingest.headerMap(lines[h])
  if (!map) die(`${path}: header not recognised. Columns were: ${lines[h].slice(0, 120)}`)

  const gt = new Map<string, AB>()
  const pos = new Map<string, { chrom: string, pos: number }>()
  const baf = new Map<string, number>()
  const cnByChrom = new Map<string, Loaded['cnByChrom'] extends Map<string, infer V> ? V : never>()
  const probeAt = new Map<string, string>()
  const byChrom = new Map()
  const bafSums = ingest.emptyBafSums()
  let first = ''
  let rows = 0

  for (let i = h + 1; i < lines.length; i += stride) {
    const r = ingest.parseRow(lines[i], map)
    if (!r) continue
    rows += 1
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    if (!parentage.isAutosome(r.chrom)) continue
    if (r.genotype !== 'NC') gt.set(r.probesetId, r.genotype as AB)
    if (!full) continue
    // Position and dosage are kept for EVERY marker, called or not. A no-call still has an
    // intensity reading, and on a collapsed chromosome those are the only readings left.
    pos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
    if (r.baf !== null && Number.isFinite(r.baf)) baf.set(r.probesetId, r.baf)
    const cn = cnByChrom.get(r.chrom) ?? []
    cn.push({ chrom: r.chrom, pos: r.pos, called: r.genotype !== 'NC', log2R: r.log2R })
    cnByChrom.set(r.chrom, cn as never)
    probeAt.set(`${r.chrom}:${r.pos}`, r.probesetId)
  }
  const id = basename(path).replace(/_\d+\.CEL\.probes$/, '')
    .replace(/\.(probes|CEL\.txt\.gz|CEL\.txt|csv\.gz|txt\.gz)$/, '')
    .replace(/^(GSM\d+)_.*$/, '$1')
  return {
    id, gt, pos, baf, cnByChrom, probeAt, rows,
    profile: ingest.finishProfile(id, byChrom as never, bafSums as never, first),
  }
}

function arraysUnder(dir: string, rel = ''): string[] {
  const outp: string[] = []
  for (const e of readdirSync(join(dir, rel)).sort()) {
    const r = rel ? join(rel, e) : e
    if (statSync(join(dir, r)).isDirectory()) outp.push(...arraysUnder(dir, r))
    else if (/\.(probes|probes\.gz|CEL\.txt\.gz|txt\.gz)$/.test(e)) outp.push(r)
  }
  return outp
}

const pct = (x: number, d = 1): string => (Number.isFinite(x) ? `${(100 * x).toFixed(d)}%` : 'n/a')

// ------------------------------------------------------------------ shared pieces

/** Opposite homozygotes. A parent and child cannot be AA and BB at the same marker. */
function oppositeHom(a: Map<string, AB>, b: Map<string, AB>): { rate: number, n: number } {
  let n = 0
  let opp = 0
  for (const [probe, ga] of a) {
    if (ga !== 'AA' && ga !== 'BB') continue
    const gb = b.get(probe)
    if (gb !== 'AA' && gb !== 'BB') continue
    n += 1
    if (ga !== gb) opp += 1
  }
  return { rate: n ? opp / n : NaN, n }
}

/** The one-parent heterozygosity call: is there a second parental contribution, or not. */
function secondParent(ref: Map<string, AB>, s: Map<string, AB>) {
  const t = obligate.emptyHet()
  for (const [probe, pg] of ref) {
    const cg = s.get(probe)
    if (cg) obligate.addOneParent(pg as never, cg as never, t as never)
  }
  return obligate.hetCall(t as never, 1) as {
    ploidy: string, fraction: number, informative: number, why: string
  }
}

/** Which parent a verdict names. The module's exhaustive mapping, never a local re-spelling. */
const parentNamed = (verdict: string, loaded: 'paternal' | 'maternal'): string | null =>
  (defects.parentNamed as (v: string, l: string) => string | null)(verdict, loaded)

/**
 * Stage thresholds as flags.
 *
 * Every constant the modules use is a flag here, and a number produced under a changed constant is
 * not comparable with one produced under the shipped configuration, so every command that has been
 * moved off a default says so in its own output.
 */
const stageOpts = () => ({
  callFloor: num('qc-call-floor', stageMod.QC_CALL_FLOOR),
  maxDiploidHet: num('max-diploid-het', stageMod.MAX_DIPLOID_HET),
  haploidMaxHet: num('haploid-max-het', stageMod.HAPLOID_MAX_HET),
  bulkHeterozygosity: num('bulk-heterozygosity', stageMod.BULK_HETEROZYGOSITY),
  bandDiploidCertain: num('band-diploid-certain', stageMod.BAND_DIPLOID_CERTAIN),
  bulkMaxBafSd: num('bulk-max-baf-sd', stageMod.BULK_MAX_BAF_SD),
})


/** One array's marker rows, as the shipped ingest parses them. */
function readRows(path: string) {
  const raw = readFileSync(path)
  const text = (path.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  if (h < 0) die(`${path}: no header found in the first 60 lines`)
  const map = ingest.headerMap(lines[h])
  if (!map) die(`${path}: header not recognised. Columns were: ${lines[h].slice(0, 120)}`)
  return function* () {
    for (let i = h + 1; i < lines.length; i += 1) {
      const r = ingest.parseRow(lines[i], map)
      if (r) yield r
    }
  }
}

const idOf = (path: string) => basename(path).replace(/_\d+\.CEL\.probes$/, '')
  .replace(/\.(probes|CEL\.txt\.gz|CEL\.txt|csv\.gz|txt\.gz)$/, '')
  .replace(/^(GSM\d+)_.*$/, '$1')

/** The parent, through the accumulator the browser fills. Read once and reused across a cohort. */
function loadParent(path: string) {
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of readRows(path)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectParentRow(r as never, acc as never)
  }
  const profile = ingest.finishProfile(idOf(path), byChrom as never, bafSums as never, first)
  return { id: idOf(path), profile, pat: score.finishParent(acc as never, profile.build.build) }
}

/**
 * One sample scored, through the module the browser runs.
 *
 * THERE IS NO SECOND IMPLEMENTATION HERE, and there must not be one. Every channel, guard and
 * refusal comes from scoreSample, so this surface cannot answer differently from the app on the
 * same file. Flags reach the science through scoreSample's own inputs, never through a fork of it.
 */
async function scoreWithParent(
  parent: ReturnType<typeof loadParent>, samplePath: string, role: 'paternal' | 'maternal',
  onLog: (tag: string, text: string) => void = () => {},
) {
  const acc = score.emptyCollected(parent.pat as never, null)
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of readRows(samplePath)()) {
    if (!first) first = r.probesetId
    ingest.accumulate(r as never, byChrom as never)
    ingest.accumulateBaf(r as never, bafSums as never)
    score.collectRow(r as never, parent.pat as never, null, acc as never)
  }
  const profile = ingest.finishProfile(idOf(samplePath), byChrom as never, bafSums as never, first)
  const result = await score.scoreSample({
    acc, profile, pat: parent.pat, mat: null, soloRole: role, sibs: [],
    sampleName: idOf(samplePath), log: onLog, stageOpts: stageOpts(),
  })
  return { id: idOf(samplePath), refId: parent.id, profile, result }
}

interface OriginRow {
  locus: string; kind: string; channel: string; origin: string | null
  verdict: string; band?: string; confidence?: number; markers?: number; why: string
}

/**
 * One row per event, from the channels scoreSample already ran.
 *
 * Reads rather than re-derives. The Mendelian channel needs no detection floor, so where it named
 * a parent it outranks the dosage channel on the same interval; where it refused, the dosage
 * channel's answer for that interval stands.
 */
function originRows(
  result: { dosageCalls?: unknown[], oneParent?: unknown[] }, role: 'paternal' | 'maternal',
): OriginRow[] {
  const byWhere = new Map<string, OriginRow>()
  for (const c of (result.dosageCalls ?? [])) {
    const d = c as {
      where: string, verdict: string, band?: string, confidence?: number, markers?: number,
      why: string, cls?: string, parent?: string, fromZygosity?: boolean,
    }
    byWhere.set(d.where, {
      locus: d.where,
      kind: d.cls ?? 'dosage',
      channel: d.fromZygosity ? 'zygosity' : 'dosage',
      origin: d.parent ?? parentNamed(d.verdict, role),
      verdict: d.verdict, band: d.band, confidence: d.confidence, markers: d.markers, why: d.why,
    })
  }
  for (const c of (result.oneParent ?? [])) {
    const g = c as {
      where: string, verdict: string, band?: string, posterior?: number, markers?: number,
      why: string,
    }
    const named = parentNamed(g.verdict, role)
    const prior = byWhere.get(g.where)
    if (!named && prior) continue
    byWhere.set(g.where, {
      locus: g.where, kind: prior?.kind ?? 'genotype', channel: 'genotype',
      origin: named ?? prior?.origin ?? null,
      verdict: g.verdict, band: g.band, confidence: g.posterior, markers: g.markers, why: g.why,
    })
  }
  return [...byWhere.values()]
}

const scoreArray = async (
  refPath: string, samplePath: string, role: 'paternal' | 'maternal',
  onLog?: (tag: string, text: string) => void,
) => scoreWithParent(loadParent(refPath), samplePath, role, onLog)


// ------------------------------------------------------------------ commands

const COMMANDS: Record<string, () => void | Promise<void>> = {
  stage() {
    const p = args.positional[0]
    if (!p) die('usage: om stage <array> [--json]')
    const a = load(p, num('stride', 1), false)
    const s = stageMod.inferStage({
      hetRate: a.profile.hetRate,
      callRate: a.profile.callRate,
      hetBand: args.bools.has('no-band') ? undefined : a.profile.hetBand,
      // The amplification axis. Omitting it here would have left this one command on the old
      // one-axis ladder while every other path used both, which is the kind of divergence a
      // literal built field by field invites.
      hetBafSd: args.bools.has('no-baf-sd') ? undefined : a.profile.hetBafSd,
    } as never, stageOpts())
    out({ id: a.id, callRate: a.profile.callRate, hetRate: a.profile.hetRate,
      hetBand: a.profile.hetBand, hetBafSd: a.profile.hetBafSd, ...s }, () => {
      process.stdout.write(`${a.id}\n`)
      process.stdout.write(`  call rate        ${pct(a.profile.callRate, 2)}\n`)
      process.stdout.write(`  heterozygous     ${pct(a.profile.hetRate, 2)}\n`)
      process.stdout.write(`  BAF band         ${pct(a.profile.hetBand, 2)}\n`)
      process.stdout.write(`  stage            ${s.stage}\n`)
      process.stdout.write(`  templates        ${s.templates}\n`)
      process.stdout.write(`  dropout          ${Number.isFinite(s.dropout) ? s.dropout.toFixed(3) : 'not assigned'} (${s.basis})\n`)
      process.stdout.write(`  marker floor     ${s.markerFloor}\n`)
      process.stdout.write(`\n  ${s.why}.\n`)
      if (s.caveat) process.stdout.write(`\n  Limits: ${s.caveat}.\n`)
    })
  },

  link() {
    const [refPath, ...rest] = args.positional
    if (!refPath || !rest.length) die('usage: om link <parent> <sample>... [--json]')
    const stride = num('stride', 1)
    const ref = load(refPath, stride, false)
    const oppMax = num('opposite-hom-max', 0.020)
    const rows = rest.map((p) => {
      const s = load(p, stride, false)
      const opp = oppositeHom(ref.gt, s.gt)
      const st = stageMod.inferStage(s.profile, stageOpts())
      const link = secondParent(ref.gt, s.gt)
      const related = opp.rate <= oppMax && opp.n >= 10_000
      return {
        id: s.id,
        stage: st.stage,
        oppositeHom: opp.rate,
        oppositeHomMarkers: opp.n,
        secondParentFraction: link.fraction,
        informative: link.informative,
        relationship: !related ? 'unrelated'
          : link.ploidy === 'biparental' ? 'child'
            : link.ploidy === 'uniparental' ? 'this parent only, a duplicate or a haploid product'
              : 'ambiguous',
        why: link.why,
      }
    })
    out({ reference: ref.id, oppositeHomMax: oppMax, samples: rows }, () => {
      process.stdout.write(`reference ${ref.id}\n\n`)
      for (const r of rows) {
        process.stdout.write(`${r.id}  ${r.relationship}\n`)
        process.stdout.write(`   opposite-hom ${r.oppositeHom.toFixed(4)} over ${r.oppositeHomMarkers}`
          + `, second contribution ${pct(r.secondParentFraction)} of ${r.informative}`
          + `, ${r.stage}\n`)
      }
    })
  },

  async origin() {
    const [refPath, samplePath] = args.positional
    if (!refPath || !samplePath) {
      die('usage: om origin <parent> <sample> [--role paternal|maternal] [--json]')
    }
    const role = (args.flags.get('role') ?? 'paternal') as 'paternal' | 'maternal'
    if (role !== 'paternal' && role !== 'maternal') die('--role must be paternal or maternal')

    // THE SAME RUN THE BROWSER PERFORMS. Every channel, every guard and every refusal comes from
    // scoreSample, so this cannot answer differently from the app on the same file.
    const verbose = args.bools.has('verbose')
    const { id, refId, result } = await scoreArray(refPath, samplePath, role, (tag, text) => {
      if (verbose && !JSON_OUT) process.stderr.write(`  [${tag}] ${text}\n`)
    })

    const st = result.stage as { stage: string, dropout: number, basis: string } | undefined
    const rows = originRows(result as never, role)

    out({
      reference: refId, sample: id, role,
      stage: st?.stage, dropout: st?.dropout, dropoutSource: st?.basis,
      originClass: (result as { originClass?: string }).originClass,
      zygosity: (result as { zygosity?: string }).zygosity,
      events: rows,
    }, () => {
      process.stdout.write(`${id} against ${refId} as the ${role} parent\n`)
      process.stdout.write(`  material ${st?.stage ?? 'unknown'}`
        + `, dropout ${Number.isFinite(st?.dropout) ? st!.dropout.toFixed(3) : 'not assigned'}`
        + ` (${st?.basis ?? 'none'})\n`)
      process.stdout.write(`  genome ${(result as { originClass?: string }).originClass}`
        + `, ${(result as { zygosity?: string }).zygosity}\n\n`)
      if (!rows.length) { process.stdout.write('  no chromosomal change found\n'); return }
      for (const r of rows) {
        process.stdout.write(`${r.locus}  ${r.kind}\n`)
        process.stdout.write(`   ${r.origin ? `${r.origin.toUpperCase()} copy lost` : r.verdict}`
          + `  [${r.channel}${r.band ? `, band ${r.band}` : ''}`
          + `${Number.isFinite(r.confidence) ? `, ${r.confidence!.toFixed(3)}` : ''}]\n`)
        process.stdout.write(`   ${r.why}\n\n`)
      }
    })
  },

  async cohort() {
    const dir = args.positional[0]
    const refPath = args.flags.get('ref')
    if (!dir || !refPath) die('usage: om cohort <dir> --ref <array> [--role R] [--json]')
    const role = (args.flags.get('role') ?? 'paternal') as 'paternal' | 'maternal'
    const other = role === 'paternal' ? 'maternal' : 'paternal'
    const screenStride = num('screen-stride', 8)
    const oppMax = num('opposite-hom-max', 0.020)
    const ref = load(refPath)
    const parent = loadParent(refPath)
    const rs = stageMod.inferStage(ref.profile, stageOpts())
    // The reference gate is CALL RATE, not the stage label. The label is derived from
    // heterozygosity against a population anchor, which measures the individual as much as the
    // material, and it refused two genuine parental arrays: an egg donor's bulk gDNA at a 0.9837
    // call rate, and a cumulus sample at 0.951, both labelled single-cell.
    const minCall = num('ref-min-call', 0.90)
    if (rs.stage === 'failed' || ref.profile.callRate < minCall) {
      die(`reference ${ref.id} is not usable: call ${ref.profile.callRate.toFixed(4)}, `
        + `het ${ref.profile.hetRate.toFixed(4)}, stage ${rs.stage}`)
    }

    const files = arraysUnder(dir)
    const results: Record<string, unknown>[] = []
    const excluded: Record<string, unknown>[] = []
    let children = 0
    for (const f of files) {
      const scr = load(join(dir, f), screenStride, false)
      const sst = stageMod.inferStage(scr.profile, stageOpts())
      if (sst.stage === 'failed') {
        excluded.push({ id: scr.id, callRate: scr.profile.callRate, hetRate: scr.profile.hetRate })
        if (!JSON_OUT) process.stderr.write(`  excluded ${scr.id}: ${sst.why}\n`)
        continue
      }
      const opp = oppositeHom(ref.gt, scr.gt)
      if (!(opp.rate <= oppMax) || opp.n < 10_000) continue
      const c = load(join(dir, f))
      const link = secondParent(ref.gt, c.gt)
      if (link.ploidy !== 'biparental') {
        if (!JSON_OUT) process.stderr.write(`  set aside ${c.id}: ${link.why}\n`)
        continue
      }
      children += 1
      // THE SAME RUN THE BROWSER PERFORMS, per sample. A cohort summary assembled from a second
      // implementation is a summary of a different tool.
      const scored = await scoreWithParent(parent, join(dir, f), role)
      const st = scored.result.stage as { stage: string, dropout: number } | undefined
      const rows = originRows(scored.result as never, role)
      if (!JSON_OUT) {
        process.stderr.write(`  ${c.id}: child, ${st?.stage}, ${rows.length} event(s)\n`)
      }
      for (const r of rows) {
        results.push({
          reference: ref.id, sample: c.id, stage: st?.stage, dropout: st?.dropout, kind: r.kind,
          locus: r.locus, origin: r.origin, verdict: r.verdict,
          posterior: r.confidence, markers: r.markers, exclusive: undefined,
          heterozygous: undefined, why: r.why,
        })
      }
    }
    const tally: Record<string, number> = {}
    for (const r of results) tally[r.verdict as string] = (tally[r.verdict as string] ?? 0) + 1
    const payload = {
      dir, reference: ref.id, role, arrays: files.length, children,
      excluded: excluded.length, events: results.length, tally, results,
    }
    const o = args.flags.get('out')
    if (o) writeFileSync(o, JSON.stringify(payload, null, 1))
    out(payload, () => {
      process.stdout.write(`\n${files.length} arrays, ${excluded.length} excluded as not-a-genome, `
        + `${children} confirmed children, ${results.length} event(s)\n`)
      for (const [k, v] of Object.entries(tally)) process.stdout.write(`  ${k}: ${v}\n`)
      for (const r of results) {
        process.stdout.write(`  ${r.sample} ${r.locus} ${r.kind} -> `
          + `${r.origin ?? r.verdict} (${r.markers} informative, ${r.exclusive} exclusive)\n`)
      }
      if (o) process.stdout.write(`\nwritten to ${o}\n`)
    })
  },

  census() {
    const dir = args.positional[0]
    if (!dir) die('usage: om census <dir> [--json]')
    const stride = num('stride', 8)
    const minProducts = num('min-products', inferredRef.MIN_PRODUCTS)
    const rows = arraysUnder(dir).map((f) => {
      const a = load(join(dir, f), stride, false)
      const s = stageMod.inferStage(a.profile, stageOpts())
      return {
        id: a.id, group: a.id.replace(/[-_].*$/, ''), stage: s.stage,
        callRate: a.profile.callRate, hetRate: a.profile.hetRate,
      }
    })
    const groups = new Map<string, typeof rows>()
    for (const r of rows) groups.set(r.group, [...(groups.get(r.group) ?? []), r])
    const summary = [...groups].map(([g, rs]) => ({
      group: g, arrays: rs.length,
      haploid: rs.filter((r) => r.stage === 'haploid').length,
      failed: rs.filter((r) => r.stage === 'failed').length,
    })).sort((a, b) => b.haploid - a.haploid)
    const usable = summary.filter((s) => s.haploid >= minProducts).map((s) => s.group)
    out({ dir, minProducts, usable, summary, rows }, () => {
      process.stdout.write('group        arrays  haploid  failed\n')
      for (const s of summary) {
        process.stdout.write(`${s.group.padEnd(12)} ${String(s.arrays).padStart(6)} `
          + `${String(s.haploid).padStart(8)} ${String(s.failed).padStart(7)}`
          + `${s.haploid >= minProducts ? '   enough to reconstruct' : ''}\n`)
      }
      process.stdout.write(`\n${usable.length} group(s) with ${minProducts}+ haploid products.\n`)
      process.stdout.write('A group holds products of BOTH parents and a reconstruction needs '
        + `${minProducts} from ONE. Sort them first with: om reconstruct <products...> --group-only\n`)
    })
  },

  reconstruct() {
    const paths = args.positional
    if (paths.length < 2) die('usage: om reconstruct <product>... [--group-only] [--out F] [--json]')
    // full, not screened: the sort needs marker POSITIONS, and a screening load does not
    // populate them. Loading them without it silently fed zero markers to every product,
    // which put four products of one father into four separate groups.
    const loaded = paths.map((p) => load(p, num('stride', 1), true))
    for (const l of loaded) {
      const s = stageMod.inferStage(l.profile, stageOpts())
      if (s.stage === 'failed') die(`${l.id} is not a genome: ${s.why}`)
      if (s.stage !== 'haploid' && !args.bools.has('allow-diploid')) {
        die(`${l.id} stages as ${s.stage}, not haploid. A reconstruction takes HAPLOID products `
          + 'of one parent; a diploid compared with a haploid has a different expected rate '
          + 'entirely and bridges unrelated groups. Pass --allow-diploid to override')
      }
    }
    // Sorting first, because a group holds both parents' products and mixing them is garbage.
    const ps = new inferredRef.ProductSet()
    const slots = loaded.map((l) => {
      const slot = ps.begin(l.id)
      const band = { inBand: 0, total: 0 }
      for (const [probe, p] of l.pos) {
        ps.add(slot, {
          probesetId: probe, chrom: p.chrom, pos: p.pos, genotype: l.gt.get(probe) ?? 'NC',
          baf: null, log2R: null, copyNumber: null,
        } as never, band)
      }
      ps.end(slot, band)
      return slot
    })
    const rates = new Map<string, { shared: number, opposite: number, rate: number }>()
    const rateOf = (a: number, b: number) => {
      const k = a < b ? `${a}:${b}` : `${b}:${a}`
      const hit = rates.get(k)
      if (hit) return hit
      const v = ps.opposite(slots[a], slots[b]) as { shared: number, opposite: number, rate: number }
      rates.set(k, v)
      return v
    }
    for (const [k, v] of [...rates]) void k, v
    for (let a = 0; a < loaded.length; a += 1) {
      for (let b = a + 1; b < loaded.length; b += 1) {
        const v = rateOf(a, b)
        if (!Number.isFinite(v.rate) || v.shared < 1_000) {
          die(`${loaded[a].id} and ${loaded[b].id} share only ${v.shared} called markers, which `
            + 'is too few to say whether they came from the same parent')
        }
      }
    }
    const groups = inferredRef.groupByParent(loaded.length, (a: number, b: number) => rateOf(a, b).rate)
    const named = groups.map((g: number[]) => g.map((i) => loaded[i].id))
    const minProducts = num('min-products', inferredRef.MIN_PRODUCTS)

    if (args.bools.has('group-only')) {
      const detail = []
      for (let a = 0; a < loaded.length; a += 1) {
        for (let b = a + 1; b < loaded.length; b += 1) {
          const v = rateOf(a, b)
          detail.push({ a: loaded[a].id, b: loaded[b].id, rate: v.rate, shared: v.shared,
            kinship: inferredRef.kinship(v.rate) })
        }
      }
      out({ groups: named, minProducts, pairs: detail }, () => {
        named.forEach((g, i) => process.stdout.write(
          `group ${i + 1} (${g.length} product${g.length === 1 ? '' : 's'})`
          + `${g.length >= minProducts ? '' : `, under the ${minProducts} a reconstruction needs`}`
          + `: ${g.join(', ')}\n`,
        ))
        if (args.bools.has('pairs')) {
          process.stdout.write('\npairwise opposite-homozygote rate, the evidence for the split\n')
          for (const d of detail) {
            process.stdout.write(`  ${d.a} vs ${d.b}  ${d.rate.toFixed(4)} over ${d.shared}`
              + `  ${d.kinship}\n`)
          }
        }
      })
      return
    }
    const biggest = named.reduce((a, b) => (b.length > a.length ? b : a), [] as string[])
    if (biggest.length < minProducts) {
      die(`the largest single-parent group has ${biggest.length} product(s), under the `
        + `${minProducts} this needs. Below that the method inverts and true offspring read as `
        + `decisively absent. Groups found: ${named.map((g) => g.length).join(', ')}`)
    }
    const keep = new Set(biggest)
    const exclude = loaded.filter((l) => !keep.has(l.id)).map((l) => l.id)
    const pick = ps.chooseM(exclude) as { mMin: number }
    const mMin = num('m-min', pick.mMin)
    const ref = ps.build(mMin, exclude) as {
      genotype: Map<string, AB>, markers: number, hRetained: number, contamination: number,
      spuriousAbsence: number, nProducts: number, meanM: number
    }
    const o = args.flags.get('out')
    if (o) {
      const lines = ['probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset']
      const posOf = new Map<string, { chrom: string, pos: number }>()
      for (const l of loaded) for (const [p, v] of l.pos) if (!posOf.has(p)) posOf.set(p, v)
      for (const [probe, g] of ref.genotype) {
        const p = posOf.get(probe)
        if (!p) continue
        lines.push(`${probe}\t${p.chrom}\t${p.pos}\t\t\t2\t${g}\t1`)
      }
      writeFileSync(o, `${lines.join('\n')}\n`)
    }
    out({
      groups: named, used: biggest, mMin, markers: ref.markers, hRetained: ref.hRetained,
      contamination: ref.contamination, spuriousAbsence: ref.spuriousAbsence, written: o ?? null,
    }, () => {
      process.stdout.write(`${named.length} parental group(s): `
        + `${named.map((g) => g.length).join(', ')} product(s)\n`)
      process.stdout.write(`reconstructed from ${biggest.length}: ${biggest.join(', ')}\n\n`)
      process.stdout.write(`  m minimum        ${mMin}\n`)
      process.stdout.write(`  markers          ${ref.markers}\n`)
      process.stdout.write(`  parent het       ${pct(ref.hRetained, 2)}\n`)
      process.stdout.write(`  contamination    ${pct(ref.contamination, 3)} of the reference is a `
        + 'heterozygous site mistaken for homozygous\n')
      process.stdout.write(`  spurious absence ${pct(ref.spuriousAbsence, 3)} a haploid true `
        + 'offspring reads at these markers\n')
      if (o) process.stdout.write(`\nwritten to ${o}, usable as the parent in om origin or om cohort\n`)
    })
  },

  enrich() {
    const regionsPath = args.positional[0]
    const trackPath = args.flags.get('track')
    if (!regionsPath || !trackPath) {
      die('usage: om enrich <regions.tsv> --track <features.json> [--permutations N] [--json]')
    }
    const track = JSON.parse(readFileSync(trackPath, 'utf8'))
    const lines = readFileSync(regionsPath, 'utf8').trim().split('\n')
    const head = lines[0].split('\t')
    const col = (n: string) => head.indexOf(n)
    const regions = lines.slice(1).map((l) => {
      const f = l.split('\t')
      return { chrom: f[col('chr')], startBp: Number(f[col('start_bp')]), endBp: Number(f[col('end_bp')]) }
    })
    const markers = new Map<string, number[]>()
    const md = args.flags.get('markers')
    if (md) {
      const a = load(md)
      for (const p of a.pos.values()) markers.set(p.chrom, [...(markers.get(p.chrom) ?? []), p.pos])
      for (const v of markers.values()) v.sort((x, y) => x - y)
    } else {
      die('--markers <array> is required: the null is matched on THIS platform\'s informative '
        + 'markers per chromosome, and a uniform null reports an enrichment for anything that '
        + 'tracks marker density, which gene density does and fragile sites therefore do')
    }
    const perms = num('permutations', features.DEFAULT_PERMUTATIONS)
    const scored = features.scoreAll(track, regions as never, markers, perms)
    out({ regions: regions.length, permutations: perms, results: scored }, () => {
      process.stdout.write(`${regions.length} regions, ${perms} permutations\n\n`)
      // The field is nullMean, not expected. It was read as `expected` here, which is undefined,
      // so this printer threw on every run and only --json ever worked.
      for (const e of scored as { feature: string, observed: number, nullMean: number, ratio: number, p: number }[]) {
        const n = (v: number) => (Number.isFinite(v) ? v.toFixed(3) : '   n/a')
        process.stdout.write(`${e.feature.padEnd(38)} obs ${n(e.observed)} `
          + `null ${n(e.nullMean)}  ratio ${Number.isFinite(e.ratio) ? e.ratio.toFixed(2) : ' n/a'}`
          + `  p ${Number.isFinite(e.p) ? e.p.toFixed(4) : '   n/a'}\n`)
      }
      process.stdout.write('\nRead the effect sizes rather than the p values. At large region '
        + 'counts the permutation null is narrow enough that a four percent difference reaches '
        + 'p 0.0005.\n')
    })
  },

  // Every class this platform can and cannot see, so the limits are as inspectable as the dials.
  // A user deciding whether this tool answers their question should not have to run it to find out
  // that their question is unanswerable on this material.
  taxonomy() {
    const rows = tax.TAXONOMY.map((t: Record<string, string>) => ({
      class: t.cls, label: t.label, detectable: t.detectable, by: t.by, origin: t.origin,
      limit: t.limit ?? null,
      originReachable: !tax.ORIGIN_UNREACHABLE.has(t.cls),
    }))
    out({ classes: rows }, () => {
      for (const r of rows) {
        const mark = r.detectable === 'yes' ? ' ' : r.detectable === 'partly' ? '~' : 'X'
        process.stdout.write(`${mark} ${r.label}\n`)
        process.stdout.write(`    by      ${r.by}\n`)
        process.stdout.write(`    origin  ${r.origin}\n`)
        if (r.limit) process.stdout.write(`    LIMIT   ${r.limit}\n`)
        process.stdout.write('\n')
      }
      const hard = rows.filter((r) => r.detectable === 'no').length
      process.stdout.write(`${rows.length} classes, ${hard} of them not answerable on this `
        + 'platform at any quality. Those are listed rather than omitted so that silence is never '
        + 'read as absence.\n')
    })
  },

  constants() {
    const rows = [
      ['stage', 'BULK_HETEROZYGOSITY', stageMod.BULK_HETEROZYGOSITY, 'diploid rate measured on bulk gDNA on this platform. A PANEL and ANCESTRY property, not a constant of nature'],
      ['stage', 'QC_CALL_FLOOR', stageMod.QC_CALL_FLOOR, 'below this no stage is inferred: a failed reaction imitates haploid material'],
      ['stage', 'MAX_DIPLOID_HET', stageMod.MAX_DIPLOID_HET, 'above this the array is not one genome. 0.168 plus drop-in at the top of the measured range gives 0.212'],
      ['stage', 'HAPLOID_MAX_HET', stageMod.HAPLOID_MAX_HET, 'below this the sample carries one genome. NOT safe in both directions'],
      ['stage', 'BAND_DIPLOID_CERTAIN', stageMod.BAND_DIPLOID_CERTAIN, 'BAF band at or above which the sample is diploid whatever its calls say. Vetoes haploid, cannot confirm it'],
      ['origin', 'DEFAULT_Q', oneParent.DEFAULT_Q, 'frequency of the allele the loaded parent lacks'],
      ['origin', 'GENOTYPE_ERROR', oneParent.GENOTYPE_ERROR, 'a homozygote read as the other homozygote'],
      ['origin', 'DROP_IN', oneParent.DROP_IN, 'false heterozygote at a truly homozygous marker, measured over 113 same-genome pairs'],
      ['origin', 'MIN_MARKERS', oneParent.MIN_MARKERS, 'informative markers a region needs before any verdict'],
      ['origin', 'MAX_REGION_HET', oneParent.MAX_REGION_HET, 'above this the region\'s genotypes are not measuring it and no origin is called'],
      ['origin', 'CALL_POSTERIOR', oneParent.CALL_POSTERIOR, 'posterior a hypothesis must reach to be named'],
      ['dosage', 'ORIGIN_WITHOUT_CLASS.trophectoderm', dosage.ORIGIN_WITHOUT_CLASS.trophectoderm, 'fraction of detected TE events that resolve an ORIGIN but not a CLASS: emitting one verdict would refuse all of them'],
      ['posterior', 'BAND_A_MIN', post.BAND_A_MIN, 'calibrated confidence for the top band, measured accuracy 0.9952-0.9980 across the four material classes'],
      ['posterior', 'BAND_B_MIN', post.BAND_B_MIN, 'confident band floor, measured 0.9448-0.9599'],
      ['posterior', 'BAND_C_MIN', post.BAND_C_MIN, 'weak band floor, measured 0.8128-0.8375; below it band D at 0.604-0.636, which still excludes chance'],
      ['posterior', 'VETO_MAX_F', post.VETO_MAX_F, 'implied GAIN fraction under which an amplified event is class-inverted: a small gain and a small loss name OPPOSITE parents, measured 0 of 34 correct on TE and 0 of 18 on blastomere'],
      ['posterior', 'BAND_ACCURACY.blastomere.A', post.BAND_ACCURACY.blastomere.A, 'why the old 0.30 fraction withhold was retired: 24.1% of blastomere events it refused sit in this band'],
      ['dosage', 'Z_DETECT', dosage.Z_DETECT, 'self-referenced |z| an imbalance must reach, two-sided at 1%'],
      ['stage', 'BULK_MAX_BAF_SD', stageMod.BULK_MAX_BAF_SD, 'BAF spread at het calls above which the BULK rung is refused: heterozygosity says how many heterozygotes survived, this says how far they scattered, and 176 of 877 arrays were called bulk on the first while scattering like amplified material'],
      ['dosage', 'MAX_HET_BAF_SD', dosage.MAX_HET_BAF_SD, 'array-level gate adopted from MoChA: BAF spread at het sites above this and the array is not analysed'],
      ['dosage', 'WINDOW_LO', dosage.WINDOW_LO, 'central window, deliberately wider than the old middle band'],
      ['dosage', 'DRIFT_TAU.blastomere', dosage.DRIFT_TAU.blastomere, 'within-array drift, the FLOOR on the standard error; it does not average down with markers'],
      ['dosage', 'VIF_CHROMOSOME.blastomere', dosage.VIF_CHROMOSOME.blastomere, 'variance inflation from spatial correlation; bulk is 0.94, amplified material is not white'],
      ['dosage', 'RESIDUAL_R.trophectoderm', dosage.RESIDUAL_R.trophectoderm, 'measured correlation between the dosage and intensity channels on THIS statistic, over 81 arrays; bulk is -0.058 and quadrature would overstate the joint z by a quarter on TE'],
      ['linkage', 'ONE_PARENT_HAPLOID_MAX', obligate.ONE_PARENT_HAPLOID_MAX, 'under this, one parent\'s genome and nothing else'],
      ['linkage', 'ONE_PARENT_DIPLOID_MIN', obligate.ONE_PARENT_DIPLOID_MIN, 'over this, a second parental contribution is present'],
      ['reconstruct', 'MIN_PRODUCTS', inferredRef.MIN_PRODUCTS, 'below this the method INVERTS and true offspring read as decisively absent'],
      ['reconstruct', 'MIN_ASCERTAINMENT', inferredRef.MIN_ASCERTAINMENT, 'recovery plateaus at m=4 while contamination keeps falling, so 0.90 sits at the knee'],
      ['reconstruct', 'SAME_PARENT_MAX', inferredRef.SAME_PARENT_MAX, 'opposite-hom rate under which two haploid products share a parent'],
      ['reconstruct', 'DIFFERENT_PARENT_MIN', inferredRef.DIFFERENT_PARENT_MIN, 'over this they do not. Between the two is ambiguous, deliberately'],
      ['enrich', 'DEFAULT_PERMUTATIONS', features.DEFAULT_PERMUTATIONS, 'null intervals drawn per feature, matched on informative markers'],
      ['enrich', 'MIN_REGIONS', features.MIN_REGIONS, 'below this an enrichment is not computed'],
    ] as [string, string, number, string][]
    out(Object.fromEntries(rows.map(([g, n, v, w]) => [n, { group: g, value: v, why: w }])), () => {
      let last = ''
      for (const [g, n, v, w] of rows) {
        if (g !== last) { process.stdout.write(`\n${g}\n`); last = g }
        process.stdout.write(`  ${n.padEnd(24)} ${String(v).padEnd(8)} ${w}\n`)
      }
      process.stdout.write('\nEvery one is a flag on the command that uses it, lower-cased with '
        + 'dashes: --max-region-het, --min-products, --drop-in. Moving one is recorded in the '
        + 'output, because a number produced under a changed constant is not comparable with the '
        + 'validation figures.\n')
    })
  },

  help() {
    process.stdout.write(`om, OriginMarker on the command line

  om stage       <array>                  material, dropout, marker floor
  om link        <parent> <sample>...     parent, duplicate, unrelated, or refused
  om origin      <parent> <sample>        which parent's copy is missing, per region
  om cohort      <dir> --ref <array>      origin for every confirmed child in a directory
  om census      <dir>                    haploid products per group
  om reconstruct <product>...             a parent's genotypes from that parent's haploid cells
  om enrich      <regions.tsv> --track T --markers A
  om taxonomy                             every abnormality class, what it rests on and its limits
  om constants                            every tunable, its value and its provenance

Common flags
  --json                machine-readable output
  --stride N            subsample every Nth marker. Screening only
  --ado N               force the dropout instead of inferring it from the stage
  --role paternal|maternal     which parent the loaded array is
  --region chr6:39302-294904   score one interval instead of the detected ones
  --out FILE            write the result, or the reconstructed array, to a file

Every constant in \`om constants\` is a flag on the command that uses it. The browser tool at
originmarker.app infers all of them and offers no knobs, on purpose. This does the opposite.
`)
  },
}

const run = COMMANDS[args.cmd] ?? COMMANDS.help
if (args.bools.has('help') && args.cmd !== 'help') { COMMANDS.help(); process.exit(0) }
await run()
