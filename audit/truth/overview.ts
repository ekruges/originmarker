/**
 * The run's overview table (web/src/overview.ts) built from per-sample records, printed as the
 * page's download prints it.
 *
 * A record is one line of audit/truth/run.ts output: `stage`, `refusals` and `card` (the result
 * card's counted rows with their pieces). Where the run grouped cells, a record also carries
 * `embryo` (its group) and each piece its `mechanism` (timing.ts), as the page attaches them.
 * One invocation is one run. audit/truth/overview_check.py recounts the same records without this
 * code.
 *
 * Run: node --experimental-strip-types audit/truth/overview.ts <records.jsonl> [more.jsonl ...]
 */
import { readFileSync } from 'node:fs'

const W = new URL('../../web/src/', import.meta.url).pathname
const overview = await import(`${W}overview.ts`)
const { normaliseTrack } = await import(`${W}comparison.ts`)

const track = normaliseTrack(JSON.parse(readFileSync(new URL('../../web/public/hg19_features.json', import.meta.url), 'utf8')))
const arms = overview.armsFrom(track!.gaps)

const samples = process.argv.slice(2).flatMap((f) => readFileSync(f, 'utf8').split('\n').filter((l) => l.startsWith('{')))
  .map((l) => JSON.parse(l)).filter((r) => r.card)
  .map((r) => overview.overviewSampleOf({
    counted: r.card.counted.map((g: any) => ({
      chrom: g.chrom, gain: g.gain, whole: g.whole, origin: g.origin, startBp: g.start, endBp: g.end,
      pieces: g.pieces.map((p: any) => ({ kind: p.kind, startBp: p.start, endBp: p.end, mechanism: p.mechanism })),
    })),
    once: Array(r.card.once).fill(null),
    seen: [],
    disturbed: r.card.disturbed ? { on: NaN, of: NaN } : undefined,
  } as never, r, r.embryo))

process.stdout.write(overview.overviewTsv(overview.overviewTable(samples, arms)))
