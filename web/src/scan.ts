// The two pieces of a Syngamy run whose cost grows with the size of the array.
//
// WHY THEY LIVE HERE RATHER THAN IN THE COMPONENT. Both of these were written inline in
// Syngamy.tsx, and one of them was quadratic: it rescanned every marker of a chromosome inside its
// own window loop, which is about 28 million operations and 216 full-length allocations on chr1
// alone. It shipped, and it locked the tab for the whole genome on every run. No check in this
// repo could have caught it, because nothing outside a browser could reach code that only exists
// inside a React component body. Moving the scaling work out is what makes `scan.check.ts` able to
// assert that these stay linear in the number of markers.
//
// Nothing here decides anything. Detection, classification and origin all read these outputs.

/**
 * Hand the event loop back so the page can paint the line just logged.
 *
 * NOT requestAnimationFrame. It does not fire AT ALL in a hidden tab, so a run the user tabbed
 * away from waits forever, which is a worse failure than the slow scan this was added to fix: at
 * least the slow one finished. A chained setTimeout is no good either, throttled to one a second
 * in a background tab. A MessageChannel message is neither paused nor throttled and still returns
 * to the event loop, so the page stays live whether or not anyone is looking at it.
 */
export const breathe = (): Promise<void> => new Promise((resolve) => {
  if (typeof MessageChannel !== 'function') { setTimeout(resolve, 0); return }
  const ch = new MessageChannel()
  ch.port1.onmessage = () => { ch.port1.close(); resolve() }
  ch.port2.postMessage(0)
})

export type CnMarker = { chrom: string; pos: number; called: boolean; log2R: number | null }
export type SelfMarker = { chrom: string; pos: number; het: boolean }
export type Window = {
  chrom: string; startBp: number; endBp: number; called: number; het: number
  logR?: number; wholeChromosome: boolean
}

const medianOf = (xs: number[]): number | undefined => (xs.length
  ? xs.slice().sort((a, b) => a - b)[xs.length >> 1] : undefined)

/** Last position seen on each chromosome. */
export const chromEnds = (cnByChrom: Map<string, CnMarker[]>): Map<string, number> => {
  const out = new Map<string, number>()
  for (const [c, ms] of cnByChrom) out.set(c, ms.reduce((a, m) => Math.max(a, m.pos), 0))
  return out
}

/**
 * Sliding windows AND whole chromosomes, so a copy-neutral event covering part of a chromosome is
 * reported at its own extent rather than as the whole thing. The whole-chromosome pass is kept
 * because it is the one with a detection floor on amplified material: a 12 Mb-scale interval has
 * none at any fraction with one parent.
 *
 * ONE PASS PER CHROMOSOME, NOT ONE PER WINDOW. Markers are sorted once, heterozygote counts come
 * from a prefix sum, and each window's intensity is read from a binary-searched slice. `onChrom`
 * lets the caller hand the page back between chromosomes; without it the whole genome is one task
 * and the tab is unresponsive for its duration however fast the arithmetic is.
 */
export const copyNeutralWindows = async (
  cnByChrom: Map<string, CnMarker[]>,
  selfMarkers: SelfMarker[],
  genomeLrr: number,
  segmentMarkers: number,
  onChrom?: (chrom: string) => Promise<void> | void,
): Promise<{ windows: Window[]; scanned: number; chromEnd: Map<string, number> }> => {
  const chromEnd = chromEnds(cnByChrom)
  const selfByChrom = new Map<string, SelfMarker[]>()
  for (const m of selfMarkers) {
    const cur = selfByChrom.get(m.chrom)
    if (cur) cur.push(m)
    else selfByChrom.set(m.chrom, [m])
  }
  const windows: Window[] = []
  let scanned = 0
  for (const [c, msRaw] of cnByChrom) {
    const self = (selfByChrom.get(c) ?? []).slice().sort((a, b) => a.pos - b.pos)
    const ms = msRaw.slice().sort((a, b) => a.pos - b.pos)
    const msPos = ms.map((m) => m.pos)
    const lrrAllChrom = ms.map((m) => m.log2R).filter((x): x is number => x !== null)
    windows.push({
      chrom: c, startBp: 0, endBp: chromEnd.get(c) ?? 0,
      called: self.length, het: self.filter((m) => m.het).length,
      logR: lrrAllChrom.length ? (medianOf(lrrAllChrom) as number) - genomeLrr : undefined,
      wholeChromosome: true,
    })
    // Prefix sum of heterozygotes, so a window's count is a subtraction.
    const hetPrefix = new Int32Array(self.length + 1)
    for (let i = 0; i < self.length; i += 1) hetPrefix[i + 1] = hetPrefix[i] + (self[i].het ? 1 : 0)
    const lowerBound = (v: number) => {
      let lo = 0
      let hi = msPos.length
      while (lo < hi) { const mid = (lo + hi) >> 1; if (msPos[mid] < v) lo = mid + 1; else hi = mid }
      return lo
    }
    const step = Math.floor(segmentMarkers / 2)
    for (let i = 0; i + segmentMarkers <= self.length; i += step) {
      const lo = self[i].pos
      const hi = self[i + segmentMarkers - 1].pos
      const a = lowerBound(lo)
      const b = lowerBound(hi + 1)
      const inWin: number[] = []
      for (let k = a; k < b; k += 1) {
        const v = ms[k].log2R
        if (v !== null && Number.isFinite(v)) inWin.push(v)
      }
      windows.push({
        chrom: c, startBp: lo, endBp: hi,
        called: segmentMarkers,
        het: hetPrefix[i + segmentMarkers] - hetPrefix[i],
        logR: inWin.length ? (medianOf(inWin) as number) - genomeLrr : undefined,
        wholeChromosome: false,
      })
      scanned += 1
    }
    if (onChrom) await onChrom(c)
  }
  return { windows, scanned, chromEnd }
}

export type GatherInput = {
  markerPos: Map<string, { chrom: string; pos: number }>
  parentGt: Map<string, string>
  myBaf: Map<string, number>
  myGt: Map<string, string>
  cnByChrom: Map<string, CnMarker[]>
}

export type Interval = { chrom: string; startBp?: number; endBp?: number }

type ChromRows = {
  pos: Int32Array
  /** Genotype channel rows, sorted by position. Shared, never mutated, never rebuilt. */
  rows: [string, number | null][]
  /** Untransmitted channel rows, at the same positions, aligned to `untPos`. */
  untPos: Int32Array
  untRows: [string, string, number | null][]
  cnPos: Int32Array
  cnLogR: Float64Array
}

/**
 * The array flattened once per sample, grouped by chromosome and sorted by position.
 *
 * THIS IS THE DIFFERENCE BETWEEN A RUN AND A WAIT. Scoring an interval used to walk four Maps and
 * allocate a fresh two-element array for every marker it saw. Twenty-five findings on 825,000
 * markers is twenty million tiny allocations that all hold the same values, and it measured nine
 * seconds. Built here once, the rows are shared: an interval takes a slice of its own chromosome
 * and the background is the other chromosomes' slices laid end to end. Nothing downstream writes
 * to a row, and everything downstream reduces them to a mean, so sharing them is safe and their
 * order within a chromosome does not matter.
 */
export type ScanIndex = { byChrom: Map<string, ChromRows>; chroms: string[] }

export const buildScanIndex = (src: GatherInput): ScanIndex => {
  type Acc = {
    pos: number[]; rows: [string, number | null][]
    untPos: number[]; untRows: [string, string, number | null][]
    cnPos: number[]; cnLogR: number[]
  }
  const acc = new Map<string, Acc>()
  const get = (c: string): Acc => {
    let a = acc.get(c)
    if (!a) {
      a = { pos: [], rows: [], untPos: [], untRows: [], cnPos: [], cnLogR: [] }
      acc.set(c, a)
    }
    return a
  }
  // Only markers the loaded parent has a genotype for reach the genotype channel at all, so the
  // rest are dropped once here rather than skipped once per finding.
  for (const [probe, p] of src.markerPos) {
    const pg = src.parentGt.get(probe)
    if (!pg) continue
    const b = src.myBaf.get(probe) ?? null
    const a = get(p.chrom)
    a.pos.push(p.pos)
    a.rows.push([pg, b])
    if (pg === 'AB') {
      a.untPos.push(p.pos)
      a.untRows.push([pg, src.myGt.get(probe) ?? 'NC', b])
    }
  }
  for (const [ch, ms] of src.cnByChrom) {
    const a = get(ch)
    for (const m of ms) {
      if (m.log2R === null || !Number.isFinite(m.log2R)) continue
      a.cnPos.push(m.pos)
      a.cnLogR.push(m.log2R)
    }
  }
  const sorted = (pos: number[]) => pos.map((_, i) => i).sort((x, y) => pos[x] - pos[y])
  const byChrom = new Map<string, ChromRows>()
  for (const [c, a] of acc) {
    const o = sorted(a.pos)
    const u = sorted(a.untPos)
    const n = sorted(a.cnPos)
    byChrom.set(c, {
      pos: Int32Array.from(o, (i) => a.pos[i]),
      rows: o.map((i) => a.rows[i]),
      untPos: Int32Array.from(u, (i) => a.untPos[i]),
      untRows: u.map((i) => a.untRows[i]),
      cnPos: Int32Array.from(n, (i) => a.cnPos[i]),
      cnLogR: Float64Array.from(n, (i) => a.cnLogR[i]),
    })
  }
  return { byChrom, chroms: [...byChrom.keys()] }
}

export type Gathered = {
  region: readonly (readonly [string, number | null])[]
  background: readonly (readonly [string, number | null])[]
  inL: number[]
  outL: number[]
  untRows: readonly (readonly [string, string, number | null])[]
}

const lowerBound = (xs: Int32Array, v: number): number => {
  let lo = 0
  let hi = xs.length
  while (lo < hi) { const mid = (lo + hi) >> 1; if (xs[mid] < v) lo = mid + 1; else hi = mid }
  return lo
}

/**
 * Everything one interval needs off the array. Called once per event, so a run with twenty
 * findings is twenty of these, which is why it reads slices of a prebuilt index rather than the
 * array itself.
 *
 * `region` and `background` are the genotype channel. Background is the rest of this same array,
 * because self-referencing is not optional: the raw one-parent null sits at -0.031 on
 * trophectoderm under no event, pointing at the parent that was NOT genotyped, which is the shift
 * a real mosaic fraction of 0.117 would produce.
 *
 * `inL`/`outL` are the intensity channel, self-referenced from the copy-number readings. It
 * informs the STATE and never the ORIGIN: on haploid pronuclei log2R cannot tell maternal from
 * paternal at all.
 *
 * `untRows` is THE UNTRANSMITTED CHANNEL, on the disjoint marker set the obligate-het path
 * discards: markers where the loaded parent is HETEROZYGOUS and this sample reads homozygous, so
 * the transmission is determined. Every marker here is informative by construction, against 32-90%
 * in the parent-homozygous window, which is where its 1.40-1.98x advantage comes from. It is also
 * the only channel that gives a single blastomere a defined floor at all.
 */
export const gatherInterval = (index: ScanIndex, iv: Interval): Gathered => {
  const lo = iv.startBp ?? Number.NEGATIVE_INFINITY
  const hi = iv.endBp ?? Number.POSITIVE_INFINITY
  const whole = !Number.isFinite(lo) && !Number.isFinite(hi)
  const here = index.byChrom.get(iv.chrom)

  const region: (readonly [string, number | null])[] = []
  const background: (readonly [string, number | null])[] = []
  const untRows: (readonly [string, string, number | null])[] = []
  const inL: number[] = []
  const outL: number[] = []

  if (here) {
    const a = whole ? 0 : lowerBound(here.pos, lo)
    const b = whole ? here.rows.length : lowerBound(here.pos, hi + 1)
    for (let i = a; i < b; i += 1) region.push(here.rows[i])
    // A SUB-CHROMOSOMAL INTERVAL LEAVES THE REST OF ITS OWN CHROMOSOME IN THE BACKGROUND. Dropping
    // it is not a rounding difference: on an event covering half a chromosome it is half the
    // markers nearest the event, and they are the ones the background most needs.
    for (let i = 0; i < a; i += 1) background.push(here.rows[i])
    for (let i = b; i < here.rows.length; i += 1) background.push(here.rows[i])

    const ua = whole ? 0 : lowerBound(here.untPos, lo)
    const ub = whole ? here.untRows.length : lowerBound(here.untPos, hi + 1)
    for (let i = ua; i < ub; i += 1) untRows.push(here.untRows[i])

    const ca = whole ? 0 : lowerBound(here.cnPos, lo)
    const cb = whole ? here.cnLogR.length : lowerBound(here.cnPos, hi + 1)
    for (let i = ca; i < cb; i += 1) inL.push(here.cnLogR[i])
    for (let i = 0; i < ca; i += 1) outL.push(here.cnLogR[i])
    for (let i = cb; i < here.cnLogR.length; i += 1) outL.push(here.cnLogR[i])
  }
  // Every other chromosome, laid end to end. The rows are the ones built with the index, so this
  // copies references and allocates nothing per marker.
  for (const c of index.chroms) {
    if (c === iv.chrom) continue
    const o = index.byChrom.get(c)!
    for (let i = 0; i < o.rows.length; i += 1) background.push(o.rows[i])
    for (let i = 0; i < o.cnLogR.length; i += 1) outL.push(o.cnLogR[i])
  }
  return { region, background, inL, outL, untRows }
}
