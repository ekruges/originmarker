/**
 * What one sample's array says, from its marker rows.
 *
 * THIS IS THE ONE IMPLEMENTATION OF THE ANSWER. The browser and the command line both call it, so
 * there is no second copy of the reasoning for them to disagree over, and a check can run a real
 * array through exactly the code either surface runs.
 *
 * Both of those are load-bearing. The two surfaces previously reported different parents for the
 * same two losses on the same file: once because the command line assembled its input without an
 * allele-frequency channel, which moved the sample across the diploid heterozygosity boundary,
 * and once because a refusal verdict was renamed and only one surface's gate still matched the
 * old string. Neither was reachable by any check, because the browser's copy lived inside a React
 * component nothing could import.
 *
 * Nothing here touches React, the DOM or the filesystem. Rows arrive through `collectRow`,
 * whatever read them, and the narrative leaves through a `LogFn` the caller supplies.
 */
import { breathe, buildScanIndex, copyNeutralWindows, gatherInterval, type Interval } from './scan.ts'
import { type AB } from './informativity.ts'
import type { ProbeRow, SampleProfile } from './ingest.ts'
import {
  classify, emptyTally, isAutosome, pct, tallyRow, type ParentageResult,
} from './parentage.ts'
import {
  scanChromosome, scanCopyNumber, externalNull, segmentCoords, type MarkerAbsence,
} from './segments.ts'
import {
  paternalShare, recentre, callGainOrigin, callHomologue, externalHetBackground,
  MIN_INFORMATIVE_DEFAULT, type DosageMarker, callLossOrigin,
} from './gainOrigin.ts'
import { addTwoParent, emptyHet, hetCall, type HetTally } from './obligateHet.ts'
import { originBlockedByClass, namedAParent } from './defects.ts'
import { callSiblingOrigin, hetRule, type AB as SibAB } from './siblingOrigin.ts'
import { callOneParentOrigin } from './oneParentOrigin.ts'
import { callDosageOrigin, materialOf, originUnreachable } from './dosageOrigin.ts'
import { uniparentalOrigin } from './uniparentalOrigin.ts'
import {
  detectLoh, detectUpd, detectTriploidy, detectComplex, runsOfHomozygosity, mergeLoh,
  LOH_SEGMENT_MARKERS,
} from './abnormalities.ts'
import { untransmittedPairs, impossibleRate, orientUntransmitted, callMechanism } from './untransmitted.ts'
import { inferStage } from './stage.ts'

/** One parent's array, reduced to what every channel below reads. */
export interface ParentIndex {
  gt: Map<string, AB>
  heterozygosity: number
  build: string | null
}

/** Where a run writes its narrative. The browser paints it, the command line prints it. */
export type LogTag = 'READ' | 'PARSE' | 'CALL' | 'WARN' | 'DONE' | 'SCAN'
export type LogFn = (tag: LogTag, text: string) => void

/**
 * The chromosomes a list of interval labels covers, collapsed.
 *
 * A refusal that applies to three hundred intervals is one fact, and a reader wants to know how
 * far it reached rather than to scroll three hundred coordinates. The full list stays in the
 * report and the export.
 */
export const listChroms = (wheres: readonly string[]): string => {
  const seen: string[] = []
  for (const w of wheres) {
    const c = /^chr([\dXY]+)/.exec(w)?.[1]
    if (c && !seen.includes(c)) seen.push(c)
  }
  if (!seen.length) return `${wheres.length} intervals`
  const n = seen.length
  return n <= 6 ? `chr${seen.join(', chr')}` : `${n} chromosomes, chr${seen[0]} to chr${seen[n - 1]}`
}

/**
 * A parent array, accumulated one row at a time.
 *
 * HERE FOR THE SAME REASON AS `collectRow`. The command line used to build its own parent index,
 * keeping only called autosomal markers, while the browser kept every row. A marker the parent
 * index omits reads as a no-call downstream, so the same array produced different tallies on the
 * two surfaces before either channel had run.
 */
export interface ParentAccum { gt: Map<string, AB>; called: number; het: number }

export const emptyParent = (): ParentAccum => ({ gt: new Map(), called: 0, het: 0 })

export function collectParentRow(r: ProbeRow, acc: ParentAccum): void {
  acc.gt.set(r.probesetId, r.genotype)
  if (r.genotype !== 'NC' && isAutosome(r.chrom)) {
    acc.called += 1
    if (r.genotype === 'AB') acc.het += 1
  }
}

/** Heterozygosity is autosomal and over CALLED markers: it is the array's own error scale. */
export const finishParent = (acc: ParentAccum, build: string | null): ParentIndex => ({
  gt: acc.gt,
  heterozygosity: acc.called ? acc.het / acc.called : NaN,
  build,
})

/**
 * Everything one sample's rows accumulate into, in a single pass.
 *
 * One pass rather than several: these arrays run past 825,000 markers and re-reading one to
 * answer a second question costs as much as the first answer did.
 */
export interface Collected {
  t: ReturnType<typeof emptyTally>
  tm: ReturnType<typeof emptyTally> | null
  dosageByChrom: Map<string, DosageMarker[]>
  hetByChrom: Map<string, { informative: number; het: number }>
  selfMarkers: { chrom: string; pos: number; het: boolean }[]
  obligateByChrom: Map<string, HetTally>
  absenceByChrom: Map<string, MarkerAbsence[]>
  cnByChrom: Map<string, { chrom: string; pos: number; called: boolean; log2R: number | null }[]>
  myGt: Map<string, string>
  myBaf: Map<string, number>
  markerPos: Map<string, { chrom: string; pos: number }>
}

/**
 * A fresh accumulator for one sample.
 *
 * The parent's assembly is stamped on the tallies here, so the pseudoautosomal boundaries are the
 * right ones. It is known before the sample streams and is the same across one experiment.
 */
export function emptyCollected(pat: ParentIndex, mat?: ParentIndex | null): Collected {
  const t = emptyTally()
  const tm = mat ? emptyTally() : null
  t.build = pat.build
  if (tm) tm.build = mat?.build ?? null
  return {
    t,
    tm,
    // Allele dosage, for saying which parent an extra copy came from. Only computable where BOTH
    // parents are known and homozygous for different alleles, so this stays empty on a
    // single-donor run and the annotation says why rather than guessing.
    dosageByChrom: new Map(),
    // Heterozygosity at the loaded parent's homozygous markers, per chromosome. In a haploid
    // product every one of these is amplification error, and a region running far above the
    // array's own rate carries the parent's OTHER homologue.
    hetByChrom: new Map(),
    // Sample-only zygosity, for the classes that need no parent: LOH, runs, ploidy.
    selfMarkers: [],
    // Parental contribution per chromosome. With both parents this is the obligate set: both
    // homozygous and OPPOSITE, where a biparental cell must read heterozygous. With one parent
    // `hetByChrom` is already the same tally in its one-parent form.
    obligateByChrom: new Map(),
    absenceByChrom: new Map(),
    // The copy-number channel: every marker on the array, called or not, with its intensity. A
    // region that is GONE stops producing calls, which the parental-absence indicator cannot see,
    // because a no-call is excluded from it as uninformative.
    cnByChrom: new Map(),
    myGt: new Map(),
    // Allele dosage per marker, kept for EVERY marker rather than every CALLED marker. A no-call
    // still has an intensity reading, and on a chromosome whose call rate has collapsed those are
    // the only readings left. This is what the dosage channel runs on.
    myBaf: new Map(),
    markerPos: new Map(),
  }
}

/**
 * One marker row, into the accumulator.
 *
 * THE SINGLE PLACE A ROW IS READ, so the two surfaces cannot disagree about what a row means.
 * Every field of `acc` is a map, an array or a tally, so destructuring and mutating through the
 * local names is the same as mutating `acc`, and the body below is the browser's original.
 */
export function collectRow(
  r: ProbeRow, pat: ParentIndex, mat: ParentIndex | null | undefined, acc: Collected,
): void {
  const {
    t, tm, myGt, myBaf, markerPos, cnByChrom, selfMarkers, hetByChrom, obligateByChrom,
    absenceByChrom, dosageByChrom,
  } = acc
  if (isAutosome(r.chrom) && r.genotype !== 'NC') {
    myGt.set(r.probesetId, r.genotype)
    markerPos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  }
  const fa = pat.gt.get(r.probesetId) ?? 'NC'
  tallyRow(fa, r, t)
  if (tm && mat) tallyRow(mat.gt.get(r.probesetId) ?? 'NC', r, tm)
  // Only a marker where the parent is homozygous and the sample is called can say
  // anything about presence, which is the same informative set the rates use.
  if (isAutosome(r.chrom)) {
    markerPos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
    if (r.baf !== null && Number.isFinite(r.baf)) myBaf.set(r.probesetId, r.baf)
    const cn = cnByChrom.get(r.chrom) ?? []
    cn.push({ chrom: r.chrom, pos: r.pos, called: r.genotype !== 'NC', log2R: r.log2R })
    cnByChrom.set(r.chrom, cn)
    // The SAMPLE's own zygosity, which is a different tally from hetByChrom above: that
    // one counts markers where the PARENT is homozygous. Copy-neutral LOH and runs of
    // homozygosity are properties of the sample alone and need no parent at all.
    if (r.genotype !== 'NC') {
      selfMarkers.push({ chrom: r.chrom, pos: r.pos, het: r.genotype === 'AB' })
    }
  }
  if ((fa === 'AA' || fa === 'BB') && r.genotype !== 'NC' && isAutosome(r.chrom)) {
    // The homologue channel: every heterozygous call here is error in a haploid.
    const h = hetByChrom.get(r.chrom) ?? { informative: 0, het: 0 }
    h.informative += 1
    if (r.genotype === 'AB') h.het += 1
    hetByChrom.set(r.chrom, h)
    if (mat) {
      const mo = mat.gt.get(r.probesetId) ?? 'NC'
      const ob = obligateByChrom.get(r.chrom) ?? emptyHet()
      addTwoParent(fa, mo, r.genotype, ob)
      obligateByChrom.set(r.chrom, ob)
      const share = paternalShare(fa, mo, r.baf)
      if (share !== null) {
        const d = dosageByChrom.get(r.chrom) ?? []
        d.push({ chrom: r.chrom, pos: r.pos, patShare: share })
        dosageByChrom.set(r.chrom, d)
      }
    }
    const arr = absenceByChrom.get(r.chrom) ?? []
    arr.push({ chrom: r.chrom, pos: r.pos,
      absent: r.genotype !== 'AB' && r.genotype !== fa })
    absenceByChrom.set(r.chrom, arr)
  }
}

/**
 * Every channel, over one sample's collected rows.
 *
 * Order is not arbitrary. The classes that need no parent run first, because they are properties
 * of the sample alone and decide whether the parental channels may run at all. A measured
 * parental answer always outranks an inherited one, so the zygosity channel is consulted last.
 */
export async function scoreSample(input: {
  acc: Collected
  profile: SampleProfile
  pat: ParentIndex
  mat: ParentIndex | null | undefined
  soloRole: 'paternal' | 'maternal'
  /** Genotypes of the other arrays in this run, for the sibling-referenced channel. */
  sibs: readonly Map<string, string>[]
  sampleName: string
  log: LogFn
  /**
   * Stage-inference overrides.
   *
   * Exists so the command line can expose its thresholds as flags without owning a second copy of
   * the inference. Undefined here means the shipped defaults, which is what the browser passes.
   */
  stageOpts?: Parameters<typeof inferStage>[1]
}): Promise<ParentageResult> {
  const { acc, profile, pat, mat, soloRole, sibs, sampleName, log } = input
  const {
    t, myGt, myBaf, markerPos, cnByChrom, selfMarkers, hetByChrom, obligateByChrom,
    absenceByChrom, dosageByChrom,
  } = acc
  const result = classify(t, pat.heterozygosity, { role: soloRole })
  // Stage from the array itself, since the dropout each stage carries is what every
  // downstream likelihood is parameterised by. Bundled into the result so every output
  // carries it, with the basis and the confounds attached to the number.
  result.stage = inferStage(profile, input.stageOpts)
  log('DONE', `stage: ${result.stage.stage}. ${result.stage.why}`)
  // One contribution or two, per chromosome. Reported, never used to admit or reject a
  // sample: the boundaries are measured for a BULK reference parent, and against a
  // single-cell reference they do not separate at all (audit section E2).
  for (const c of result.chroms) {
    const two = mat !== null && mat !== undefined
    const tally = two
      ? obligateByChrom.get(c.chrom) ?? emptyHet()
      : hetByChrom.get(c.chrom) ?? emptyHet()
    const call = hetCall(tally, two ? 2 : 1)
    if (call.ploidy !== 'uncalled') c.contribution = call
  }
  const uni = result.chroms.filter((c) => c.contribution?.ploidy === 'uniparental')
  if (uni.length) {
    log('WARN', `one parental contribution only on ${uni.map((c) => `chr${c.chrom}`)
      .join(', ')}: ${uni[0].contribution!.why}`
      + (uni[0].contribution!.provisional
        ? '. One parent loaded, so this boundary is provisional' : ''))
  }
  // Segments, after the per-chromosome verdicts, because a chromosome whose calls are not
  // measuring it must not be scanned either: the same broken calls would produce a
  // confident segment inside it.
  const measured = new Set(result.chroms
    .filter((c) => c.verdict !== 'not_measured').map((c) => c.chrom))
  // Copy number first: a chromosome already called whole is not scanned for a segment
  // inside it, since that is the same event described twice.
  const whole = new Set(result.chroms.filter((c) => c.aneuploidy).map((c) => c.chrom))
  const noCall = new Map<string, [number, number]>()
  for (const [c, ms] of cnByChrom) {
    noCall.set(c, [ms.length, ms.filter((m) => !m.called).length])
  }
  const lrrAll = [...cnByChrom.values()].flat()
    .map((m) => m.log2R).filter((x): x is number => x !== null).sort((a, b) => a - b)
  const genomeLrr = lrrAll.length ? lrrAll[lrrAll.length >> 1] : 0
  const copy = [...cnByChrom].filter(([c]) => !whole.has(c))
    .flatMap(([c, ms]) => scanCopyNumber(ms, externalNull(noCall, c), genomeLrr))
  result.segments = [
    ...copy,
    ...[...absenceByChrom].filter(([c]) => measured.has(c))
      .flatMap(([c, ms]) => scanChromosome(ms, externalNull(t.byChrom, c))),
  ].sort((a, b) => b.score - a.score)

  // THE TAXONOMY'S CLASSES, run on the same pass and entering the same defect list.
  //
  // These need no parent at all, which is why they can run before any parental channel is
  // consulted: copy-neutral loss of heterozygosity, runs of homozygosity, ploidy and a
  // genome too disturbed to reference against are all properties of the sample alone.
  // Their PARENTAL ORIGIN, where one exists, is then scored through exactly the same
  // posterior and bands as every older event, because a reader compares rows.
  {
    // Both of the passes whose cost is the size of the array live in `scan.ts`, where a
    // check can reach them and assert they stay linear. The version here rescanned every
    // marker of a chromosome inside its own window loop, which locked the tab for the
    // whole genome on every run. See the note at the top of that file.
    const { windows, scanned, chromEnd } = await copyNeutralWindows(
      cnByChrom, selfMarkers, genomeLrr, LOH_SEGMENT_MARKERS,
      // Hand the page back between chromosomes.
      async () => { await breathe() },
    )
    log('SCAN', `copy-neutral scan: ${windows.length} windows over ${cnByChrom.size} `
      + `chromosomes (${scanned} sliding, ${cnByChrom.size} whole-chromosome)`)

    // AN ARRAY THAT IS NOT MEASURING A GENOME PRODUCES NO FINDINGS.
    //
    // Every detector below reads a genome's own statistics against itself, and that is only
    // meaningful where the array is reading a genome at all. The stage inference already
    // decides this and says so in those words: one example array reads 31.2% heterozygous
    // where a single diploid tops out near 25% on this platform, and its verdict is "no
    // genome reads this way, so this is not a stage". It went on to yield 22 findings.
    //
    // Detecting structure in an array that failed its own quality inference is not a
    // conservative reading of weak data, it is reading noise as biology. The event list
    // stays empty and the reason is the stage's own sentence.
    if (result.stage?.stage === 'failed') {
      log('WARN', 'no chromosomal changes are looked for: this array did not resolve to a '
        + `stage. ${result.stage.why}`)
      result.findings = []
    } else {
    log('SCAN', 'looking for copy-neutral loss of heterozygosity')
    await breathe()
    const findings = [
      // Overlapping windows report the same event several times, and a whole chromosome
      // reports it again, so the redundancy is collapsed: the widest interval covering a
      // position wins, and a segment that is merely its chromosome restated is dropped.
      // A GENOME WITH ONE PARENTAL CONTRIBUTION HAS NO HETEROZYGOSITY TO LOSE, so the
      // copy-neutral detector is not run on it at all. Its relative-depletion test divides
      // by the array's own mean heterozygosity, which on such a genome is near zero.
      ...mergeLoh(detectLoh(windows, { zygosity: result.zygosity })),
      // Same guard as the copy-neutral detector beside it: a genome with one parental
      // contribution is homozygous by construction, so its runs are that one call
      // restated, not separate events.
      ...detectUpd(runsOfHomozygosity(selfMarkers, { chromEndBp: chromEnd }),
        { zygosity: result.zygosity }),
    ]
    log('SCAN', `runs of homozygosity and ploidy over ${selfMarkers.length} called markers`)
    await breathe()
    const tri = detectTriploidy([...myBaf.values()])
    if (tri) findings.push(tri)
    // A genome with too little undisturbed remainder cannot self-reference, which blocks
    // every origin call on the array rather than only on the affected chromosomes.
    const deviant = new Set([...whole, ...result.segments.map((sg) => sg.chrom)]).size
    // Genome call rate from the markers themselves. StageCall does not carry one, and the
    // stage's own inference took it as an input rather than storing it.
    const allMarkers = [...cnByChrom.values()].flat()
    const callRate = allMarkers.length
      ? allMarkers.filter((m) => m.called).length / allMarkers.length : NaN
    const cx = detectComplex(deviant, cnByChrom.size, callRate)
    if (cx) findings.push(cx)
    log('SCAN', `taxonomy: ${findings.length} finding`
      + `${findings.length === 1 ? '' : 's'} across `
      + `${new Set(findings.map((f) => f.chrom)).size} chromosome(s)`)
    await breathe()
    result.findings = findings
    }
    // Which walls stand depends on what the user actually supplied, so it is recorded here
    // rather than guessed at display time. Reading these as undefined would tell a
    // two-parent run that heterodisomy is unreachable when that run has already cleared it.
    result.twoParents = !!mat
    result.units = 1
    // ONE LINE PER KIND OF FINDING, not per finding. A single run produced 183
    // copy-neutral regions whose lines differed only in their measurements, and reading
    // the 183rd told you nothing the first had not. The first of each kind is printed in
    // full, then one line saying how many more there were and where. Every finding keeps
    // its own numbers in the genome viewer, the report and the export, which is where a
    // reader goes to compare them; the log is a narrative of the run.
    const byShape = new Map<string, {
      cls: string; chrom: string; evidence: string; blocked: boolean; at: string[]
    }>()
    for (const f of (result.findings ?? [])) {
      // Numbers masked, so findings differing only in their measurements group together.
      const shape = `${f.cls}|${f.evidence.replace(/[\d.]+/g, '#')}`
      const cur = byShape.get(shape)
      if (cur) cur.at.push(`chr${f.chrom}`)
      else {
        byShape.set(shape, {
          cls: f.cls, chrom: f.chrom, evidence: f.evidence,
          blocked: !!f.originBlocked, at: [`chr${f.chrom}`],
        })
      }
    }
    for (const g of byShape.values()) {
      const tag = g.blocked ? 'WARN' : 'DONE'
      log(tag, `${g.cls} ${g.chrom}: ${g.evidence}`)
      if (g.at.length > 1) {
        log(tag, `and ${g.at.length - 1} more ${g.cls} finding`
          + `${g.at.length === 2 ? '' : 's'} reading the same way, on ${listChroms(g.at)}`)
      }
    }
  }
  // --- where each extra copy came from ---------------------------------------------
  //
  // Two different questions and the tool must not confuse them. With both parents loaded
  // the cell can be biparental and allele dosage says WHICH parent supplied the extra
  // copy. With one parent the cell is uniparental by construction, so the parent is not
  // in question and the answer worth having is whether the extra copy is that parent's
  // other homologue, which is meiotic, or a duplicate of the same one, which is invisible.
  const allDosage = [...dosageByChrom.values()].flat()
  // The sample's own centre, never the theoretical 0.5: against a reconstructed parent
  // the theoretical value is biased toward maternal by up to 0.077, which is enough to
  // invert a call outright.
  const centre = allDosage.length ? recentre(allDosage) : NaN
  // Background per chromosome, excluding the one being judged. A genome-wide rate is
  // lifted by the very events under test where several chromosomes carry one.
  const annotate = (
    where: string, kind: 'whole chromosome' | 'segment', chrom: string,
    from: number, to: number, event: 'gain' | 'loss' = 'gain',
  ) => {
    if (mat && allDosage.length) {
      const inside = (dosageByChrom.get(chrom) ?? [])
        .filter((d) => d.pos >= from && d.pos <= to)
      // A loss reads the same share in the opposite direction, so the two must not share
      // one call: the under-represented parent is the one LOST, not the one that gained.
      const c = event === 'loss'
        ? callLossOrigin(inside, centre, MIN_INFORMATIVE_DEFAULT)
        : callGainOrigin(inside, centre, MIN_INFORMATIVE_DEFAULT)
      return {
        where,
        kind,
        origin: c.origin,
        why: c.why,
        called: c.origin !== 'unclear',
        confidence: c.confidence,
        band: c.band,
      }
    }
    const h = hetByChrom.get(chrom) ?? { informative: 0, het: 0 }
    const background = externalHetBackground(hetByChrom, chrom)
    const c = callHomologue(
      h.informative ? h.het / h.informative : NaN, h.informative, background,
    )
    return {
      where,
      kind,
      origin: c.verdict === 'other homologue'
        ? 'this parent, other homologue (meiotic)' : 'not determinable',
      why: mat ? c.why
        : `${c.why}. Which PARENT supplied it is not in question: this cell carries one `
          + 'parent\'s genome, so the extra copy is that parent\'s. Naming a parent for a '
          + 'gain needs a biparental cell and both parents loaded.',
      called: c.verdict === 'other homologue',
    }
  }
  result.gains = [
    ...result.chroms.filter((c) => c.aneuploidy === 'gain')
      .map((c) => annotate(`chr${c.chrom}`, 'whole chromosome', c.chrom, 0, Infinity)),
    ...result.segments.filter((sg) => sg.kind === 'copy-gain')
      .map((sg) => annotate(
        `chr${sg.chrom} ${(sg.startBp / 1e6).toFixed(1)}-${(sg.endBp / 1e6).toFixed(1)}Mb`,
        'segment', sg.chrom, sg.startBp, sg.endBp,
      )),
  ]
  // The same for losses, which is the direction the whole question usually turns on:
  // a lost chromosome or segment belonged to one parent, and which one is the answer.
  result.losses = [
    ...result.chroms.filter((c) => c.aneuploidy === 'loss')
      .map((c) => annotate(`chr${c.chrom}`, 'whole chromosome', c.chrom, 0, Infinity,
        'loss')),
    ...result.segments.filter((sg) => sg.kind === 'copy-loss')
      .map((sg) => annotate(
        `chr${sg.chrom} ${(sg.startBp / 1e6).toFixed(1)}-${(sg.endBp / 1e6).toFixed(1)}Mb`,
        'segment', sg.chrom, sg.startBp, sg.endBp, 'loss',
      )),
  ]
  // THE FEATURE COMPARISON NO LONGER RUNS HERE, and that is deliberate rather than a
  // regression. It answers a different question from the rest of this run: everything else
  // asks WHOSE a change is, while that asks whether the change sits where the genome breaks
  // anyway. Those are not competing answers, but a reader who meets them in one
  // undifferentiated report treats the second as evidence about the first, and it is not:
  // the fragile compartment is established on BOTH parental genomes from the first cell
  // cycle. It is now an addon, run on request, with its own report.
  //
  // The marker positions it needs are kept here so the addon costs nothing to run later:
  // recomputing them would mean streaming every array a second time.
  {
    const byChrom = new Map<string, number[]>()
    for (const [c, ms] of absenceByChrom) {
      byChrom.set(c, ms.map((m) => m.pos).sort((a, b) => a - b))
    }
    result.markerPositions = byChrom
  }
  // ONE-PARENT ORIGIN. With a single parent loaded, a marker where that parent is
  // homozygous and the sample carries the allele it does NOT have is Mendelian evidence
  // that the parent's copy is absent: dropout removes alleles and never invents one.
  // Validated on a real CEPH trio with the mother hidden, 12 of 12 correct across dropout
  // rates from 0.05 to 0.45, with 5,130-5,172 exclusive markers when the loaded parent's
  // copy was removed and exactly 0 when the other parent's was.
  // WHOLE-CHROMOSOME ORIGIN, FROM DOSAGE. The genotype channel below cannot answer these
  // and the reason is structural: a whole chromosome is DETECTED by the collapse of its
  // genotype call rate, so on exactly those events its evidence is already gone. Dosage is
  // read whether or not a genotype is emitted. Measured on a public series, all four
  // segmental losses scored from genotypes and all three whole-chromosome losses refused.
  // ANYTHING TO SCORE, not only whole chromosomes. This guard used to read `whole.size`
  // alone, so the taxonomy findings nested inside it were scored only when the sample also
  // happened to carry a whole-chromosome aneuploidy. A run whose changes were all
  // copy-neutral events or runs of homozygosity got no origin on any of them, which is
  // most of what the taxonomy detects.
  if (!mat && (whole.size || (result.findings?.length ?? 0) > 0)) {
    // Background from everything OUTSIDE the chromosome under test, which is the same
    // external-null rule the region scan uses: a chromosome cannot set its own baseline.
    // The BAF spread this block used to compute for the array gate now lives on the
    // profile, measured in the streaming pass, and is shown in the stage panel and the
    // quality table. Recomputing it here served only the gate that no longer exists.

    const material = materialOf(result.stage?.stage ?? 'unknown')

    // ONE SCORER, USED TWICE. The whole-chromosome calls here and the taxonomy findings
    // below are the same measurement over different intervals, so the interval is a
    // predicate rather than a chromosome name. Duplicating it for the new classes would
    // have let the two drift apart, and a reader compares their confidences directly.
    // ONE INDEX FOR THE SAMPLE, not one walk of the array per finding. See scan.ts.
    const scanIndex = buildScanIndex(
      { markerPos, parentGt: pat.gt, myBaf, myGt, cnByChrom })
    const scoreInterval = (
      label: string,
      iv: Interval,
      wholeChromosome: boolean,
      state: 'loss' | 'gain' | 'cnn-loh' = 'loss',
      /** The taxonomy class this came from, so the log can group by it. */
      cls = 'dosage',
    ) => {
      // Background is every OTHER chromosome of this same array. Self-referencing is not
      // optional: the raw one-parent null sits at -0.031 on trophectoderm under no event,
      // pointing at the parent that was NOT genotyped, which is the shift a real mosaic
      // fraction of 0.117 would produce.
      // ASK WHETHER ANY ARRAY OF THIS KIND COULD ANSWER BEFORE READING THE ARRAY.
      // callDosageOrigin asks this first too, and says so, but by the time it is called
      // the whole background has already been assembled. The answer depends only on the
      // material, the class, the width and how many parents are loaded, so it is a table
      // lookup, and on single-cell material every segment is undefined at every fraction.
      // Gathering a background to be told that is the bulk of a run on that material.
      const unreachable = originUnreachable(
        material, state, wholeChromosome, mat ? 2 : 1)
      const { region, background, inL, outL, untRows } =
        gatherInterval(scanIndex, iv, { regionOnly: unreachable })
      const mean = (xs: number[]) => (xs.length
        ? xs.reduce((a, x) => a + x, 0) / xs.length : NaN)
      const sdOf = (xs: number[], mu: number) => (xs.length > 1
        ? Math.sqrt(xs.reduce((a, x) => a + (x - mu) ** 2, 0) / (xs.length - 1)) : NaN)
      const muOut = mean(outL)
      const lrrShift = mean(inL) - muOut
      const lrrSe = sdOf(outL, muOut) / Math.sqrt(Math.max(1, inL.length))
      const intensityZ = Number.isFinite(lrrShift) && lrrSe > 0
        ? lrrShift / lrrSe : undefined

      // Spread of the per-window log2R on this chromosome, which decides whether the
      // CLASS can be separated from its nearest feasible alternative. Almost never on
      // amplified material, which is why the origin is emitted without it rather than
      // withheld along with it.
      const inSorted = [...inL].sort((a, b) => a - b)
      const q1 = inSorted[Math.floor(inSorted.length * 0.25)] ?? NaN
      const q3 = inSorted[Math.floor(inSorted.length * 0.75)] ?? NaN
      const windowLogRSd = Number.isFinite(q3 - q1) ? (q3 - q1) / 1.349 : undefined
      const unt = untransmittedPairs(untRows as never)
      const untOriented = unt.pairs.map(orientUntransmitted)
      const untMean = untOriented.length
        ? untOriented.reduce((a, x) => a + x, 0) / untOriented.length : NaN
      // Only once a gain is established: on a euploid chromosome this answers SPH by
      // exclusion every time, which is confidently wrong.
      // The mechanism question only exists once copy number three is established, and
      // that is a whole-chromosome property. Asking it of a sub-chromosomal interval would
      // return "not both homologues" by exclusion on every one of them, which is the trap
      // this gate was added for in the first place.
      const mech = callMechanism(unt.pairs as never,
        { copyNumberThree: wholeChromosome && result.chroms.some(
          (x: { chrom: string, aneuploidy?: string }) =>
            x.chrom === iv.chrom && x.aneuploidy === 'gain',
        ) })

      const c = callDosageOrigin(region as never, background as never, material, {
        wholeChromosome,
        // The array-level refusal is now the structural one: a genome with no undisturbed
        // remainder cannot self-reference. The BAF spread is still measured and reported,
        // it just no longer refuses, because noise already reaches the answer through the
        // standard error and earns a lower band rather than a silence.
        noSelfReference: !!result.findings?.some((f) => f.cls === 'complex'),
        intensityZ,
        parents: mat ? 2 : 1,
        windowLogRSd,
        // The class decides which floor applies, and the floors differ enormously: a
        // copy-neutral event on single-cell material has one at 0.399 where a loss has
        // 0.348 and a segment has none at any fraction.
        state,
      })
      // WHERE THE DOSAGE CHANNEL CANNOT REACH, ASK THE ZYGOSITY. A uniparental sample has
      // one parental genome in it, so every change in it belongs to that parent by
      // construction and needs no detection floor. Only consulted when dosage returned
      // nothing, so a measured answer always wins over an inherited one.
      // Gated on "did not name a parent", not on one refusal verdict: see namedAParent.
      const namedByDosage = namedAParent(c.verdict)
      const zyg = !namedByDosage
        ? uniparentalOrigin({
          originClass: result.originClass, zygosity: result.zygosity, role: soloRole,
          genomeRate: result.genomeRate, explainable: result.explainable,
          hetBand: result.hetBand,
        })
        : null
      // EVERY ROW LEAVES WITH A GRADE. Where neither the measured channel nor the
      // zygosity can answer, the row is graded F rather than dropped: the direction the
      // displacement leans is named, and the grade says there is nothing behind it. A
      // reader scanning a column of grades sees F and knows the row is unusable, which an
      // empty field never conveyed as directly. An F is not a result and nothing should
      // aggregate or act on one.
      // THE LAST RUNG, for an interval that carried no usable marker of its own. There is
      // no direction to measure there, so the only honest source left is the array-wide
      // one: absence of the loaded parent measured across the whole genome. It says
      // nothing about THIS interval and the reason says so. Graded F like the rest.
      const genomeLean: 'loaded-parent' | 'other-parent' | null =
        Number.isFinite(result.genomeRate) && Number.isFinite(result.explainable)
          ? (result.genomeRate > result.explainable ? 'other-parent' : 'loaded-parent')
          : null
      const unnamed = !namedByDosage
      // MEASURED AT 0.27 TO 0.44, so this grade names no parent. See the note on
      // BAND_F_MIN: every call landing here has an unresolved class, and with the class
      // unresolved a gain inverts the sign that loss and copy-neutral share, so the
      // direction of the shift carries no information about which parent it was.
      const graded = !zyg && unnamed && !c.lean && genomeLean
        ? {
          verdict: 'not-evaluable' as const,
          confidence: undefined,
          band: 'F' as const,
          why: `${c.why}. Graded F: this interval carried no usable marker of its own, `
            + 'so nothing here bears on which parent it was. The parent is NOT named. An '
            + 'injection series on real arrays recovers the parent 0.51 to 0.56 of the '
            + 'time in this band, at chance, because every call reaching it has an '
            + 'unresolved copy-number class and a gain inverts the sign that loss and '
            + 'copy-neutral share',
        }
        : !zyg && unnamed && c.lean
        ? {
          verdict: 'not-evaluable' as const,
          confidence: undefined,
          band: 'F' as const,
          why: `${c.why}. Graded F: the displacement over ${c.lean.markers} usable `
            + `marker${c.lean.markers === 1 ? '' : 's'} does not determine a parent. The `
            + 'copy-number class is unresolved here, and a gain inverts the sign that loss '
            + 'and copy-neutral share, so the direction of the shift carries no parental '
            + 'information at all. Measured on real arrays this band recovers the parent '
            + '0.51 to 0.56 of the time, which is chance, so no parent is named',
        }
        : null
      return {
        where: label,
        verdict: zyg?.verdict ?? graded?.verdict ?? c.verdict,
        // NO NUMBER FROM THE ZYGOSITY CHANNEL. It emits a verdict and the margin behind
        // it; a confidence there would be a reparameterisation of the margin wearing three
        // decimal places. See the note on UniparentalCall.
        confidence: zyg ? undefined : (graded?.confidence ?? c.posterior?.confidence),
        band: zyg?.band ?? graded?.band ?? c.posterior?.band,
        inheritedMargin: zyg?.foldOverCeiling,
        fromZygosity: !!zyg,
        gradedOnly: !!graded,
        parent: zyg?.parent,
        limitedBy: c.posterior?.limitedBy,
        uncalibrated: c.posterior?.uncalibrated,
        classVerdict: c.classVerdict,
        classWhy: c.classWhy,
        untransmittedMarkers: unt.pairs.length,
        untransmittedAmbiguous: unt.ambiguous,
        untransmittedShare: untMean,
        untransmittedImpossible: impossibleRate(unt.pairs as never),
        mechanism: mech.mechanism,
        mechanismWhy: mech.why,
        shift: c.shift,
        z: c.z,
        impliedF: c.impliedF,
        window: c.window,
        markers: c.markers,
        material: c.material,
        floor: c.floor,
        why: zyg ? zyg.why : (graded?.why ?? c.why),
        cls,
      }
    }

    if (whole.size) {
      log('SCAN', `scoring ${whole.size} whole-chromosome event`
        + `${whole.size === 1 ? '' : 's'} from allele dosage`)
    }
    result.dosageCalls = []
    for (const chrom of whole) {
      result.dosageCalls.push(scoreInterval(`chr${chrom}`, { chrom }, true))
      await breathe()
    }

    // THE TAXONOMY'S FINDINGS GO THROUGH THE SAME SCORER. A copy-neutral event and an
    // isodisomy carry a parental origin exactly as a deletion does, and it must be the
    // same posterior, the same four bands and the same class-inversion veto, or a reader
    // comparing two rows of one table is comparing two different kinds of number.
    //
    // Classes whose origin is blocked BY THE CLASS are skipped rather than scored and
    // discarded: a triploidy has no parental origin at any quality, and running the
    // statistic on it would produce a number that means nothing.
    const scorable = (result.findings ?? [])
      .filter((f) => !originBlockedByClass(f.cls) && f.chrom !== 'genome')
    if (scorable.length) {
      log('SCAN', `scoring parental origin on ${scorable.length} finding`
        + `${scorable.length === 1 ? '' : 's'}`)
    }
    let scoredSoFar = 0
    for (const f of scorable) {
      // THE FINDING'S OWN CLASS DECIDES ITS FLOOR, and passing the default instead cost
      // most of them an answer. A copy-neutral event is the LARGEST-signal class, floor
      // 0.399 on single-cell material against a loss's 0.348 and a segment's none at all.
      const call = scoreInterval(
        `chr${f.chrom} ${(f.startBp / 1e6).toFixed(1)}-${(f.endBp / 1e6).toFixed(1)}Mb`,
        { chrom: f.chrom, startBp: f.startBp, endBp: f.endBp },
        f.wholeChromosome,
        f.cls === 'cnn-loh' ? 'cnn-loh'
          : f.cls === 'segmental-duplication' ? 'gain' : 'loss',
        f.cls,
      )
      result.dosageCalls.push(call)
      scoredSoFar += 1
      // One yield every couple of findings, so the lines above appear as they are decided
      // rather than all at once when the loop ends.
      if (scoredSoFar % 2 === 0) await breathe()
    }
    // A NAMED PARENT IS ALWAYS PRINTED IN FULL. It is the answer, and two of them that
    // happen to read alike are still two answers.
    //
    // A REFUSAL IS PRINTED ONCE PER DISTINCT REASON. A run on single-cell material
    // produced 319 of these and every one was the same sentence: no array of this kind at
    // this width can answer, whatever the data says. Printing it once per interval told a
    // reader nothing the first one had not, buried the handful of real answers, and cost a
    // render each. The intervals are still listed, and the full table is still in the
    // report and the export.
    const refusals = new Map<string, { cls: string; why: string; at: string[] }>()
    for (const c of result.dosageCalls) {
      // 'refused' is the GENOTYPE channel's vocabulary and never appears here, so this
      // line logged DONE for every outcome including a withheld parent.
      const named = c.verdict === 'loaded-parent' || c.verdict === 'other-parent'
      const cls = (c as { cls?: string }).cls ?? 'dosage'
      if (named) { log('DONE', `${cls} origin ${c.where}: ${c.why}`); continue }
      // Grouped by CLASS as well as reason, so the summary says what it stood for rather
      // than lumping copy-neutral events together with uniparental disomy.
      const key = `${cls}|${c.why}`
      const at = refusals.get(key)
      if (at) at.at.push(c.where)
      else refusals.set(key, { cls, why: c.why, at: [c.where] })
    }
    for (const { cls, why, at } of refusals.values()) {
      if (at.length === 1) { log('WARN', `${cls} origin ${at[0]}: ${why}`); continue }
      log('WARN', `${at.length} ${cls} intervals not evaluable for the same reason, on `
        + `${listChroms(at)}: ${why}`)
    }
  }
  // GATED ON THERE BEING AN EVENT, NOT ON THERE BEING A SEGMENT.
  //
  // This block scores whole chromosomes and segments alike, but the guard asked only about
  // segments, so a sample whose only changes were whole-chromosome losses skipped it entirely and
  // those events fell through to the dosage channel, which cannot reach them on this material.
  // The clearest event the tool sees was the one guaranteed to miss its best channel.
  //
  // AND NOT RUN AT ALL ON A UNIPARENTAL GENOME. This channel reads "the loaded parent's copy is
  // absent here", which on a genome carrying one parental contribution is true on every
  // chromosome by construction and says nothing about any interval. Worse, it names the wrong
  // parent: on a gynogenetic genome the copy that went missing from a chromosome was the maternal
  // one, because that was the only copy there. Same guard, and same reason, as the copy-neutral
  // and runs-of-homozygosity detectors above.
  const uniparental = result.zygosity?.startsWith('uniparental') ?? false
  const mendelEvents = !mat && !uniparental ? [
    ...(result.chroms ?? []).filter((c) => c.aneuploidy)
      .map((c) => ({ chrom: c.chrom, start: 0, end: Number.MAX_SAFE_INTEGER,
        label: `chr${c.chrom}` })),
    ...result.segments.map((sg) => {
      const co = segmentCoords(sg)
      return { chrom: sg.chrom, start: co.start, end: co.end,
        label: `chr${sg.chrom} ${(co.start / 1e6).toFixed(1)}-${(co.end / 1e6).toFixed(1)}Mb` }
    }),
  ] : []
  if (mendelEvents.length) {
    // THE OBVIOUS EVENTS GET THE BEST CHANNEL, which is where they were not getting it.
    //
    // A whole chromosome gained or lost is the clearest event this tool sees, and until now
    // only the dosage channel was asked about its parent, which on amplified material
    // returns band D or F. The Mendelian channel was run on SEGMENTS alone.
    //
    // On a biparental sample that channel does not need a detection floor at all: at a
    // marker where the loaded parent is homozygous, losing that parent's copy leaves an
    // allele the parent does not have, and dropout removes alleles without inventing one.
    // Measured by removing one parent's copy across a chromosome on real biparental arrays,
    // both directions from the same array: 92 of 100 calls correct, per-array 0.920 +/-
    // 0.100, on arrays that resolve to a stage. The same experiment on arrays their own
    // inference rejects returns 0.650, which is why those are excluded rather than
    // reported weakly.
    result.oneParent = mendelEvents.map((sg) => {
      const co = { start: sg.start, end: sg.end }
      const pairs: [string, string][] = []
      for (const [probe, gt] of myGt) {
        const q = markerPos.get(probe)
        if (!q || q.chrom !== sg.chrom || q.pos < co.start || q.pos > co.end) continue
        const pg = pat.gt.get(probe)
        if (pg) pairs.push([pg, gt])
      }
      // Dropout comes from the inferred stage. A failed array carries no usable figure
      // rather than a flattering one, so it gets the most conservative dropout measured
      // on any stage instead of NaN, which would silently void every likelihood below.
      const ado = Number.isFinite(result.stage!.dropout) ? result.stage!.dropout : 0.308
      const c = callOneParentOrigin(pairs as never, ado)
      return {
        where: sg.label,
        verdict: c.verdict,
        posterior: c.posterior,
        band: c.band,
        markers: c.markers,
        exclusive: c.exclusive,
        why: c.why,
      }
    })
    for (const c of result.oneParent) {
      log(c.verdict === 'refused' ? 'WARN' : 'DONE', `origin ${c.where}: ${c.why}`)
    }
  }
  // Sibling-referenced call, where the run holds other cells of the same embryo. Reports
  // whether a copy is genuinely missing rather than which parent's: the observations here
  // are UNPHASED, so a side call would be counting a per-marker reference allele across
  // markers where it denotes different parental sides, which cancels. Naming needs phase
  // and then one anchor per donor group; the module refuses rather than guessing.
  // Siblings are supplied by the caller: which arrays are units of one embryo is a property of
    // the run, not of this sample.
  if (sibs.length >= 2 && result.segments.length) {
    const need = hetRule(sibs.length)
    result.siblingCalls = result.segments.map((sg) => {
      const co = segmentCoords(sg)
      const obs: SibAB[] = []
      for (const [probe, gt] of myGt) {
        const pos = markerPos.get(probe)
        if (!pos || pos.chrom !== sg.chrom || pos.pos < co.start || pos.pos > co.end) continue
        let het = 0
        for (const sib of sibs) if (sib.get(probe) === 'AB') het += 1
        if (het >= need) obs.push(gt as SibAB)
      }
      const call = callSiblingOrigin(obs, sibs.length, profile.nocallRate)
      return {
        where: `chr${sg.chrom} ${(co.start / 1e6).toFixed(1)}-${(co.end / 1e6).toFixed(1)}Mb`,
        hypothesis: call.hypothesis,
        posterior: call.posterior,
        markers: call.markers,
        phi: call.phi,
        why: call.why,
      }
    })
    for (const c of result.siblingCalls) {
      log(c.hypothesis === 'refused' ? 'WARN' : 'DONE', `sibling call ${c.where}: ${c.why}`)
    }
  }
  for (const g of result.gains) {
    log(g.called ? 'DONE' : 'WARN', `gain ${g.where}: ${g.origin}. ${g.why}`)
  }
  for (const l of result.losses) {
    log(l.called ? 'DONE' : 'WARN', `loss ${l.where}: ${l.origin} copy missing. ${l.why}`)
  }

  if (result.segments.length) {
    log('WARN', `${sampleName}: ${result.segments.length} segment(s) where the `
      + `${soloRole} `
      + `genome is missing: ${result.segments.map((x) => `chr${x.chrom} `
        + `${(segmentCoords(x).spanBp / 1e6).toFixed(1)}Mb `
        + `${segmentCoords(x).localised ? segmentCoords(x).interval : '(not localised)'} `
        + `at ${pct(x.rate, 1)}`).join(', ')}`)
  }

  return result
}
