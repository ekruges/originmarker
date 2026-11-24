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
  classify, emptyTally, isAutosome, knownBiparental, pct, tallyRow, PANEL_HET_FLOOR,
  type ParentageResult,
} from './parentage.ts'
import {
  scanChromosome, scanCopyNumber, scanCopyGain, highLrrProfile, externalNull,
  segmentCoords, type MarkerAbsence,
  type Segment, RECIPROCAL_OVERLAP, MIN_SEGMENT_SPAN_BP, canDetectLoss,
} from './segments.ts'
import { readLostChromosome, RESTATED_SHARE, type PartLost } from './partLost.ts'
import { armSegment, readArms, type ArmReading } from './wholeArms.ts'
import {
  paternalShare, recentre, callGainOrigin, callHomologue, externalHetBackground, scanAlleleGains,
  type ChromShare,
  MIN_INFORMATIVE_DEFAULT, type DosageMarker, callLossOrigin, ONE_PARENT_ORIENT_Z,
} from './gainOrigin.ts'
import {
  addTwoParent, emptyHet, hetCall, type HetTally, MIXED_SHARE, HAPLOID_MAX, nameAbsentParent,
  lossSignalRate, whyNotInformative,
  DIPLOID_MIN,
} from './obligateHet.ts'
import { originBlockedByClass, namedAParent } from './defects.ts'
import { changeLogLines } from './pageSays.ts'
import { callSiblingOrigin, hetRule, type AB as SibAB } from './siblingOrigin.ts'
import { callOneParentOrigin } from './oneParentOrigin.ts'
import { callDosageOrigin, materialOf, originUnreachable } from './dosageOrigin.ts'
import { uniparentalOrigin } from './uniparentalOrigin.ts'
import { accumulateSex, emptySex, sexCall, type SexCall, type SexTally } from './sexing.ts'
import type { InferredBanner } from './inferredArray.ts'
import { callSexChromosomes } from './sexChromosome.ts'
import { sexMaterial } from './sexMaterial.ts'
import { PROVEN_CHIP_GENOTYPES_Y } from './platform.ts'
import { reconcileParentSex } from './parentSanity.ts'
import {
  relate, OPPOSITE_HOM_MAX, relatednessAssessable, REL,
} from './relatedness.ts'
import { callBothParentsOrigin } from './bothParentsOrigin.ts'
import {
  otherRates, scanOtherParent, ownRates, scanOwnCopy, coveredShare, sameGenomeShare, presentScore,
  OTHER_MIN_SPAN_BP, OTHER_LRT, SAME_GENOME_SHARE, type OtherMarker, type OtherRates,
} from './otherParentLoss.ts'
import {
  detectLoh, detectUpd, detectTriploidy, detectComplex, runsOfHomozygosity, mergeLoh, COMPLEX_CALL_RATE,
  LOH_SEGMENT_MARKERS,
} from './abnormalities.ts'
import { untransmittedPairs, impossibleRate, orientUntransmitted, callMechanism } from './untransmitted.ts'
import {
  emptyHomTally, addHomozygous, informationRatio, whyNoInformation, type HomTally,
  NO_INFORMATION_MIN,
} from './informationRatio.ts'
import {
  emptySameHom, addSameHom, dropIn, MIN_SAME_HOM, onParentScale, type SameHomTally,
  type PanelCalibration, type FounderPanel,
} from './answerKeys.ts'
import type { PhysicalCause } from './refusals.ts'
import {
  emptyShareTally, addShare, shareCall, normalFraction, NORMAL_MAJORITY_MIN, DISOMIC_LOW,
  type ShareTally, type ShareCall,
} from './oneParentState.ts'
import {
  gainReadAsLoss, emptyWeightedShare, addWeightedShare, weightedShareR, hetCallShare,
  otherCopyFromHetCalls, type WeightedShare,
} from './oneParentGain.ts'
import { motherFromRun, routeGains, type RunMarker, type RouteChrom } from './motherFromRun.ts'
import {
  emptyCarryTally, addCarry, parentGenome, parentGenomeSentence, parentGenomeNotice, UNREAD,
  type CarryTally,
} from './parentGenome.ts'
import {
  inferStage, stageDefaults, type Stage, PANEL_ANCHOR_CALL_RATE, parentOverCeiling, MAX_DIPLOID_HET,
  panelScale, hetSpread, HET_SPREAD_LOCALISED, allelePloidy,
} from './stage.ts'
import { nullScale, calibratedZ, median, intensityDetects, STATE_GAP, Z_WINDOW } from './intensityNull.ts'
import { reconcileStage } from './declaredStage.ts'
import { assessIntegrity } from './integrity.ts'

/** One parent's array, reduced to what every channel below reads. */
export interface ParentIndex {
  gt: Map<string, AB>
  /**
   * The parent's heterozygosity over called autosomal markers.
   *
   * ON A RECONSTRUCTED REFERENCE THIS CANNOT BE MEASURED FROM THE FILE, which is homozygous at
   * every marker by construction and therefore reads exactly 0.0%. When `inferred` is set this
   * carries the banner's `hRetained` instead, estimated from product disagreement at build time.
   * `measuredHeterozygosity` keeps what the file itself says, so the substitution is visible.
   */
  heterozygosity: number
  /** Called fraction of this array's autosomal markers. */
  callRate: number
  build: string | null
  /**
   * What chromosome Y says about this array. Read on the PARENT, which is the slot an operator can
   * put the wrong file into: see parentSanity.ts.
   */
  sex: SexCall
  /**
   * This array's banner, when it is a reconstruction rather than a measurement.
   *
   * Present means every figure keyed to the reference's own imperfection is available: the
   * absence floor it contributes and the parental heterozygosity the genotypes cannot show. A
   * caller that leaves this unset on a reconstructed file is scoring against a parent whose
   * heterozygosity reads zero, which withholds the second-parent axis for the whole run.
   */
  inferred?: InferredBanner
  /** Heterozygosity as the FILE reads, before any banner substitution. Equal to
   *  `heterozygosity` on a measured array; 0 on a reconstruction. */
  measuredHeterozygosity: number
  /** No parental array was loaded: see `noParent`. */
  absent?: boolean
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
export interface ParentAccum {
  gt: Map<string, AB>; called: number; het: number
  /** Autosomal markers seen, so the parental array's own call rate is known. */
  markers: number
  /** Chromosome Y, accumulated in the same pass. See parentSanity.ts for what it is for. */
  sex: SexTally
}

export const emptyParent = (): ParentAccum => ({
  gt: new Map(), called: 0, het: 0, markers: 0, sex: emptySex(),
})

/**
 * Autosomes that must read both parents at the opposite-homozygote markers before any chromosome
 * of the same genome is read as missing one. One would do in principle; two keeps a single noisy
 * chromosome from admitting a uniparental genome.
 */
export const BIPARENTAL_AUTOSOMES_MIN = 2

export function collectParentRow(r: ProbeRow, acc: ParentAccum): void {
  acc.gt.set(r.probesetId, r.genotype)
  // ON THE SAME PASS, because the alternative is a second read of a 790,000-row file to answer a
  // question that costs two counters. This is what lets the tool notice that the array in the
  // paternal slot is not male.
  accumulateSex(r, acc.sex)
  if (isAutosome(r.chrom)) {
    acc.markers += 1
    if (r.genotype !== 'NC') {
      acc.called += 1
      if (r.genotype === 'AB') acc.het += 1
    }
  }
}

/**
 * Heterozygosity is autosomal and over CALLED markers: it is the array's own error scale.
 *
 * EXCEPT ON A RECONSTRUCTION, where that quantity is zero by construction and says nothing. Pass
 * the file's banner and the parent's real heterozygosity comes from there instead. The guard in
 * `parentage.classify` that withholds the second-parent axis below PANEL_HET_FLOOR is correct and
 * is left alone: a reconstruction genuinely cannot supply that axis from its genotypes, and the
 * repair is to give the guard the figure the reconstruction DID establish, not to relax it.
 */
export const finishParent = (
  acc: ParentAccum, build: string | null, inferred?: InferredBanner,
): ParentIndex => {
  const measured = acc.called ? acc.het / acc.called : NaN
  return {
    gt: acc.gt,
    // The banner only substitutes when it carries a usable figure; an unreadable banner leaves
    // the measured zero in place, which withholds rather than guesses.
    heterozygosity: inferred && Number.isFinite(inferred.hRetained) && inferred.hRetained > 0
      ? inferred.hRetained : measured,
    measuredHeterozygosity: measured,
    // Whether this array is in a state to say anything about its panel: see PANEL_ANCHOR_CALL_RATE.
    callRate: acc.markers ? acc.called / acc.markers : NaN,
    build,
    sex: sexCall(acc.sex),
    ...(inferred ? { inferred } : {}),
  }
}

/**
 * The parental slot with nothing in it. A run against it reads what needs no parent: the stage, the
 * sex chromosomes, and copy number from the call rate and the intensity. Every tally keyed to a
 * parental genotype stays empty, so no chromosome is given a parent.
 */
export const noParent = (): ParentIndex => ({ ...finishParent(emptyParent(), null), absent: true })

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
  hetByChrom: Map<string, { offered: number; informative: number; het: number }>
  /** The same tally at markers where the loaded parent is HETEROZYGOUS, which the one-parent
   *  contribution call reads instead. See hetCallOneParent. */
  parentHetByChrom: Map<string, { informative: number; het: number }>
  selfMarkers: { chrom: string; pos: number; het: boolean }[]
  obligateByChrom: Map<string, HetTally>
  matAbsenceByChrom: Map<string, MarkerAbsence[]>
  matAbsenceTally: Map<string, [number, number]>
  absenceByChrom: Map<string, MarkerAbsence[]>
  cnByChrom: Map<string, { chrom: string; pos: number; called: boolean; log2R: number | null;
    baf: number | null }[]>
  myGt: Map<string, string>
  myBaf: Map<string, number>
  markerPos: Map<string, { chrom: string; pos: number }>
  /** Chromosome Y of the SAMPLE. What one X means depends on it; see sexChromosome.ts. */
  sex: SexTally
  /**
   * Shared autosomal homozygotes against each MEASURED parent, for the information ratio.
   *
   * Null where that parent is a reconstruction: a reconstructed reference is homozygous at every
   * marker by construction, and one built from this sample's own relatives agrees with it for a
   * reason that is not array quality.
   */
  informationPat: HomTally | null
  informationMat: HomTally | null
  /** Markers where both MEASURED parents are the same homozygote: this cell's answer key for drop-in
   *  and for an absent allele. Null unless both parents are measured arrays. See answerKeys.ts. */
  sameHom: SameHomTally | null
  /**
   * The parent-free answer key. A panel of unrelated founders' shared homozygotes stands in for a
   * pair of parents that always agree: drop-in is counted where the panel says the sample can only
   * be one homozygote, and the information ratio is taken against the panel as the reference. So
   * both answer keys exist with no parent loaded, or with one reconstructed from the children.
   */
  panel: ReadonlyMap<string, 'AA' | 'BB'> | null
  panelSameHom: SameHomTally | null
  informationPanel: HomTally | null
  /** The one-parent copy-state tallies per autosome, where the panel carries founder frequencies.
   *  See oneParentState.ts. */
  shareByChrom: Map<string, ShareTally> | null
  /** The parent-genome tallies per autosome against each loaded parent, where the panel carries
   *  founder frequencies. See parentGenome.ts. */
  carryByChrom: Map<string, CarryTally> | null
  carryMatByChrom: Map<string, CarryTally> | null
  /** Per autosome, the markers the likelihood scans of otherParentLoss.ts read against each loaded
   *  parent. Null where that parent is not loaded or the panel carries no founder frequencies. */
  otherByChrom: Map<string, OtherMarker[]> | null
  otherMatByChrom: Map<string, OtherMarker[]> | null
  /** One parent loaded: the gain test's share with a rebuilt parent's markers weighted (null with a
   *  measured parent), and the share at heterozygous calls. See oneParentGain.ts. */
  shareWByChrom: Map<string, WeightedShare> | null
  hetShareByChrom: Map<string, number[]> | null
  /** Where the run's other draws read the other parent as the loaded one's opposite homozygote
   *  (motherFromRun.ts). Set by the scoring pass for the cell figures; nothing scores on it. */
  otherOpposite?: Set<string> | null
}

/**
 * A fresh accumulator for one sample.
 *
 * The parent's assembly is stamped on the tallies here, so the pseudoautosomal boundaries are the
 * right ones. It is known before the sample streams and is the same across one experiment.
 */
export function emptyCollected(
  pat: ParentIndex, mat?: ParentIndex | null, panel?: ReadonlyMap<string, 'AA' | 'BB'> | null,
): Collected {
  const t = emptyTally()
  const tm = mat ? emptyTally() : null
  t.build = pat.build
  if (pat.absent) t.noParent = true
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
    parentHetByChrom: new Map(),
    // Sample-only zygosity, for the classes that need no parent: LOH, runs, ploidy.
    selfMarkers: [],
    // Parental contribution per chromosome. With both parents this is the obligate set: both
    // homozygous and OPPOSITE, where a biparental cell must read heterozygous. With one parent
    // `hetByChrom` is already the same tally in its one-parent form.
    obligateByChrom: new Map(),
    matAbsenceByChrom: new Map(),
    matAbsenceTally: new Map(),
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
    sex: emptySex(),
    informationPat: pat.inferred ? null : emptyHomTally(),
    informationMat: mat && !mat.inferred ? emptyHomTally() : null,
    sameHom: !pat.inferred && mat && !mat.inferred ? emptySameHom() : null,
    panel: panel ?? null,
    panelSameHom: panel ? emptySameHom() : null,
    informationPanel: panel ? emptyHomTally() : null,
    shareByChrom: (panel as FounderPanel | null | undefined)?.freqB ? new Map() : null,
    carryByChrom: (panel as FounderPanel | null | undefined)?.freqB && !pat.absent ? new Map() : null,
    carryMatByChrom: (panel as FounderPanel | null | undefined)?.freqB && mat ? new Map() : null,
    otherByChrom: (panel as FounderPanel | null | undefined)?.freqB && !pat.absent ? new Map() : null,
    otherMatByChrom: (panel as FounderPanel | null | undefined)?.freqB && mat ? new Map() : null,
    shareWByChrom: (panel as FounderPanel | null | undefined)?.freqB && !mat && pat.inferred
      && Number.isFinite(pat.inferred.mMin) ? new Map() : null,
    hetShareByChrom: (panel as FounderPanel | null | undefined)?.freqB && !mat && !pat.absent ? new Map() : null,
  }
}

/** One called marker into a parent's list for the likelihood scans. See otherParentLoss.OtherMarker. */
const pushOther = (
  by: Map<string, OtherMarker[]>, parent: string, r: ProbeRow, freqB: ReadonlyMap<string, number>,
): void => {
  if (r.genotype === 'NC' || (parent !== 'AA' && parent !== 'AB' && parent !== 'BB')) return
  const fq = parent === 'AB' ? undefined : freqB.get(r.probesetId)
  const q = parent === 'AB' ? null : fq === undefined ? 0 : parent === 'AA' ? fq : 1 - fq
  if (q !== null && !(q > 0)) return
  let ms = by.get(r.chrom)
  if (!ms) by.set(r.chrom, ms = [])
  ms.push({ id: r.probesetId, pos: r.pos, q,
    seen: parent === 'AB' ? r.genotype === 'AB' : r.genotype !== parent,
    opp: parent !== 'AB' && r.genotype !== 'AB' && r.genotype !== parent })
}

const carryTally = (m: Map<string, CarryTally>, chrom: string): CarryTally => {
  let t = m.get(chrom)
  if (!t) m.set(chrom, t = emptyCarryTally())
  return t
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
    t, tm, myGt, myBaf, markerPos, cnByChrom, selfMarkers, hetByChrom,
    obligateByChrom,
    matAbsenceByChrom, matAbsenceTally,
    absenceByChrom, dosageByChrom,
  } = acc
  // On the same pass, for the same reason the parental arrays do it: a second read of an
  // 800,000-row file to count one chromosome costs more than the counters do.
  accumulateSex(r, acc.sex)
  if (isAutosome(r.chrom) && r.genotype !== 'NC') {
    myGt.set(r.probesetId, r.genotype)
    markerPos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
  }
  const fa = pat.gt.get(r.probesetId) ?? 'NC'
  tallyRow(fa, r, t)
  if (tm && mat) tallyRow(mat.gt.get(r.probesetId) ?? 'NC', r, tm)
  // Shared homozygotes, for the information ratio. On this pass because a second read of an
  // 800,000-row file to count four numbers costs more than the counters do.
  if (isAutosome(r.chrom)) {
    if (acc.informationPat) addHomozygous(acc.informationPat, r.genotype, fa)
    if (acc.informationMat && mat) {
      addHomozygous(acc.informationMat, r.genotype, mat.gt.get(r.probesetId) ?? 'NC')
    }
    if (acc.sameHom && mat) addSameHom(acc.sameHom, r.genotype, fa, mat.gt.get(r.probesetId) ?? 'NC')
    const g = acc.panel?.get(r.probesetId)
    if (g) {
      addSameHom(acc.panelSameHom!, r.genotype, g, g)
      addHomozygous(acc.informationPanel!, r.genotype, g)
    }
    // The one-parent copy state reads BAF, at no-calls too, against the founders' frequency.
    if (acc.shareByChrom) {
      const fq = (acc.panel as FounderPanel).freqB!.get(r.probesetId)
      if (fq !== undefined) {
        let st = acc.shareByChrom.get(r.chrom)
        if (!st) acc.shareByChrom.set(r.chrom, st = emptyShareTally())
        addShare(st, fa, r.baf, fq, g === fa)
        if (acc.shareWByChrom) {
          let w = acc.shareWByChrom.get(r.chrom)
          if (!w) acc.shareWByChrom.set(r.chrom, w = emptyWeightedShare())
          addWeightedShare(w, fa, r.baf, fq, pat.inferred?.mMin)
        }
        if (acc.hetShareByChrom) {
          const y = hetCallShare(fa, r.genotype, r.baf, fq)
          if (y !== null) {
            let ys = acc.hetShareByChrom.get(r.chrom)
            if (!ys) acc.hetShareByChrom.set(r.chrom, ys = [])
            ys.push(y)
          }
        }
        if (acc.carryByChrom) addCarry(carryTally(acc.carryByChrom, r.chrom), fa, r.genotype, fq)
        if (acc.carryMatByChrom && mat) {
          addCarry(carryTally(acc.carryMatByChrom, r.chrom), mat.gt.get(r.probesetId), r.genotype, fq)
        }
      }
    }
    if (acc.otherByChrom) pushOther(acc.otherByChrom, fa, r, (acc.panel as FounderPanel).freqB!)
    if (acc.otherMatByChrom && mat) {
      pushOther(acc.otherMatByChrom, mat.gt.get(r.probesetId) ?? 'NC', r, (acc.panel as FounderPanel).freqB!)
    }
  }
  // Only a marker where the parent is homozygous and the sample is called can say
  // anything about presence, which is the same informative set the rates use.
  if (isAutosome(r.chrom)) {
    markerPos.set(r.probesetId, { chrom: r.chrom, pos: r.pos })
    if (r.baf !== null && Number.isFinite(r.baf)) myBaf.set(r.probesetId, r.baf)
    const cn = cnByChrom.get(r.chrom) ?? []
    // The B-allele frequency rides along: the gain arm reads it to ask whether a raised stretch
    // holds two different alleles (parentage.GAIN_OFF_POLE_MIN).
    cn.push({ chrom: r.chrom, pos: r.pos, called: r.genotype !== 'NC', log2R: r.log2R, baf: r.baf })
    cnByChrom.set(r.chrom, cn)
    // The SAMPLE's own zygosity, which is a different tally from hetByChrom above: that
    // one counts markers where the PARENT is homozygous. Copy-neutral LOH and runs of
    // homozygosity are properties of the sample alone and need no parent at all.
    if (r.genotype !== 'NC') {
      selfMarkers.push({ chrom: r.chrom, pos: r.pos, het: r.genotype === 'AB' })
    }
  }
  // THE CLASS THE ONE-PARENT CHANNEL EXCLUDES, which carries four times the signal for the
  // question of whether a second contribution is present. See hetCallOneParent.
  if (fa === 'AB' && r.genotype !== 'NC' && isAutosome(r.chrom)) {
    const ph = acc.parentHetByChrom.get(r.chrom) ?? { informative: 0, het: 0 }
    ph.informative += 1
    if (r.genotype === 'AB') ph.het += 1
    acc.parentHetByChrom.set(r.chrom, ph)
  }
  // MARKERS THE PARENTS OFFER, counted whether or not this sample managed to call them. The
  // homologue block below requires a call, so the offered count has to be taken here as well or a
  // chromosome the array failed on would report that the parents offered nothing. The two numbers
  // together are what separates a marker-poor chromosome from an array that did not read one:
  // `addTwoParent` counts the offer first and returns on a no-call. See whyNotInformative.
  if ((fa === 'AA' || fa === 'BB') && r.genotype === 'NC' && isAutosome(r.chrom)) {
    const h = hetByChrom.get(r.chrom) ?? { offered: 0, informative: 0, het: 0 }
    h.offered += 1
    hetByChrom.set(r.chrom, h)
    if (mat) {
      const ob = obligateByChrom.get(r.chrom) ?? emptyHet()
      addTwoParent(fa, mat.gt.get(r.probesetId) ?? 'NC', 'NC', ob)
      obligateByChrom.set(r.chrom, ob)
    }
  }
  if ((fa === 'AA' || fa === 'BB') && r.genotype !== 'NC' && isAutosome(r.chrom)) {
    // The homologue channel: every heterozygous call here is error in a haploid.
    const h = hetByChrom.get(r.chrom) ?? { offered: 0, informative: 0, het: 0 }
    h.offered += 1
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
  // THE SAME QUESTION ASKED OF THE OTHER PARENT. The series above is gated on the FIRST parent
  // being homozygous and asks whether that parent's allele is missing, so a deletion of the second
  // parent's copy leaves it reading present at every marker and is never scanned for. On the
  // hand-karyotyped cells that asymmetry showed as 11 of 17 paternal segmental events found and 0
  // of 10 maternal. The mother is homozygous at a different set of markers than the father, so
  // this cannot sit inside the gate above and carries its own null.
  if (mat && isAutosome(r.chrom) && r.genotype !== 'NC') {
    const mo = mat.gt.get(r.probesetId) ?? 'NC'
    if (mo === 'AA' || mo === 'BB') {
      const ma = matAbsenceByChrom.get(r.chrom) ?? []
      const absent = r.genotype !== 'AB' && r.genotype !== mo
      ma.push({ chrom: r.chrom, pos: r.pos, absent })
      matAbsenceByChrom.set(r.chrom, ma)
      const d = matAbsenceTally.get(r.chrom) ?? [0, 0]
      d[0] += 1
      if (absent) d[1] += 1
      matAbsenceTally.set(r.chrom, d)
    }
  }
}

/**
 * motherFromRun.ts wired in: the run's other arrays read for the other parent's genotype, and the
 * two-parent test run on it where the first build named no parent. Names gains only; every other verdict
 * stands. The limits say how many draws of the other parent were used, or why none were.
 */
function readMotherFromRun(x: {
  acc: Collected; pat: ParentIndex; sibs: readonly Map<string, string>[]; shareCalls: Map<string, ShareCall>
  result: ParentageResult; soloRole: 'paternal' | 'maternal'; otherRole: 'paternal' | 'maternal'; log: LogFn
}): void {
  const { acc, pat, sibs, shareCalls, result, soloRole, otherRole, log } = x
  const her = otherRole === 'maternal' ? 'mother' : 'father'
  if (!sibs.length) {
    result.limits.push(`The run's other arrays were not used to read the ${her}: none were loaded with this one.`)
    return
  }
  const freqB = (acc.panel as FounderPanel).freqB!
  const markers: RunMarker[] = []
  for (const [id, where] of acc.markerPos) {
    const g = pat.gt.get(id)
    const f = freqB.get(id)
    if ((g !== 'AA' && g !== 'BB') || f === undefined) continue
    const q = g === 'AA' ? f : 1 - f
    if (q > 0) markers.push({ id, chrom: where.chrom, pos: where.pos, parent: g, q })
  }
  const run = motherFromRun({ markers, sample: acc.myGt, baf: (id) => acc.myBaf.get(id), others: sibs,
    parentKey: pat, rebuilt: !!pat.inferred, mMin: pat.inferred?.mMin })
  result.motherFromRun = { status: run.status, draws: run.draws, drawArrays: run.drawArrays }
  acc.otherOpposite = run.opposite ?? null
  if (run.status !== 'used') {
    result.limits.push(`The run's other arrays were not used to read the ${her}: ${
      run.status === 'sample-unfit' ? 'this array\'s noise is too high to place it in an embryo'
      : run.status === 'sample-not-a-draw' ? 'this array\'s own calls are too noisy, or show too little of the other '
        + 'parent, to tell which arrays share its mother'
      : `${run.draws} other draw${run.draws === 1 ? '' : 's'} of the ${her} (an embryo or a haploid product of hers), `
        + 'and 2 are needed'}.`)
    return
  }
  result.limits.push(`The ${her}'s genotype was read from ${run.draws} other draws of the ${her} in this run `
    + '(embryos and haploid products sharing her copies), each weighted by its own drop-in and dropout.')
  const named = new Set(result.chroms.filter((c) => c.alleleGain?.verdict === 'gain').map((c) => c.chrom))
  const chroms: RouteChrom[] = result.chroms.filter((c) => isAutosome(c.chrom)).map((c) => {
    const sh = run.shares.get(c.chrom)
    return { chrom: c.chrom, share: sh?.share ?? NaN, weight: sh?.weight ?? 0, state: shareCalls.get(c.chrom)?.state,
      measured: c.verdict !== 'not_measured', intensityGain: c.aneuploidy === 'gain', lrrShift: c.lrrShift ?? NaN }
  })
  for (const g of routeGains(chroms, named)) {
    const c = result.chroms.find((y) => y.chrom === g.chrom)
    if (!c) continue
    const parent = g.loadedParent ? soloRole : otherRole
    const paternalShift = soloRole === 'paternal' ? g.shift : -g.shift
    c.alleleGain = {
      chrom: c.chrom, verdict: 'gain', parent, shift: paternalShift, z: soloRole === 'paternal' ? g.z : -g.z,
      scale: g.scale, centre: soloRole === 'paternal' ? g.centre : 1 - g.centre, informative: Math.round(g.weight),
      lrrShift: c.lrrShift ?? NaN, lrrZ: g.lrrZ,
      why: `against the ${her}'s genotype read from ${run.draws} other draws of her in this run, the ${parent} `
        + `alleles are over-represented by ${Math.abs(g.shift).toFixed(3)} against the ${g.references} other `
        + `chromosomes reading both parents (${Math.abs(g.z).toFixed(1)} of their spread), and the intensity `
        + `reads ${g.lrrZ.toFixed(1)} of the cell's own spread above them, so the ${parent} parent supplied an `
        + 'extra copy',
    }
    log('WARN', `chr${c.chrom}: extra copy, ${parent}. ${c.alleleGain.why}.`)
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
  /**
   * What the operator says this material is, overriding the inference.
   *
   * The inference still runs and the disagreement is reported: see declaredStage.ts.
   */
  declaredStage?: Stage
}): Promise<ParentageResult> {
  const { acc, profile, pat, mat, soloRole, sibs, sampleName, log } = input
  const {
    t, myGt, myBaf, markerPos, cnByChrom, selfMarkers, hetByChrom,
    obligateByChrom,
    matAbsenceByChrom, matAbsenceTally,
    absenceByChrom, dosageByChrom,
  } = acc
  // THE PANEL THE SAMPLE IS ON, from a parental array in a state to say so. An array that is
  // itself heavily dropped out cannot define a panel, and the case that proves it is the operator
  // dropping the SAMPLE into the parent slot: the sample then anchors its own ceiling, and an
  // implausibly heterozygous array certifies itself as ordinary. See PANEL_ANCHOR_CALL_RATE.
  const parentHets = [pat, mat]
    .filter((x): x is ParentIndex => !!x && x.callRate >= PANEL_ANCHOR_CALL_RATE)
    .map((x) => x.heterozygosity)
    .filter((x) => Number.isFinite(x) && x > 0)
  // THE PARENT-FREE ANSWER KEYS GO TO THE STAGE, and only those: the stage must come out the same
  // whatever parents are loaded, and a reading against the panel never sees one.
  const panelDropIn = acc.panelSameHom && acc.panelSameHom.n >= MIN_SAME_HOM
    ? onParentScale(dropIn(acc.panelSameHom), (acc.panel as { calibration?: PanelCalibration } | null)?.calibration)
    : NaN
  const panelInformation = acc.informationPanel ? informationRatio(acc.informationPanel) : NaN
  // A RECONSTRUCTED PARENT DOES NOT SCALE THE STAGE. Its heterozygosity is an estimate from product
  // disagreement over the markers its depth retained, and it moves with which products went in: one
  // father reads 15.5% rebuilt from five of his pronuclei and 19.3% from all eight, which raised the
  // ceiling by 1.15 and moved a cell's rung. A measured array's heterozygosity is the platform's.
  const stageAnchors = [pat, mat]
    .filter((x): x is ParentIndex => !!x && !x.inferred && x.callRate >= PANEL_ANCHOR_CALL_RATE)
    .map((x) => x.heterozygosity)
    .filter((x) => Number.isFinite(x) && x > 0)
  // THE ALLELES SAY HOW MANY GENOMES, where a parent is loaded: see stage.allelePloidy. Read from the
  // same tallies as the copy state below, against the father in either configuration, so one
  // reconstructed parent and two measured ones give one stage.
  // WHOSE GENOME IS THIS, before any allele is read against a parent. See parentGenome.ts.
  const genomeOf = (m: Map<string, CarryTally> | null) => (m
    ? parentGenome([...m].filter(([c]) => isAutosome(c)).map(([, t]) => t)) : UNREAD)
  const patGenome = genomeOf(acc.carryByChrom)
  const matGenome = mat ? genomeOf(acc.carryMatByChrom) : undefined
  // The copy states count genomes only where the parent's is among them: against an adult whose
  // genome is not in the array they read that adult's copy lost whatever the array holds.
  const parentShown = patGenome.verdict === 'carries' || patGenome.verdict === 'unread'
  // Whether the alleles can be read against a parent at all: one is loaded, and its genome is here.
  const allelesReadable = parentShown && !pat.absent
  const copyStates = acc.shareByChrom && parentShown
    ? [...acc.shareByChrom].filter(([c]) => isAutosome(c)).map(([, t]) => shareCall(t).state)
    : undefined
  const stageAgreement = reconcileStage(
    inferStage({ ...profile, dropIn: panelDropIn, information: panelInformation, copyStates }, {
      panelHeterozygosity: stageAnchors.length ? Math.min(...stageAnchors) : undefined,
      ...input.stageOpts,
    }), input.declaredStage, stageDefaults,
  )
  // The material goes in because the whole-chromosome magnitude gate depends on it. Bulk here means
  // the stage reads bulk AND the array calls at a bulk rate: a single cell amplified by MDA can
  // read bulk heterozygosity through drop-in, and its call rate is what it cannot hide. See
  // BULK_SPREAD_MAX for the third condition, which classify checks itself.
  const bulkMaterial = stageAgreement.used.stage === 'bulk'
    && profile.callRate >= PANEL_ANCHOR_CALL_RATE
  // THE REFERENCE'S OWN FALSE ABSENCE GOES IN HERE. A reconstructed parent is homozygous at
  // markers where she is really heterozygous, and a haploid product of hers carries the other
  // allele at half of those, which reads as an opposite homozygote through no fault of the
  // sample. `classify` has taken that floor since parentage.check.ts section 16 and this caller
  // never supplied it, so every row scored against a reconstruction was judged by the ceiling
  // calibrated for a MEASURED array. Zero for a measured parent, which contributes none.
  const result = classify(t, pat.heterozygosity, {
    role: soloRole,
    bulk: bulkMaterial,
    dropIn: panelDropIn,
    ...(pat.inferred && Number.isFinite(pat.inferred.spuriousAbsence)
      ? { spuriousAbsence: pat.inferred.spuriousAbsence } : {}),
    ...(pat.inferred ? { inferredParent: true } : {}),
    yBearing: sexCall(acc.sex).yBearing,
  })
  result.refusals = []
  result.parentGenome = patGenome
  if (pat.absent) {
    result.noParent = true
    result.limits = result.limits.filter((l) => !l.startsWith('No markers where the parent is homozygous'))
    result.limits.push('No parental array is loaded. The stage, the sex chromosomes and copy number from '
      + 'the call rate and the intensity are read; which parent a chromosome came from is not.')
  }
  if (matGenome) result.parentGenomeOther = matGenome
  for (const [g, role] of [[patGenome, soloRole], [matGenome, soloRole === 'paternal' ? 'maternal' : 'paternal']] as const) {
    if (!g || g.verdict === 'carries' || g.verdict === 'unread') continue
    const who = role === 'paternal' ? 'father' : 'mother'
    const said = parentGenomeSentence(g, who)
    log('WARN', said)
    result.limits.push(said)
    result.parentGenomeNotices = [...result.parentGenomeNotices ?? [], parentGenomeNotice(g, who)!]
  }
  // ZYGOSITY FROM THE ALLELES, where they decide: one genome or two, read as the stage reads it
  // (allelePloidy). classify read it from the heterozygous band, which drop-in fills and allele
  // dropout empties; the other parent's alleles in the calls are touched by neither.
  const allelesSay = allelePloidy(copyStates)
  if (allelesSay) {
    result.zygosity = allelesSay === 'two genomes' ? 'diploid' : 'uniparental_homozygous'
    result.limits = result.limits.filter((l) => !l.includes('Zygosity is '))
    result.notes.push(`Zygosity is read from the alleles against the loaded parent: ${allelesSay === 'two genomes'
      ? 'most chromosomes carry a second genome' : 'most chromosomes carry one parent\'s copy alone'}.`)
  }
  // A reconstruction whose banner could not be read is NOT a measured array, and proceeding as
  // though it were is the failure this whole path exists to stop. Say so where the operator reads.
  if (pat.inferred && !Number.isFinite(pat.inferred.hRetained)) {
    result.limits.push(
      'The parental file is a reconstructed reference whose banner could not be read, so the '
      + 'parent heterozygosity it established and the absence floor it contributes are both '
      + 'unavailable. Every call here is scored as though the reference were measured, which it '
      + 'is not. Rebuild the reference, or load the parent\'s own array.',
    )
  }
  // Stage from the array itself, since the dropout each stage carries is what every
  // downstream likelihood is parameterised by. Bundled into the result so every output
  // carries it, with the basis and the confounds attached to the number.
  // THE OPERATOR'S LABEL WINS, AND THE DISAGREEMENT IS REPORTED. The inference is always run, even
  // when it is overridden, because a declaration that contradicts the array is the single most
  // useful thing this tool can tell somebody: either the tube is mislabelled or the reaction
  // failed, and both are worth knowing before any parental call is read.
  result.stage = stageAgreement.used
  result.stageAgreement = stageAgreement
  log('DONE', `stage: ${result.stage.stage}. ${result.stage.why}`)
  if (!stageAgreement.agrees) {
    log('WARN', stageAgreement.notice)
  }
  // One contribution or two, per chromosome. The boundaries are measured for a BULK reference
  // parent, and against a single-cell reference they do not separate at all (audit section E2),
  // so with one parent, or an amplified one, this is reported and never used to admit or reject.
  //
  // WITH TWO BULK PARENTS IT IS ALSO A DETECTOR, and the one that reaches a whole chromosome on
  // amplified material. The call-rate and intensity gates in classify() miss about half of the
  // whole-chromosome losses a published method called by hand on single cells: 53 of 111 sit
  // inside 0.2 log2R of the genome, where no intensity floor can reach them. The alleles do not
  // depend on intensity. A chromosome reading one parent's allele alone at the markers where the
  // parents are homozygous for opposite alleles, on a genome carrying both parents elsewhere, is
  // missing the other parent's copy. Dropout cannot produce that pattern: it removes a
  // heterozygote's alleles at random rather than one parent's allele at every marker. Measured
  // against the same hand calls through this pipeline: 78 of 85 monosomies on the arrays it
  // scores name the same parent, 1 the other, and the Mendelian channel agrees on all 78
  // (obligateHet.SURVIVOR_MIN, audit/zuccaro-table-s3.ts).
  //
  // WHAT THE ALLELES DO NOT SAY IS THE COPY NUMBER. One parent's allele alone is what a monosomy
  // reads, and it is also what an isodisomy reads: the other parent's copy is absent in both, and
  // only the intensity separates one remaining copy from two. So this sets `lostParent`, which is
  // the fact the alleles establish, and leaves `aneuploidy` to the call rate and the intensity.
  // The page reports a chromosome with a lost parent and no intensity verdict as a parental
  // absence with the one-or-two question stated, which is also as far as a hand karyotype of the
  // same cell can go. The Mendelian channel then runs on the chromosome like any detected event,
  // so the row carries that channel's posterior and band, and a disagreement between the two
  // readings is reported as unclear rather than resolved in favour of either.
  const twoBulk = !!mat && pat.heterozygosity >= PANEL_HET_FLOOR
    && mat.heterozygosity >= PANEL_HET_FLOOR
  const bulkReference = pat.heterozygosity >= PANEL_HET_FLOOR
    && (!mat || mat.heterozygosity >= PANEL_HET_FLOOR)
  const otherRole = soloRole === 'paternal' ? 'maternal' : 'paternal'
  const two = mat !== null && mat !== undefined
  // A PARENTAL ARRAY OVER THE DIPLOID CEILING is reported, never withheld: its homozygous calls
  // still make the opposite-homozygote markers informative, which is all the allele count reads,
  // and what it cannot anchor is anything read from its heterozygosity. One rule for the
  // application and the harnesses alike: see stage.parentOverCeiling.
  if (parentHets.length === 2) {
    const panelHet = Math.min(...parentHets)
    for (const [p, role] of [[pat, soloRole], [mat as ParentIndex, otherRole]] as const) {
      if (!parentOverCeiling(p.heterozygosity, panelHet)) continue
      const text = `The ${role} array is ${pct(p.heterozygosity, 1)} heterozygous, over the `
        + `${pct(MAX_DIPLOID_HET * panelScale(panelHet), 0)} a single diploid reaches on this `
        + 'panel: two people or a contaminated reaction. Its homozygous calls still inform the '
        + 'allele count; anything resting on its heterozygosity is provisional.'
      result.limits.push(text)
      log('WARN', text)
    }
  }
  // This cell's drop-in, from markers where both measured parents are the same homozygote. Only
  // with two measured parents: with one, the informative set is defined by one array and a marker
  // where the other parent is heterozygous would read as drop-in that is not.
  const cellDropIn = two && acc.sameHom && acc.sameHom.n >= MIN_SAME_HOM ? dropIn(acc.sameHom) : 0
  for (const c of result.chroms) {
    const tally = two
      ? obligateByChrom.get(c.chrom) ?? emptyHet()
      : { ...emptyHet(), ...(hetByChrom.get(c.chrom) ?? {}) }
    const call = hetCall(tally, two ? 2 : 1, cellDropIn)
    // Kept wherever the markers were enough, between the boundaries included, so the counts a
    // reader checks by hand are the counts the verdict rests on.
    // A CHROMOSOME THAT WAS BARELY GENOTYPED IS NOT A CHROMOSOME WITH AN EXTREME VALUE, and the
    // informative count cannot tell them apart. One chromosome of a laboratory series carries
    // 1,605 informative markers, clears the absolute floor, and reports 55.1 percent heterozygous
    // where its siblings read 0.07: it was called at 0.175 of the markers the array holds for it.
    // The call rate separates that from everything else with nothing in between. See
    // `informativeEnough`.
    const stats = profile.byChrom.get(c.chrom)
    const chromCallRate = stats && stats.markers > 0 ? stats.called / stats.markers : NaN
    // A SEX CHROMOSOME IS NOT A CHROMOSOME WITH NO MARKERS ON IT. The opposite-homozygote channel
    // counts autosomes only, so the offered count for chrX is zero by construction, and the
    // offered-count wording would state as a fact about the parents something that is a fact about
    // which channel ran. chrX is scored by callSexChromosomes, on intensity and heterozygosity,
    // and saying so is the honest refusal.
    const withheld = isAutosome(c.chrom)
      ? whyNotInformative(call.informative, chromCallRate,
        { offered: tally.offered, genomeCallRate: profile.callRate })
      : 'the opposite-homozygote channel counts autosomes only, so it makes no call on a sex '
        + 'chromosome whatever the array quality. chrX and chrY are scored separately, on '
        + 'intensity and heterozygosity rather than on which parent\'s allele is present'
    if (withheld) c.contributionWithheld = withheld
    else c.contribution = call
  }
  // THE CHROMOSOMES THAT NEVER REACHED THE MEASUREMENT LIST. `chroms` excludes any chromosome under
  // 200 informative markers, correctly: its rate is not a measurement and would widen every null
  // taken from the spread of the others. But it was leaving the report as well as the measurement,
  // so an autosome the array failed to genotype produced no row, no call and no reason, which is
  // the one outcome a reader cannot distinguish from a chromosome that was fine. They are recorded
  // here instead, with the same two-way reason the withheld chromosomes get.
  const inChroms = new Set(result.chroms.map((c) => c.chrom))
  for (const [chrom, tally] of [...hetByChrom].sort((a, b) => Number(a[0]) - Number(b[0]))) {
    if (inChroms.has(chrom) || !isAutosome(chrom)) continue
    const stats = profile.byChrom.get(chrom)
    const cr = stats && stats.markers > 0 ? stats.called / stats.markers : NaN
    result.unmeasured.push({
      chrom,
      why: whyNotInformative(tally.informative, cr,
        { offered: tally.offered, genomeCallRate: profile.callRate })
        // Reached only if the chromosome clears both conditions here yet was still excluded, which
        // means the two tallies disagree. Said plainly rather than papered over.
        ?? `${tally.informative.toLocaleString()} informative markers on this chromosome, which `
          + 'clears the floor for a contribution call but not the one the per-chromosome '
          + 'measurements need',
    })
  }
  // "ON A GENOME THAT CARRIES BOTH PARENTS ELSEWHERE" IS A CONDITION, AND IT IS MEASURED HERE
  // RATHER THAN TAKEN FROM THE GENOME-LEVEL CLASS. A bulk array of the father loaded as a sample is
  // diploid and androgenetic: every chromosome reads his allele alone, the zygosity gate passes it,
  // and naming the maternal copy absent chromosome by chromosome put eight events on a genome that
  // has none. The genome-level class is no gate either, in the other direction: a cell that has
  // lost one parent's copy on half its chromosomes reads gynogenetic or androgenetic as a whole,
  // and excluding it by class threw away thirteen correct calls on the Zuccaro cells. What
  // separates the two is whether ANY autosome still reads both parents at these markers. An
  // androgenetic genome cannot reach 45 percent heterozygous at opposite-homozygote markers on
  // even one chromosome; a chaotic cell reads it on every chromosome it kept.
  const bothParentsSomewhere = result.chroms
    .filter((c) => isAutosome(c.chrom) && c.contribution?.ploidy === 'biparental').length
    >= BIPARENTAL_AUTOSOMES_MIN
  // AND WITH BOTH PARENTS LOADED THAT IS TWO CONTRIBUTIONS, read from the alleles with drop-in
  // removed. Where the band could not say, because drop-in filled it, this does: both parents'
  // alleles on that many chromosomes is a genome carrying both. Not with one parent, where no drop-in
  // is removed from this reading and drop-in alone would pass for the second parent.
  if (two && cellDropIn > 0 && result.zygosity === 'unknown' && bothParentsSomewhere) {
    result.zygosity = 'diploid'
    result.limits = result.limits.filter((l) => !l.includes('Zygosity is '))
    result.notes.push('Zygosity is read from the alleles: both parents\' alleles are present, with '
      + `drop-in removed, on ${BIPARENTAL_AUTOSOMES_MIN} or more autosomes.`)
  }
  for (const c of result.chroms) {
    const call = c.contribution
    if (!twoBulk || !bothParentsSomewhere || !knownBiparental(result.zygosity)
      || !isAutosome(c.chrom) || c.verdict === 'not_measured' || !call) continue
    if (call.ploidy === 'uniparental' && call.survivor) {
      // The first parental slot holds the array whose role is soloRole; the survivor's allele is
      // what remains, so the copy that is missing is the other parent's.
      const lost = call.survivor === 'a' ? otherRole : soloRole
      c.lostParent = lost
      if (c.aneuploidy === 'loss') c.aneuploidyParent = lost === soloRole ? 'this' : 'other'
    } else if (call.ploidy === 'uncalled' && call.corrected.fraction >= HAPLOID_MAX
      && call.corrected.fraction < DIPLOID_MIN
      && Math.max(call.corrected.onlyA, call.corrected.onlyB) >= MIXED_SHARE) {
      // BETWEEN THE BOUNDARIES WITH ONE PARENT DOMINANT is said out loud rather than passed over.
      // Read on the drop-in corrected shares, as the boundaries are. See MIXED_SHARE.
      c.mixed = true
    }
  }
  // WITH ONE PARENT, WHICH COPY IS GONE, FROM THAT PARENT AND THE POPULATION. The two-parent reading
  // above needs markers where the parents are opposite homozygotes, so with one parent it never runs,
  // and one parent, reconstructed from its products, is the laboratory's whole configuration. See
  // oneParentState.ts: against a father rebuilt from five pronuclei it matched the measured father on
  // 634 of 638 chromosomes, and against Table S3 it is 120 of 129 abnormal states with the rebuild,
  // 121 with the measured father alone and 123 with both parents. The same guard as the two-parent
  // reading: a genome must carry both parents somewhere before a chromosome can have lost one.
  const shareCalls = new Map<string, ShareCall>()
  if (acc.shareByChrom) {
    for (const c of result.chroms) {
      if (!isAutosome(c.chrom)) continue
      const t = acc.shareByChrom.get(c.chrom)
      if (!t) continue
      const sc = shareCall(t)
      shareCalls.set(c.chrom, sc)
      // Kept for the normal-majority test below either way; shown as a state only against a parent.
      if (parentShown) c.copyState = sc
    }
  }
  // Against no parent, or an adult whose genome is not in the array, the states show no majority.
  const normal = allelesReadable ? normalFraction([...shareCalls.values()]) : NaN
  /** Whether this cell's own intensity level is a normal chromosome, by the test the intensity
   *  reference is already gated on. Rule 4 of oneParentGain.ts is scoped to it. */
  const trustedLevel = Number.isFinite(normal) && normal >= NORMAL_MAJORITY_MIN
  // BOTH PARENTS SOMEWHERE, FROM ONE PARENT'S ALLELES. Where neither the band nor the genome-wide
  // reading decides, as on a chaotic genome that lost chromosomes of both parents, chromosomes reading
  // as a second copy reads say a second contribution exists. The copy state cancels drop-in, so this
  // holds with one parent as with two; the same count guards the one-parent loss reading below.
  const secondCopies = [...shareCalls.values()].filter((x) => x.state === 'disomic' && x.r >= DISOMIC_LOW).length
  if (parentShown && result.zygosity === 'unknown' && allelesSay !== 'one genome' && secondCopies >= BIPARENTAL_AUTOSOMES_MIN) {
    result.zygosity = 'diploid'
    result.limits = result.limits.filter((l) => !l.includes('Zygosity is '))
    result.notes.push(`Zygosity is read from the alleles: ${secondCopies} autosomes carry a second copy `
      + 'against the loaded parent, drop-in removed.')
  }
  if (parentShown && !mat && [...shareCalls.values()].filter((x) => x.state === 'disomic').length >= BIPARENTAL_AUTOSOMES_MIN) {
    for (const c of result.chroms) {
      const sc = shareCalls.get(c.chrom)
      if (!sc || c.verdict === 'not_measured') continue
      const lost = sc.state === 'other-lost' ? otherRole : sc.state === 'parent-lost' ? soloRole : undefined
      // A lost copy cannot raise intensity: on a chromosome whose intensity reads more DNA than this
      // cell's own level, the loss state is the extra copy's share, and the gain test reads it. The
      // intensity does not have to clear the bar a gain CALL needs; see oneParentGain.ts rule 1.
      if (!lost) continue
      // AND ONLY WHERE THIS CELL'S OWN LEVEL MEANS SOMETHING. Rule 4 reads the intensity against
      // the cell's own level, so it is worth no more than that level is. On a cell that has lost a
      // third of its chromosomes the level is dragged down and a one-copy chromosome reads ABOVE
      // it, up to +0.48 measured on this corpus, which is where a third copy reads. The test for
      // that is the one the tool already applies before reading copy number at all.
      if (gainReadAsLoss(sc.state, c.aneuploidy === 'gain', trustedLevel ? c.lrrShift : NaN)) {
        if (c.aneuploidy !== 'gain') {
          result.limits.push(`chr${c.chrom}: the alleles read only one parent's copy here, but this `
            + `chromosome's intensity sits ${c.lrrShift.toFixed(2)} log2 ABOVE the chromosomes this `
            + 'array is measured against, and a copy that is missing cannot do that. No lost parent '
            + 'is reported on it. The intensity does not clear the bar an extra copy is called at, '
            + 'so no gain is reported either.')
        }
        continue
      }
      c.lostParent = lost
      if (c.aneuploidy === 'loss') c.aneuploidyParent = lost === soloRole ? 'this' : 'other'
    }
  }
  // AN EXTRA CHROMOSOME, FROM THE SAME ALLELES. See gainOrigin.scanAlleleGains. It runs under the
  // conditions the allele count runs under, and is judged against the cell's own biparental
  // chromosomes, so a chaotic cell supplies a wide scale and refuses rather than over-calls.
  let chromShares: ChromShare[] = []
  if (twoBulk && bothParentsSomewhere && knownBiparental(result.zygosity)) {
    const shares = result.chroms.filter((c) => isAutosome(c.chrom)).map((c) => {
      const d = dosageByChrom.get(c.chrom) ?? []
      return {
        chrom: c.chrom,
        share: d.length ? median(d.map((m) => m.patShare)) : NaN,
        informative: d.length,
        biparental: c.contribution?.ploidy === 'biparental' && c.verdict !== 'not_measured',
        lrrShift: c.lrrShift,
      }
    })
    chromShares = shares
    for (const g of scanAlleleGains(shares)) {
      const c = result.chroms.find((x) => x.chrom === g.chrom)
      if (c) c.alleleGain = g
    }
  } else if (!two && shareCalls.size && knownBiparental(result.zygosity)) {
    // AND WITH ONE PARENT, the laboratory's configuration. The copy state reads, at markers where the
    // loaded parent is homozygous, the share r of the allele it lacks, without selecting on this
    // cell's calls. An extra copy of the loaded parent's chromosome leaves the other parent's alleles
    // a third of the copies, r = 2/3 of the cell's disomic level; an extra copy of the other's, 4/3.
    // So the loaded parent's share, 1 - r / (2 x that level), is 1/2, 2/3 or 1/3: the two-parent
    // share's units, read by the same test against the same cell's disomic chromosomes. Whole
    // chromosomes only. See audit/one-parent-gain-preregistration.txt.
    //
    // A REBUILT parent's homozygotes are weighted by the chance each is real (oneParentGain rule 2),
    // in this share only; the copy state keeps the unweighted r. A chromosome intensity calls gained
    // is tested even where its copy state reads a lost parent (rule 1), and is never a reference.
    const rOf = (chrom: string): number => {
      const sc = shareCalls.get(chrom)
      if (!sc) return NaN
      return acc.shareWByChrom
        ? weightedShareR(acc.shareByChrom!.get(chrom)!, acc.shareWByChrom.get(chrom)) : sc.r
    }
    // Rules 1 and 2 decide only on chromosomes the intensity scan calls gained; everywhere else the
    // verdict is the unweighted test's (rule 4).
    const sharesOf = (rFor: (chrom: string) => number, rules: boolean) => {
      const level = median(result.chroms.filter((c) => shareCalls.get(c.chrom)?.state === 'disomic')
        .map((c) => rFor(c.chrom)).filter((r) => Number.isFinite(r)))
      if (!Number.isFinite(level) || level <= 0) return null
      return result.chroms.filter((c) => isAutosome(c.chrom)).map((c) => {
        const sc = shareCalls.get(c.chrom)
        const r = rFor(c.chrom)
        const own = Number.isFinite(r) ? 1 - r / (2 * level) : NaN
        const measured = c.verdict !== 'not_measured'
        return {
          chrom: c.chrom,
          share: soloRole === 'paternal' ? own : 1 - own,
          informative: sc?.n2 ?? 0,
          biparental: measured && (sc?.state === 'disomic'
            || (rules && gainReadAsLoss(sc?.state, c.aneuploidy === 'gain',
              trustedLevel ? c.lrrShift : NaN))),
          reference: measured && sc?.state === 'disomic',
          lrrShift: c.lrrShift,
        }
      })
    }
    const shares = sharesOf((chrom) => shareCalls.get(chrom)?.r ?? NaN, false)
    if (shares) {
      for (const g of scanAlleleGains(shares)) {
        const c = result.chroms.find((x) => x.chrom === g.chrom)
        if (c) c.alleleGain = g
      }
      const ruled = sharesOf(rOf, true)
      if (ruled) {
        for (const g of scanAlleleGains(ruled, ruled.filter((x) => x.reference))) {
          const c = result.chroms.find((x) => x.chrom === g.chrom)
          if (c && c.aneuploidy === 'gain') c.alleleGain = g
        }
      }
      // ORIENTED WHERE INTENSITY HAS ALREADY CALLED THE GAIN, AND ONLY THERE. The test above can detect
      // a gain from the alleles alone at ALLELE_GAIN_Z, a bar set on two parents' share. One parent's
      // share is noisier: on half a cell's markers it put a gain, with a parent, on a chromosome Table
      // S3 publishes normal (GSM4552372 chr16) where intensity called none. So with one parent
      // intensity detects and the alleles only say whose, at the bar that keeps direction honest.
      for (const c of result.chroms) {
        const g = c.alleleGain
        if (g?.verdict !== 'gain' || c.aneuploidy === 'gain') continue
        c.alleleGain = { ...g, verdict: 'none', parent: undefined,
          why: 'the alleles lean to one parent, but intensity calls no extra copy here, and with one '
            + 'parent the alleles name the parent of a gain intensity has called, not detect one' }
      }
      for (const c of result.chroms) {
        const g = c.alleleGain
        if (c.aneuploidy !== 'gain' || !g || g.verdict !== 'none' || !Number.isFinite(g.z)
          || Math.abs(g.z) < ONE_PARENT_ORIENT_Z) continue
        const parent = g.z > 0 ? 'paternal' : 'maternal'
        c.alleleGain = { ...g, verdict: 'gain', parent,
          why: `intensity reads an extra copy, and against the loaded ${soloRole === 'paternal' ? 'father' : 'mother'} `
            + `the alleles read the ${parent} share ${Math.abs(g.z).toFixed(1)} of this cell's spread `
            + 'from its disomic chromosomes, toward the extra copy being that parent\'s' }
      }
      // THE OTHER PARENT'S COPY FROM THE HETEROZYGOUS CALLS, UPWARD ONLY (oneParentGain rule 3), where
      // intensity calls the gain and the share above named no parent.
      if (acc.hetShareByChrom) {
        const refs = result.chroms.filter((c) => shareCalls.get(c.chrom)?.state === 'disomic'
          && c.verdict !== 'not_measured')
        for (const c of result.chroms) {
          if (!isAutosome(c.chrom) || c.aneuploidy !== 'gain' || c.alleleGain?.verdict === 'gain') continue
          const call = otherCopyFromHetCalls(acc.hetShareByChrom.get(c.chrom) ?? [],
            refs.filter((x) => x.chrom !== c.chrom).map((x) => acc.hetShareByChrom!.get(x.chrom) ?? []))
          if (!call?.named) continue
          const g = c.alleleGain ?? { chrom: c.chrom, shift: NaN, z: NaN, scale: NaN, centre: NaN,
            informative: shareCalls.get(c.chrom)?.n2 ?? 0, lrrShift: c.lrrShift }
          c.alleleGain = { ...g, verdict: 'gain', parent: otherRole,
            why: `intensity reads an extra copy, and where this cell is heterozygous at markers the loaded `
              + `${soloRole === 'paternal' ? 'father' : 'mother'} is homozygous for, the allele that parent lacks `
              + `holds ${call.median.toFixed(3)} of the signal against ${call.centre.toFixed(3)} on its other `
              + `two-parent chromosomes, ${call.z.toFixed(1)} of their spread above. Only an extra copy of the `
              + `${otherRole} parent lifts it` }
        }
      }
    }
  }
  {
    for (const c of result.chroms) {
      const g = c.alleleGain
      if (!g || (g.verdict !== 'gain' && g.verdict !== 'imbalance')) continue
      log('WARN', `chr${c.chrom}: ${g.verdict === 'gain' ? `extra copy, ${g.parent}` : 'allele imbalance, no parent named'}. ${g.why}.`)
    }
  }
  for (const c of result.chroms) {
    if (!c.mixed || !c.contribution) continue
    const k = c.contribution
    const [father, mother] = soloRole === 'paternal' ? [k.onlyA, k.onlyB] : [k.onlyB, k.onlyA]
    log('WARN', `chr${c.chrom}: mixed signal, not called. Of ${k.informative} markers where the `
      + `parents are opposite homozygotes, ${pct(k.fraction, 1)} read both alleles, `
      + `${pct(father, 1)} the father's alone and ${pct(mother, 1)} the mother's alone: between `
      + 'the boundaries with one parent dominant. A mosaic, a partial loss or noise, not '
      + 'separable here.')
  }
  for (const c of result.chroms) {
    // The opposite-homozygote count needs both parents' arrays.
    if (!c.lostParent || !c.contribution || !two) continue
    const k = c.contribution
    const [father, mother] = soloRole === 'paternal' ? [k.onlyA, k.onlyB] : [k.onlyB, k.onlyA]
    log('DONE', `chr${c.chrom}: the ${c.lostParent} copy is absent. Of ${k.informative} markers `
      + `where the parents are opposite homozygotes, ${pct(k.fraction, 1)} read both alleles, `
      + `${pct(father, 1)} the father's alone and ${pct(mother, 1)} the mother's alone`
      + (c.aneuploidy ? '.' : '. Whether one copy remains or two is for the intensity, which did '
        + 'not reach a verdict on this chromosome.'))
  }
  const uni = result.chroms.filter((c) => c.contribution?.ploidy === 'uniparental' && !c.lostParent)
  if (uni.length) {
    log('WARN', `one parental contribution only on ${uni.map((c) => `chr${c.chrom}`)
      .join(', ')}: ${uni[0].contribution!.why}`
      + (uni[0].contribution!.provisional
        ? '. One parent loaded, so this boundary is provisional' : ''))
  }
  // AN ARRAY THAT IS NOT MEASURING A GENOME PRODUCES NO EVENTS OF ANY KIND.
  //
  // The stage inference decides this and says so in those words: one example array reads 31.2%
  // heterozygous where a single diploid tops out near 25% on this platform, and its verdict is
  // "no genome reads this way, so this is not a stage".
  //
  // THE SAME CONDITION `integrity.ts` ALREADY CALLS UNINTERPRETABLE, and its limits are the ones
  // enforced here rather than a second opinion beside them: "Nothing on this sample is reportable.
  // Do not read its changes, its parental calls or its counts". This gate used to empty the
  // taxonomy findings alone, so whole-chromosome aneuploidy and copy-number segments walked past
  // it and shipped as calls: 30 aneuploidy calls and 86 segment calls across 186 arrays that had
  // failed their own quality inference. They are found by the same measurement the tool has just
  // declared broken.
  //
  // POSITION DOES NOT SURVIVE HERE, though it survives the other two uninterpretable verdicts. On
  // a low call rate or a genome with no undisturbed remainder the array IS reading a genome, so
  // the coordinates of what it found stand and only the parent is withheld. A failed array is the
  // case where there is no genome being read, so there is no event to have a position.
  //
  // Clearing the aneuploidy here rather than filtering it downstream is what makes the rest of
  // this function fall silent on its own: `whole` empties, the segment scans are skipped, the
  // dosage and Mendelian channels have nothing to walk, and `assessIntegrity` counts zero events
  // under the same headline.
  // A CHAOTIC GENOME IS A GENOME. The ceiling at MAX_DIPLOID_HET refuses an array reading above
  // what a single diploid can reach, and its own caveat says the causes behind that are not
  // separable. One of them is a genome with many chromosomes at odd copy number, which raises the
  // apparent heterozygous rate on those chromosomes and not on the others. The other causes are
  // elevated everywhere. So a LOCALISED excess is a genome whose copy numbers are wrong, and a
  // uniform one is a mixture, a contaminated reaction or a destroyed array.
  //
  // WHAT THAT OPENS, AND ONLY THAT. This said the allele channels "need informative markers and a
  // call rate and never the intensity baseline", and that is false: the absence scan runs only on
  // chromosomes in `measured`, whose verdict depends on `clustered`, which reads the per-chromosome
  // extreme-BAF share, and BAF is intensity-derived. An external review traced it and the chain is
  // absence scan -> measured -> verdict -> clustered -> BAF.
  //
  // THE TRUE DISTINCTION IS NARROWER AND IT IS STILL ENOUGH. The absence channel needs BAF
  // readable PER CHROMOSOME; `clustered` is a within-chromosome share with no cross-chromosome
  // normalisation, so a genome with no normal majority costs it nothing. The copy-number channel
  // needs BAF and intensity COMPARABLE ACROSS chromosomes, because it normalises against this
  // sample's own per-chromosome median, and a genome without a normal majority has no such median.
  // That is what licenses gating them separately. "Never consults the intensity" does not.
  //
  // So the absence scan runs here and the copy-number scan does not. The whole-chromosome verdicts
  // stay withheld too: two further gates stand in their way, and neither has been measured on this
  // material.
  //
  // Measured over the public corpus: the four hand-called segmental deletions on cells refused
  // this way are recovered, each naming the parent the hand karyotype names.
  const spread = hetSpread(profile.byChrom)
  const chaoticGenome = result.stage?.stage === 'failed'
    && Number.isFinite(spread) && spread >= HET_SPREAD_LOCALISED
  const uninterpretable = result.stage?.stage === 'failed'
  // THE NUMBERS EVERY STAGE REFUSAL FIRED ON, so each record carries its measurement.
  const stageGate = `stage.${result.stage?.refusedBy ?? 'unknown'}`
  // A STAGE THAT NAMED A PHYSICAL FAILURE is that failure for every claim it clears. The ceiling, a
  // declared failure and the call-rate floor name none, and stay proxies.
  const stageCause: PhysicalCause | 'proxy' = result.stage?.refusedBy === 'no-genotype-signal'
    || result.stage?.refusedBy === 'second-genome' ? result.stage.refusedBy : 'proxy'
  const stageMeasured = { heterozygosity: profile.hetRate, callRate: profile.callRate }
  // THE ANSWER KEYS, reported whatever uses them, so the parent-free and parent-based readings can be
  // set side by side on any cell that has both.
  result.answerKeys = {
    dropIn: {
      parents: acc.sameHom && acc.sameHom.n >= MIN_SAME_HOM ? dropIn(acc.sameHom) : NaN,
      panel: acc.panelSameHom && acc.panelSameHom.n >= MIN_SAME_HOM ? dropIn(acc.panelSameHom) : NaN,
    },
    information: {
      parents: ((rs: number[]) => (rs.length ? Math.min(...rs) : NaN))([acc.informationPat, acc.informationMat]
        .map((t) => (t ? informationRatio(t) : NaN)).filter((r) => Number.isFinite(r))),
      panel: acc.informationPanel ? informationRatio(acc.informationPanel) : NaN,
    },
    heterozygositySpread: spread,
    normalFraction: normal,
  }
  // DOES THIS CELL'S GENOTYPE CHANNEL WORK? The stage refusal fires on heterozygosity, which drop-in
  // raises without moving which parent's allele is present: on the cells it refuses in GSE148488 the
  // allele calls agree with the published karyotype at 133 of 138. The information ratio asks the
  // question directly, against a measured parent, and is the test for a channel carrying nothing.
  // Where it says the calls carry information, the allele claims stand through the stage refusal.
  const genotypeSignal = (() => {
    const rs = [acc.informationPat, acc.informationMat]
      .map((t) => (t ? informationRatio(t) : NaN)).filter((r) => Number.isFinite(r))
    return rs.length ? Math.min(...rs) < NO_INFORMATION_MIN : false
  })()
  // ONLY THE HETEROZYGOSITY CEILING. That branch is the proxy: drop-in trips it while the alleles
  // are fine. A declared failure is the operator's instruction and a collapsed call rate is a
  // different question, so both still clear every claim.
  const alleleClaims = new Set<string>(
    genotypeSignal && result.stage?.refusedBy === 'heterozygosity-ceiling' ? ['lostParent', 'mixed'] : [])
  if (uninterpretable) {
    log('WARN', 'no chromosomal changes are reported: this array did not resolve to a '
      + `stage. ${result.stage.why}`)
    // RECORDED BEFORE THEY ARE CLEARED, so the census can score what this gate threw away. Every
    // field cleared below has its record here; refusals.check.ts pins the two lists to each other.
    for (const c of result.chroms) {
      for (const claim of ['aneuploidy', 'lostParent', 'aneuploidyParent', 'mixed', 'alleleGain'] as const) {
        if (c[claim] === undefined || alleleClaims.has(claim)) continue
        result.refusals!.push({ claim, chrom: c.chrom, gate: stageGate, cause: stageCause,
          measured: stageMeasured, withheld: c[claim] })
      }
    }
    for (const c of result.chroms) {
      c.aneuploidy = undefined
      if (!alleleClaims.has('lostParent')) c.lostParent = undefined
      c.aneuploidyParent = undefined
      // The mixed signal and the allele gain are findings of the same measurement the array has
      // just been refused for, and leave with the rest.
      if (!alleleClaims.has('mixed')) c.mixed = undefined
      c.alleleGain = undefined
    }
  }
  // ONE GENOME WITH HEAVY DROP-IN: THE ALLELES ARE READ, THE INTENSITY IS NOT. The stage removed the
  // drop-in and found one genome (stage.ts), so the allele claims stand: drop-in adds heterozygous
  // calls without adding an absent parent's allele. Copy number is a different question. It is read
  // against this array's own per-chromosome median, which needs a normal majority, and on arrays like
  // this the test that would show one cannot: the localised-excess spread reads 0.07 to 0.18 on every
  // array over the ceiling in GSE148488, chaotic or not, against 0.50. So it is withheld, recorded,
  // under its physical cause, until a normal majority can be shown from the alleles.
  const dropInCell = result.stage?.dropInCorrected !== undefined
  const dropInMeasured = { dropIn: result.stage?.dropInCorrected ?? NaN, heterozygositySpread: spread,
    normalFraction: normal }
  // A NORMAL MAJORITY, SHOWN FROM THE ALLELES, IS WHAT COPY NUMBER NEEDS, on every array and not only
  // on a drop-in one. The intensity is read against this array's own per-chromosome median, which is a
  // normal chromosome only when most are. The disomic fraction asks that of the alleles directly: 0.000
  // to 0.727 on cells published chaotic, 0.864 to 1.000 on the rest, 0 on arrays whose calls are noise.
  // Where it cannot be computed, a drop-in array has shown none and every other array is read as before.
  // With no alleles to ask, the majority is the one intensityNull.majorityScale finds or fails to.
  const normalMajority = Number.isFinite(normal) ? normal >= NORMAL_MAJORITY_MIN
    : !allelesReadable || !dropInCell
  const intensityWithheld = !uninterpretable && !normalMajority
  if (intensityWithheld) {
    for (const c of result.chroms) {
      for (const claim of ['aneuploidy', 'aneuploidyParent', 'alleleGain'] as const) {
        if (c[claim] === undefined) continue
        result.refusals!.push({ claim, chrom: c.chrom, gate: 'normal-majority',
          cause: 'no-intensity-reference', measured: dropInMeasured, withheld: c[claim] })
      }
      c.aneuploidy = undefined
      c.aneuploidyParent = undefined
      c.alleleGain = undefined
    }
  }
  // THE OTHER PARENT, READ OFF THE OTHER DRAWS OF HER THIS RUN HOLDS (motherFromRun.ts), where one parent
  // is loaded and shown, and the cell has a two-copy majority to read intensity against.
  if (!two && parentShown && shareCalls.size && !intensityWithheld
    && (acc.panel as FounderPanel | null)?.freqB) {
    readMotherFromRun({ acc, pat, sibs, shareCalls, result, soloRole, otherRole, log })
  }
  // WHETHER THE CALLS THEMSELVES SAY ANYTHING, which is the one thing the stage refusal cannot
  // tell a reader. The ceiling it fires on cannot separate a mixture from a destroyed array; the
  // information ratio can, because a mixture reads below 1 and noise reads at it.
  //
  // THE MOST FAVOURABLE PARENT DECIDES. Scored against each measured parent and judged on the
  // LOWEST reading, so a second array that is itself poor cannot put this on a sample whose own
  // calls are fine.
  //
  // REASONS ONLY. It goes to the array-level limits, so no chromosome moves between called and
  // refused and every per-chromosome reason stands underneath unchanged.
  {
    const scored = ([[acc.informationPat, soloRole], [acc.informationMat, otherRole]] as const)
      .flatMap(([tally, role]) => (tally ? [{ r: informationRatio(tally), role }] : []))
      .filter(({ r }) => Number.isFinite(r))
    if (scored.length) {
      const best = scored.reduce((a, b) => (b.r < a.r ? b : a))
      const why = whyNoInformation(best.r, `the ${best.role} array`)
      if (why) {
        const text = `${why[0].toUpperCase()}${why.slice(1)}.`
        result.limits.push(text)
        log('WARN', text)
      }
    }
  }

  // IS THE ARRAY IN THE PARENTAL SLOT THE KIND OF PERSON THAT SLOT SAYS IT IS?
  //
  // Checked on every parental array loaded, against the role the operator declared for it. The
  // mistake this catches makes no error anywhere else: a valid, high-quality array of a real
  // person, put in the wrong slot, produces arithmetic that is entirely correct under a label that
  // is entirely wrong. See parentSanity.ts for the measurement and for why silence is not a
  // conflict.
  result.parentSanity = reconcileParentSex(soloRole, pat.sex)
  const matSanity = mat ? reconcileParentSex(otherRole, mat.sex) : undefined
  if (matSanity) result.parentSanityOther = matSanity
  // EITHER slot being wrong poisons the naming, because with both arrays loaded the answer comes
  // from whichever side reported its own copy absent, and a swap makes that side the wrong one.
  // AND IS THE ARRAY IN THE PARENTAL SLOT A PARENT OF THIS SAMPLE AT ALL?
  //
  // THE SEX CHECK ABOVE CATCHES THE WRONG KIND OF PERSON. This catches the wrong person. Measured
  // on 8 children of the published trios, each scored against five wrong files, the tool named a
  // parent on 2 runs against an unrelated adult and 4 against a full sibling, at the confidence of
  // a correct run and with nothing on the page to say otherwise.
  //
  // NOTHING HERE IS NEW ARITHMETIC. `om link` has separated a parent from a sibling and from an
  // unrelated adult since the segmental test was measured: over the published trios, 181 verified
  // parent-child pairs reach a window spread of at most 0.0238 and 166 constructed full-sibling
  // pairs at least 0.0312. Those functions lived in the command-line file, so the surface most lab
  // staff use could not run them. They now live in relatedness.ts and both surfaces call them.
  //
  // A REFUSAL IS ONLY TAKEN ON A DECISIVE VERDICT. `firstDegree` means the segmental test had no
  // power here, not that the pair failed it, and refusing on it would refuse whenever the panel is
  // thin. Only `unrelated` and `sibling` are grounds to withhold, and `oneParent` is left alone
  // because a haploid product genuinely reads that way against its own parent.
  //
  // WHAT THAT COSTS, MEASURED RATHER THAN ASSUMED. On amplified single-cell material the segmental
  // test usually has no power: a real trophectoderm biopsy against its own sperm donor reads
  // opposite homozygotes 0.00680 over 521,965 markers, a window spread of 0.0250 over 2,719
  // windows, and comes back `firstDegree`. So on exactly the material this tool is for, the
  // UNRELATED half of the test does the work and the SIBLING half often cannot. An unrelated adult
  // is caught by the rate, which is decisive; a sibling of the same embryo series may not be, and
  // is then caught only if the sexes happen to differ. That gap is real and is not papered over.
  //
  // IT COSTS ABOUT THREE AND A HALF SECONDS PER SAMPLE on the slowest machine this was measured on,
  // dominated by sorting half a million markers into genomic order. Paid once per sample, against
  // a report that would otherwise name a parent that is not one.
  //
  // NOT ON AN ARRAY THE TOOL HAS ALREADY REFUSED. A failed array reports nothing, so a relationship
  // verdict on it describes a genome that is not being read, at three and a half seconds a sample.
  // It is also where the statistic is least meaningful: the arrays that fail this way read 55 to 62
  // percent heterozygous, which is a mixed or contaminated sample rather than a genome, and they
  // read as unrelated to everyone including themselves.
  if (uninterpretable) {
    result.refusals!.push({ claim: 'relatedness', gate: stageGate, cause: stageCause,
      measured: stageMeasured })
  }
  const rel = uninterpretable || pat.absent ? null : relate(
    { gt: pat.gt, pos: markerPos },
    { gt: myGt as Map<string, AB>, pos: markerPos },
    OPPOSITE_HOM_MAX,
  )
  if (rel) {
    // The verdict is withheld where the material cannot separate a parent from a stranger (see
    // RELATEDNESS_SEPARABLE_MATERIALS). The rate and marker count are measurements and stay.
    const assessable = relatednessAssessable(result.stage?.stage)
    result.relationship = {
      verdict: assessable ? rel.relationship : REL.notAssessable,
      oppositeHomRate: rel.opp.rate,
      oppositeHomMarkers: rel.opp.n,
      windowSpread: rel.win.spread,
      windows: rel.win.windows,
    }
  }
  // A UNIPARENTAL GENOME IS NOT AN OPERATOR ERROR, AND IT LOOKS EXACTLY LIKE ONE HERE. A
  // gynogenote carries no paternal contribution at all, so against the true sperm donor its
  // opposite-homozygote rate is that of two unrelated people. That is the FINDING, not a wrong
  // file, and it is the one case this tool most needs to be able to report: withholding on it
  // silenced "the paternal copy is absent" on a genome whose paternal copy is absent.
  //
  // `relate` already has a verdict for this and cannot reach it here: it tests ploidy only AFTER
  // the opposite-homozygote gate, which a uniparental sample fails first. So the sample's own
  // zygosity, which is measured independently of any parental array, is what separates the two.
  // AND IT DOES NOT WITHHOLD, because on this material no threshold exists that could.
  //
  // WIRING IT AS A GATE WITHHELD EVERY NAMED PARENT ON EIGHT OF EIGHT CORRECT RUNS. The first real
  // single cell it saw, a blastomere against its genetically confirmed father, read 0.0276 opposite
  // homozygotes over 524,821 markers and came back `unrelated`. Amplification drops alleles; a
  // heterozygous marker that loses one allele is CALLED homozygous, and if what it kept is the
  // allele the parent lacks, the pair reads as opposite homozygotes. The statistic is inflated by
  // exactly what single-cell material has most of, in the direction that makes a parent a stranger.
  // The 0.020 gate is sound where it was measured, on bulk adult DNA for `om link`.
  //
  // A PER-MATERIAL THRESHOLD WAS THEN MEASURED RATHER THAN ASSUMED, in
  // audit/relatedness-by-material.ts: every usable trio of GSE148488 against its confirmed father,
  // its confirmed mother, a second array of the father as a control, and an unrelated donor. Of 76
  // trios the tool refuses 36 outright as `failed`, so those never reach this code and are
  // excluded. On the 40 that remain:
  //
  //   material        true parents      an unrelated donor   separable?
  //   esc-line        0.0013 to 0.0110  0.0531 to 0.0562     yes, a five-fold gap, n=16
  //   esc-single      0.0077 to 0.0311  0.0592 to 0.0803     yes, but n=4
  //   blastomere      0.0183 to 0.0658  0.0693 to 0.1161     a 5 percent margin, n=8
  //   trophectoderm   0.0035 to 0.0831  0.0567 to 0.0954     NO, they overlap by 0.0263
  //
  // The shipped 0.020 fires on 8 of 8 true blastomere fathers, 2 of 4 esc-single and 4 of 12
  // trophectoderm. Only esc-line passes it. And trophectoderm, the commonest PGT material, has no
  // threshold at all: a true parent reaches 0.0831 where a stranger starts at 0.0567.
  //
  // So this cannot gate a parental call on the material this tool is for, and pretending otherwise
  // would refuse correct runs. Revisiting it needs more trios per material, trophectoderm most.
  // The verdict is still computed and still reported, because a reader deciding whether to trust
  // a run should see what the loaded array measures against the sample. It just does not gate
  // anything until a per-material threshold exists. See audit/relatedness-by-material.ts.

  // ONE REASON STRING, because three different checks can withhold and every message that quotes
  // one of them has to quote the one that actually fired. Reading the sex conflict when the
  // relationship test was what refused threw `Cannot read properties of undefined`, on a run whose
  // only fault was a wrong file.
  const sexConflict = result.parentSanity.conflict ? result.parentSanity
    : matSanity?.conflict ? matSanity : undefined
  // One parent's genome absent is a finding, a uniparental genome or the other parent's product;
  // both absent is two arrays that are not this sample's parents.
  const neitherParent = patGenome.verdict === 'none' && matGenome?.verdict === 'none'
  const withholdReason = sexConflict
    ? `the parental array does not match the role it was loaded under. ${sexConflict.why}`
    : neitherParent
      ? 'neither loaded parental array is a parent of this sample: no chromosome read carries the '
        + 'genome of either.'
      : ''
  /** The same fact in one clause, for the per-row strings. See the alert for the whole reason. */
  const shortWithholdReason = sexConflict
    ? `the array loaded as the ${sexConflict.declared} parent is not the sex that slot requires.`
    : neitherParent ? 'neither loaded parental array carries its genome into this sample.' : ''
  const parentSlotWrong = !!withholdReason
  if (parentSlotWrong) {
    log('WARN', `NO PARENT IS NAMED ON THIS RUN: ${withholdReason}`)
  }
  /**
   * Strip the parent off a dosage row when the parental slot is wrong.
   *
   * THE VERDICT CHANGES AND NOT ONLY THE FIELD. `originRows` and the report both fall back to
   * naming a parent from the VERDICT when the explicit field is absent, so clearing `parent` alone
   * would withhold nothing: `loaded-parent` still reads as the loaded parent's name. The row keeps
   * its location, its shift and its evidence, because those do not depend on which slot the
   * parental array went into.
   */
  const withheldIfSlotWrong = <T extends { verdict: string; why: string; parent?: unknown }>(
    call: T,
  ): T => (parentSlotWrong
    ? {
      ...call,
      verdict: 'array-excluded',
      parent: null,
      // SHORT, because the full reason is on the run's own alert and this string is repeated on
      // every row. Printing a seven-hundred-character paragraph beside each event buried the
      // measurement it was supposed to qualify.
      why: `no parent is named on this run: ${shortWithholdReason} See the notice on this sample `
        + `for what to do. The measurement itself stands: ${call.why}`,
    }
    : call)
  // Segments, after the per-chromosome verdicts, because a chromosome whose calls are not
  // measuring it must not be scanned either: the same broken calls would produce a
  // confident segment inside it.
  const measured = new Set(result.chroms
    .filter((c) => c.verdict !== 'not_measured').map((c) => c.chrom))
  // Copy number first: a chromosome already called whole is not scanned for a segment
  // inside it, since that is the same event described twice.
  // A chromosome the alleles read as missing one parent is whole as well: a segment found inside
  // it would be the same absence described twice.
  const whole = new Set(result.chroms.filter((c) => c.aneuploidy || c.lostParent).map((c) => c.chrom))
  const lostByAlleles = new Set(result.chroms.filter((c) => c.lostParent).map((c) => c.chrom))
  // LOSSES ONLY, for the run-of-homozygosity filter below. A chromosome at ONE copy reads
  // homozygous end to end because it is hemizygous, and that artefact has to go. A chromosome at
  // THREE copies has no such excuse: a trisomy that IS homozygous end to end is an isodisomic
  // trisomy, which is the trisomy-rescue and imprinting mechanism and a separate fact from the
  // gain. Filtering on `whole` deleted it, and the gain row left behind does not carry it.
  const lost = new Set(result.chroms.filter((c) => c.aneuploidy === 'loss').map((c) => c.chrom))
  const noCall = new Map<string, [number, number]>()
  for (const [c, ms] of cnByChrom) {
    noCall.set(c, [ms.length, ms.filter((m) => !m.called).length])
  }
  const allCn = [...cnByChrom.values()].flat()
  // The segment scans read against the chromosomes the whole-chromosome test does, those sharing
  // the array's commonest copy state (parentage.copyReference): absent chromosomes would otherwise
  // pull the baseline under the present ones.
  const referenceChroms = new Set(result.copyReference ?? [])
  const referenceCn = referenceChroms.size >= 8
    ? new Map([...cnByChrom].filter(([c]) => referenceChroms.has(c))) : cnByChrom
  const lrrAll = [...referenceCn.values()].flat()
    .map((m) => m.log2R).filter((x): x is number => x !== null).sort((a, b) => a - b)
  const genomeLrr = lrrAll.length ? lrrAll[lrrAll.length >> 1] : 0
  // Genome call rate from the markers themselves. StageCall does not carry one, and the stage's
  // own inference took it as an input rather than storing it.
  const genomeCallRate = allCn.length ? allCn.filter((m) => m.called).length / allCn.length : NaN
  // WHAT A FILE WITHOUT INTENSITY CANNOT SHOW, stated rather than left to read as a clean genome.
  if (!lrrAll.length) {
    result.limits.push(
      'This file carries genotype calls and no intensity (log2R). A gained chromosome or segment '
      + 'cannot be seen at all, because a third copy leaves the calls unchanged. A lost one is seen '
      + 'only as heterozygosity gone, which is what isodisomy looks like too, so those findings are '
      + 'reported as either. A chromosome absent altogether still reads through its call rate.',
    )
  }
  // --- what the sex chromosomes are ---------------------------------------------------
  //
  // The profile's sex call reads chrX heterozygosity, which is the same on a 46,XY male and on a
  // 45,X female: one X, no heterozygous calls outside the pseudoautosomal region. Copy number and
  // chrY together separate them, where the panel carries chrY at all.
  //
  // Haploid material is not asked, because one X is what it is supposed to carry, and neither is
  // an array whose stage did not resolve: this reads the same intensity every other call rests on.
  // With no parent, two copies must also be shown past drop-in: see sexMaterial.ts.
  const xChrom = result.chroms.find((c) => c.chrom === 'X' || c.chrom === '23')
  const material = sexMaterial({
    stage: result.stage, hetRate: profile.hetRate, dropIn: panelDropIn, parentLoaded: !pat.absent,
  })
  if (pat.absent) result.sexMaterial = material
  result.sexChromosome = callSexChromosomes({
    xShift: xChrom?.lrrShift ?? NaN,
    xZ: xChrom && Number.isFinite(xChrom.lrrZ) ? xChrom.lrrZ : undefined,
    y: acc.sex,
    bulk: bulkMaterial,
    diploid: material.diploid,
    material,
    twoInsideNoise: pat.absent,
    hetSex: profile.sex,
    hetRatio: profile.chrXHetRatio,
    // Every file scored here is the proven chip: the readers refuse any other (platform.ts).
    chipGenotypesY: PROVEN_CHIP_GENOTYPES_Y,
    xCallFraction: xChrom?.callFraction,
    // A row needs 200 chrX markers, its level 200 with intensity. The pseudoautosomal region is
    // left out only where the tally knows the assembly, which comes from the parent.
    xMissing: !xChrom && !pat.absent
      ? 'the loaded parent gives fewer than 200 chromosome X markers to read against, so no chrX '
        + 'row was built'
      : `${!xChrom ? 'this file has fewer than 200 chromosome X markers'
      : 'fewer than 200 chromosome X markers carry intensity on this file'}${t.build
      ? ' outside the pseudoautosomal region' : ''}`,
  })
  if (result.sexChromosome.aneuploid) {
    log('WARN', `sex chromosomes: ${result.sexChromosome.constitution}. `
      + `${result.sexChromosome.why}.`)
  } else if (result.sexChromosome.xCopies === 1
    && result.sexChromosome.constitution === 'unresolved') {
    const yWhy = sexCall(acc.sex).why
    result.limits.push(
      'This sample carries ONE copy of chromosome X. Whether that is an ordinary male karyotype '
      + `or a monosomy X cannot be said from this file, ${yWhy ? `because ${yWhy}`
        : 'which carries too few chromosome Y markers to test for a Y'}. Both are consistent with `
      + 'everything measured here.',
    )
  }
  // The intensity scan is withheld on a chaotic genome: it normalises against this sample's own
  // per-chromosome median and there is no normal majority here to be that median.
  //
  // TWO NOMINATION CHANNELS, because a gain and a loss do not announce themselves the same way.
  // scanCopyNumber proposes a window where the array STOPPED GENOTYPING, which an extra copy does
  // not do: measured inside hand-called gains the scan found, the no-call rate is 74% to 99%
  // against an array baseline near 40%, which is amplification damage, while inside the gains it
  // never proposed it is 40% to 58%, barely over baseline. Of 26 hand-called gains that channel
  // loses 8 to nomination and 0 to the intensity gate, so a lower gate could never have reached
  // them. scanCopyGain proposes on intensity instead and reaches 25 of the 26, and 0 of 67 losses,
  // the second by construction since a lost copy cannot read intensity-high.
  const gainProfile = highLrrProfile(referenceCn)
  // A chromosome above the majority and not reported as gained is not reported in pieces either.
  const unresolved = new Set(result.copyUnresolved ?? [])
  // THE COPY-NUMBER SCAN READS EVERY CHROMOSOME INTENSITY HAS NOT ALREADY CALLED WHOLE, whatever the
  // alleles say. It used to skip one the alleles name as missing a parent too, and those names are what
  // differ with how many parents are loaded, so the scan itself depended on the parents. A segment that
  // restates an allele-named loss across most of the chromosome is dropped below, at reporting, not here.
  // A WHOLE CALL NEEDS EACH ARM: see wholeArms.ts. Read before 3.4a, so a chromosome withdrawn here goes
  // to both piece scans and is never read again as a whole loss. A chromosome the alleles also read
  // whole is the alleles' call as well and is left as it is.
  const withdrawn: { chrom: string; kind: 'loss' | 'gain'; read: ArmReading }[] = []
  if (result.copyNull && profile.build?.build === 'GRCh37') {
    for (const c of result.chroms) {
      if (!c.aneuploidy || !isAutosome(c.chrom) || c.lostParent
        || c.alleleGain?.verdict === 'gain' || c.alleleGain?.verdict === 'imbalance') continue
      const ms = cnByChrom.get(c.chrom)
      if (!ms?.length) continue
      const read = readArms(ms, c.aneuploidy, c.lrrShift, result.copyNull)
      if (read.whole) continue
      withdrawn.push({ chrom: c.chrom, kind: c.aneuploidy, read })
      c.aneuploidy = undefined
      c.aneuploidyParent = undefined
    }
  }
  const wholeByIntensity = new Set(result.chroms.filter((c) => c.aneuploidy).map((c) => c.chrom))
  // A CHROMOSOME CALLED LOST WHOLE IS READ AGAIN IN PIECES: see partLost.ts. Where its loss pieces
  // account for the shift the pieces are reported and the whole call is not; where the rest is lost
  // as well the whole call stands, with any piece a further copy down beside it.
  const partLost = new Map<string, PartLost>()
  if (result.copyNull) {
    for (const c of result.chroms) {
      if (c.aneuploidy !== 'loss' || !isAutosome(c.chrom)) continue
      const ms = cnByChrom.get(c.chrom)
      if (ms) partLost.set(c.chrom, readLostChromosome(ms, externalNull(noCall, c.chrom), genomeLrr, result.copyNull))
    }
  }
  const scannedCopyRaw = [...[...partLost.values()].flatMap((p) => p.pieces), ...[...cnByChrom]
    .filter(([c]) => !wholeByIntensity.has(c))
    .flatMap(([c, ms]) => {
      const byNoCall = scanCopyNumber(ms, externalNull(noCall, c), genomeLrr)
      return [
        ...byNoCall,
        // `byNoCall` is passed so a region both channels see is reported once, not twice.
        ...(unresolved.has(c) ? []
          : scanCopyGain(ms, externalNull(gainProfile.byChrom, c), genomeLrr, gainProfile, byNoCall,
            result.copyOffPole ?? NaN)),
      ]
    })]
  // A WITHDRAWN WHOLE CALL IS NOT SILENCE: where the scans place nothing of its kind on the arm that
  // carries it, that arm is the call.
  for (const w of withdrawn) {
    if (!w.read.arm) continue
    const sg = armSegment(cnByChrom.get(w.chrom)!, w.read.arm, w.kind === 'gain' ? 'copy-gain' : 'copy-loss')
    if (!scannedCopyRaw.some((s) => s.chrom === sg.chrom && s.kind === sg.kind
      && s.startBp <= sg.endBp && s.endBp >= sg.startBp)) scannedCopyRaw.push(sg)
  }
  for (const [chrom, p] of partLost) {
    if (p.whole) continue
    const c = result.chroms.find((x) => x.chrom === chrom)!
    c.aneuploidy = undefined
    c.aneuploidyParent = undefined
  }
  // ONE EVENT, TWO LINES OF EVIDENCE. A copy loss covering 80% or more of a chromosome the alleles already
  // name as missing a parent's copy is that same loss, so the chromosome row carries it and the segment
  // is not printed a second time. A partial segment inside it is a different claim and stands.
  const spanOf = (c: string): number => {
    const ms = cnByChrom.get(c) ?? []
    let lo = Infinity; let hi = -Infinity
    for (const m of ms) { if (m.pos < lo) lo = m.pos; if (m.pos > hi) hi = m.pos }
    return hi > lo ? hi - lo : NaN
  }
  const alleleLost = new Set(result.chroms.filter((c) => c.lostParent).map((c) => c.chrom))
  const scannedCopy = scannedCopyRaw.filter((sg) => !(sg.kind === 'copy-loss' && alleleLost.has(sg.chrom)
    && (sg.endBp - sg.startBp) >= RESTATED_SHARE * spanOf(sg.chrom)))
  // RUN EVEN WHERE WITHHELD, and recorded rather than reported, so the census can score what the
  // withholding cost. Nothing below reads `scannedCopy` except through `copy`.
  const copyWithheld = uninterpretable || chaoticGenome || intensityWithheld
  // On a drop-in array a gained segment is printed only where the alleles show a normal majority: it
  // is the one copy claim resting on a few megabases of intensity alone.
  const gainSegmentsWithheld = !copyWithheld && dropInCell && !allelesReadable
  const copy = copyWithheld ? []
    : gainSegmentsWithheld ? scannedCopy.filter((sg) => sg.kind !== 'copy-gain') : scannedCopy
  if (gainSegmentsWithheld) {
    for (const sg of scannedCopy) {
      if (sg.kind !== 'copy-gain') continue
      result.refusals!.push({ claim: 'segment', chrom: sg.chrom, gate: 'drop-in.no-alleles',
        cause: 'no-intensity-reference', measured: dropInMeasured, withheld: sg })
    }
  }
  if (copyWithheld) {
    for (const sg of scannedCopy) {
      result.refusals!.push(intensityWithheld
        ? { claim: 'segment', chrom: sg.chrom, gate: 'normal-majority', cause: 'no-intensity-reference',
          measured: dropInMeasured, withheld: sg }
        : { claim: 'segment', chrom: sg.chrom,
          gate: chaoticGenome ? 'copy-scan.chaotic-genome' : stageGate, cause: chaoticGenome ? 'proxy' : stageCause,
          measured: { ...stageMeasured, heterozygositySpread: spread }, withheld: sg })
    }
  }
  // THE LIKELIHOOD SCANS, otherParentLoss.ts. Each loaded parent's markers on the chromosomes the
  // one-parent copy state says carry both parents, on a genome read as two parents'. Those
  // chromosomes are the reference the array's own rates are fitted on, and the ones scanned.
  const twoParents = knownBiparental(result.zygosity)
  const carryingBoth = (by: Map<string, OtherMarker[]> | null): Map<string, OtherMarker[]> => new Map(
    [...(by ?? [])]
      .filter(([c]) => shareCalls.get(c)?.state === 'disomic' && measured.has(c) && !whole.has(c))
      .map(([c, ms]) => [c, [...ms].sort((x, y) => x.pos - y.pos)] as const))
  const ratesOf = (ms: Map<string, OtherMarker[]>) => (
    ms.size >= BIPARENTAL_AUTOSOMES_MIN ? otherRates([...ms.values()]) : null)
  const matShown = matGenome !== undefined && (matGenome.verdict === 'carries' || matGenome.verdict === 'unread')
  const patMarkers = parentShown && twoParents ? carryingBoth(acc.otherByChrom) : new Map<string, OtherMarker[]>()
  const matMarkers = mat && matShown && twoParents ? carryingBoth(acc.otherMatByChrom) : new Map<string, OtherMarker[]>()
  const patRates = ratesOf(patMarkers)
  const matRates = ratesOf(matMarkers)
  // Mean log R across a run less the median of the rest of its chromosome.
  const lrrShiftOf = (c: string, startBp: number, endBp: number): number => {
    const inside: number[] = []
    const rest: number[] = []
    for (const m of cnByChrom.get(c) ?? []) {
      if (m.log2R !== null) (m.pos >= startBp && m.pos <= endBp ? inside : rest).push(m.log2R)
    }
    return inside.length && rest.length ? inside.reduce((x, y) => x + y, 0) / inside.length - median(rest) : NaN
  }
  // How far a run's mean log R sits from stretches of the same number of markers on the other
  // chromosomes that carry both parents, in units of their spread. See intensityNull.ts.
  const lrrSorted = new Map<string, number[]>()
  const lrrOf = (c: string): number[] => {
    let v = lrrSorted.get(c)
    if (!v) {
      v = [...(cnByChrom.get(c) ?? [])].sort((x, y) => x.pos - y.pos)
        .map((m) => m.log2R).filter((x): x is number => x !== null)
      lrrSorted.set(c, v)
    }
    return v
  }
  const mean = (xs: readonly number[], from = 0, to = xs.length): number => {
    let sum = 0
    for (let i = from; i < to; i += 1) sum += xs[i]
    return sum / (to - from)
  }
  const lrrZOf = (run: Segment, reference: Iterable<string>): number | undefined => {
    const inside = (cnByChrom.get(run.chrom) ?? [])
      .filter((m) => m.log2R !== null && m.pos >= run.startBp && m.pos <= run.endBp).map((m) => m.log2R as number)
    if (!inside.length) return undefined
    const units: number[] = []
    for (const c of reference) {
      if (c === run.chrom) continue
      const v = lrrOf(c)
      for (let i = 0; i + inside.length <= v.length; i += inside.length) units.push(mean(v, i, i + inside.length))
    }
    return calibratedZ(mean(inside), nullScale(units))
  }
  // ANOTHER ARRAY OF THE SAME EMBRYO carries the same two parental copies, so what the parents share
  // reads as missing in it too and a loss in this array alone does not. Each other array of the run
  // is read against the loaded parent at this array's markers, kept if it reads as the same genome,
  // and asked whether it shows the other parent's alleles across a run.
  let sisters: { byChrom: Map<string, OtherMarker[]>; rates: OtherRates }[] | undefined
  const shownBySister = (run: Segment): boolean => {
    sisters ??= sibs.flatMap((sib) => {
      const byChrom = new Map<string, OtherMarker[]>()
      const where: { q: number; seen: boolean }[] = []
      for (const [c, ms] of patMarkers) {
        const theirs: OtherMarker[] = []
        for (const m of ms) {
          const g = sib.get(m.id)
          const p = pat.gt.get(m.id)
          if (!g || g === 'NC' || !p) continue
          const seen = m.q === null ? g === 'AB' : g !== p
          theirs.push({ id: m.id, pos: m.pos, q: m.q, seen, opp: m.q !== null && g !== 'AB' && g !== p })
          if (m.q !== null && m.seen) where.push({ q: m.q, seen })
        }
        byChrom.set(c, theirs)
      }
      const rates = otherRates([...byChrom.values()])
      return rates && sameGenomeShare(where, rates) >= SAME_GENOME_SHARE ? [{ byChrom, rates }] : []
    })
    return sisters.some((s) => presentScore(s.byChrom.get(run.chrom) ?? [], s.rates, run.startBp, run.endBp) >= OTHER_LRT)
  }
  // AN EXTRA COPY READS LIKE A LOSS IN THE CALLS: it outnumbers the other parent's allele, which
  // then drops out of them. A lost copy cannot raise intensity, so a run is refused where a gain
  // segment covers it or its own intensity stands a copy state above the rest of its chromosome.
  // The gain scan alone is not enough: it asks a larger rise than one extra copy gives a single cell.
  const underGain = (run: Segment, shift: number): boolean => shift >= STATE_GAP || coveredShare(
    run.startBp, run.endBp,
    scannedCopy.filter((sg) => sg.chrom === run.chrom && sg.kind === 'copy-gain')) >= RECIPROCAL_OVERLAP
  // WITH ONE PARENT, THE OTHER PARENT'S COPY. The absence scan below reads the loaded parent's own
  // alleles, so a stretch missing the other parent's copy reads present throughout. Not where
  // intensity reads a gain.
  const otherSegments: Segment[] = []
  if (acc.otherByChrom && !mat && parentShown && twoParents) {
    if (patRates) {
      // A run under OTHER_MIN_SPAN_BP reads the same as a haplotype the parents share, in the
      // alleles. It is a finding where intensity falls across it as well, by the calibrated test,
      // since a shared haplotype is at two copies; or where another array of the same embryo shows
      // the other parent's alleles across it. Otherwise it is withheld, and its place recorded
      // and stated; the other parent's array is what would settle it.
      const short: Segment[] = []
      for (const [c, ms] of patMarkers) {
        for (const run of scanOtherParent(c, ms, patRates, MIN_SEGMENT_SPAN_BP)) {
          const lrrShift = lrrShiftOf(c, run.startBp, run.endBp)
          if (underGain(run, lrrShift)) continue
          const sg: Segment = { ...run, absentParent: otherRole, lrrShift }
          if (sg.spanBp >= OTHER_MIN_SPAN_BP) { otherSegments.push(sg); continue }
          sg.lrrZ = lrrZOf(sg, patMarkers.keys())
          sg.shortBy = sg.lrrZ !== undefined && sg.lrrZ <= -Z_WINDOW ? 'intensity'
            : shownBySister(sg) ? 'sister' : undefined
          ;(sg.shortBy ? otherSegments : short).push(sg)
        }
      }
      for (const sg of otherSegments) {
        if (!sg.shortBy) continue
        const mb = (bp: number): string => (bp / 1e6).toFixed(1)
        result.notes.push(
          `chr${sg.chrom} ${mb(sg.startBp)}-${mb(sg.endBp)} Mb is under ${OTHER_MIN_SPAN_BP / 1e6} Mb, where a stretch both `
          + `parents share reads the same in the alleles. It is reported as a loss of the ${otherRole} copy because `
          + (sg.shortBy === 'intensity' ? 'intensity falls across it as well.'
            : `another array of this run, read as the same embryo, shows the ${otherRole} alleles across it.`))
      }
      for (const sg of short) {
        result.refusals!.push({ claim: 'segment', chrom: sg.chrom, gate: 'other-parent.under-span',
          cause: 'parent-not-loaded',
          measured: { spanBp: sg.spanBp, score: sg.score, markers: sg.markers, ...(sg.lrrZ === undefined ? {} : { lrrZ: sg.lrrZ }) },
          withheld: sg })
      }
      if (short.length) {
        const mb = (bp: number): string => (bp / 1e6).toFixed(1)
        result.limits.push(
          `No ${otherRole} allele is seen on ${short.map((sg) => `chr${sg.chrom} ${mb(sg.startBp)}-${mb(sg.endBp)} Mb`).join(', ')}. `
          + `Under ${OTHER_MIN_SPAN_BP / 1e6} Mb that is reported as a loss of the ${otherRole} copy only `
          + 'with a second line, and here there is none: intensity does not measurably fall across it, '
          + 'and no other array of this run that reads as the same embryo shows those alleles there. '
          + 'A stretch where both parents carry the same DNA reads the same in the alleles with one '
          + 'parent loaded. The other parent\'s array, or another cell of the embryo, settles it.',
        )
      }
      log('SCAN', `${otherRole} copy scan: ${patMarkers.size} chromosomes carrying both parents, `
        + `clustering ${patRates.inflation.toFixed(1)}, ${otherSegments.length} run(s)`
        + `${short.length ? `, ${short.length} under ${OTHER_MIN_SPAN_BP / 1e6} Mb withheld` : ''}`)
    } else {
      result.limits.push(
        `A stretch missing the ${otherRole} copy was not looked for. With one parent loaded it is read `
        + 'against this array\'s own chromosomes that carry both parents, and too few of them do, or '
        + 'they carry too few markers, to give a rate to read it against.',
      )
    }
  }
  // The window scan, per loaded parent, against that parent's own homozygous markers. A chromosome
  // the alleles read as missing one parent is not scanned for an absence segment inside it: that
  // would be the same absence described twice.
  const windowPat = [...absenceByChrom].filter(([c]) => measured.has(c) && !lostByAlleles.has(c))
    .flatMap(([c, ms]) => scanChromosome(ms, externalNull(t.byChrom, c))
      .map((sg): Segment => (mat ? { ...sg, absentParent: soloRole } : sg)))
  // Without the second run a deletion of the other parent's copy is invisible: the first parent's
  // allele reads present at every marker of it.
  const windowMat = mat ? [...matAbsenceByChrom].filter(([c]) => measured.has(c) && !lostByAlleles.has(c))
    .flatMap(([c, ms]) => scanChromosome(ms, externalNull(matAbsenceTally, c))
      .map((sg): Segment => ({ ...sg, absentParent: otherRole }))) : []
  // A PARENT'S OWN COPY, WHERE THE WINDOW SCAN HAS TOO FEW MARKERS: a short loss, or a parent
  // rebuilt from few cells. A run the window scan already reports stands as that scan gave it.
  const ownSegments: Segment[] = []
  for (const [markers, rates, role, windowed] of [
    [patMarkers, patRates, soloRole, windowPat], [matMarkers, matRates, otherRole, windowMat],
  ] as const) {
    const own = rates && ownRates([...markers.values()], rates.a)
    if (!own) continue
    const before = ownSegments.length
    for (const [c, ms] of markers) {
      for (const run of scanOwnCopy(c, ms, own, MIN_SEGMENT_SPAN_BP)) {
        const lrrShift = lrrShiftOf(c, run.startBp, run.endBp)
        if (underGain(run, lrrShift)) continue
        if (windowed.some((w) => w.chrom === c && segmentCoords(w).start <= run.endBp
          && segmentCoords(w).end >= run.startBp)) continue
        ownSegments.push({ ...run, absentParent: role, lrrShift })
      }
    }
    log('SCAN', `${role} own-copy scan: ${markers.size} chromosomes carrying both parents, clustering `
      + `${own.inflation.toFixed(1)}, ${ownSegments.length - before} run(s) the window scan does not give`)
  }
  const scanFound = new Set([...otherSegments, ...ownSegments])
  const allSegments = [
    ...copy,
    ...scanFound,
    ...windowPat,
    ...windowMat,
  ].sort((a, b) => b.score - a.score)
  const segmentsWithheld = uninterpretable && !chaoticGenome
  if (segmentsWithheld) {
    // `copy` is already empty here and its segments already recorded, so these are the absence
    // segments alone and nothing is counted twice.
    for (const sg of allSegments) {
      result.refusals!.push({ claim: 'segment', chrom: sg.chrom, gate: stageGate, cause: stageCause,
        measured: stageMeasured, withheld: sg })
    }
  }
  result.segments = segmentsWithheld ? [] : allSegments
  // A WINDOW CALLED BY BOTH RUNS NAMES NO PARENT. The two runs read different marker sets, so
  // the same region can clear the threshold in both, and that is not an event: one lost copy
  // takes one parent's alleles with it. Both depleted over one window is either a region missing
  // both copies, which the intensity channel calls as a copy loss on its own, or dropout wide
  // enough to clear the window minimum. Printing the pair would show a reader two contradictory
  // findings on one region; dropping it would lose a real position. So the region is kept and the
  // parent is withheld.
  //
  // THE TEST IS RECIPROCAL OVERLAP, NOT INTERSECTION. A bare intersection test was what shipped
  // first, and at the near-chromosome window lengths this corpus produces, any two events on one
  // chromosome intersect: one arm losing the father's copy and the other arm losing the mother's
  // is two findings, and intersection collapses them into one contested region with no parent.
  // The edges are `segmentCoords`, refined where refinement localised the event, which is this
  // module's own answer to what a segment's coordinates are.
  for (const a of result.segments) {
    if (a.kind !== 'parental-absence' || !a.absentParent) continue
    const ca = segmentCoords(a)
    for (const b of result.segments) {
      if (b === a || b.kind !== 'parental-absence' || !b.absentParent) continue
      if (b.chrom !== a.chrom || b.absentParent === a.absentParent) continue
      const cb = segmentCoords(b)
      const inter = Math.min(ca.end, cb.end) - Math.max(ca.start, cb.start)
      if (inter <= 0) continue
      if (inter / (Math.max(ca.end, cb.end) - Math.min(ca.start, cb.start)) < RECIPROCAL_OVERLAP) {
        continue
      }
      a.contested = true
      break
    }
  }
  if (chaoticGenome) {
    result.limits.push(
      'This array reads above the heterozygosity a single genome can reach, and the departure '
      + `from a uniform genome is LOCALISED: the spread of heterozygosity across the chromosomes `
      + `that are being genotyped is ${spread.toFixed(2)} against a uniform reading of under `
      + `${HET_SPREAD_LOCALISED.toFixed(2)}. A mixture of two people or a contaminated reaction `
      + 'raises the rate at every marker on every chromosome, because two genomes differ '
      + 'everywhere. A genome with chromosomes at one copy has no heterozygotes to call on those '
      + 'and a normal rate on the rest. Chromosomes whose call rate has collapsed are left out of '
      + 'that measurement, because a region that failed to amplify loses its heterozygotes too '
      + 'and is not evidence of a copy number. So the regions below are '
      + 'reported from the ALLELES, which need no baseline. Copy number is withheld throughout: '
      + 'the intensity is normalised against this sample\'s own per-chromosome median, and a '
      + 'genome without a normal majority has no such median. The whole-chromosome verdicts are '
      + 'withheld as well, on gates that have not been measured on material like this.',
    )
  }
  // WHERE THE SCAN HAS GONE BLIND, SAY SO. The null each chromosome is scored against is the
  // median absence rate over this sample's OTHER autosomes, which holds only while most of them
  // are unaffected. Past that the bar rises toward the rate a real event produces and no window
  // can clear the threshold, which is indistinguishable in the output from a chromosome that
  // carries nothing. The rate a true single-parent loss produces is computed from the two
  // parental arrays, so the test references the parents and not this sample's other chromosomes.
  if (mat) {
    // Whichever parent's loss is harder to see sets the bound, so the smaller of the two rates.
    const lossRate = Math.min(
      lossSignalRate(pat.gt, mat.gt),
      lossSignalRate(mat.gt, pat.gt),
    )
    const blind: string[] = []
    for (const c of result.chroms) {
      if (!isAutosome(c.chrom) || !measured.has(c.chrom)) continue
      const nPat = externalNull(t.byChrom, c.chrom)
      const nMat = externalNull(matAbsenceTally, c.chrom)
      const worst = Math.max(nPat, nMat)
      if (!canDetectLoss(worst, lossRate)) blind.push(c.chrom)
    }
    if (blind.length) {
      result.limits.push(
        `The segment scan cannot detect a lost parental copy on ${blind.length} of the autosomes `
        + `scanned here (chr${blind.join(', chr')}), so silence on them is not evidence of their `
        + 'being intact. Each chromosome is scored against the median absence rate of this '
        + "sample's OTHER chromosomes, and on a genome where most of them carry events that "
        + `median is itself an affected chromosome. A true lost copy reads ${pct(lossRate, 1)} `
        + 'here, and the bar on those chromosomes has risen far enough that no window of the '
        + 'minimum size could clear the calling threshold even if the whole window were missing.',
      )
    }
  }
  const contested = result.segments.filter((sg) => sg.contested)
  if (contested.length) {
    const where = [...new Set(contested.map((sg) => `chr${sg.chrom}`))].join(', ')
    result.limits.push(
      `On ${where} the absence scan called an overlapping region against BOTH parents. A single `
      + 'lost copy takes one parent\'s alleles with it, so both cannot be the missing ones over '
      + 'one window: either the region is missing both copies, which the intensity channel calls '
      + 'as a copy loss on its own evidence, or dropout there is wide enough to clear the window '
      + 'minimum. The regions are reported with no parent named.',
    )
  }

  // WHICH PARENT A SEGMENT'S CHANGE FELL ON. The same statistic as the whole-chromosome scan,
  // computed over the informative markers that fall INSIDE the called window instead of across
  // the chromosome, and judged against the same cell-own null of biparental whole chromosomes so
  // that a segment and a chromosome are measured on one scale. A whole chromosome dilutes an
  // arm-sized event across its full length, which is why an arm reads near zero there and can
  // read plainly here.
  if (chromShares.length) {
    const segShares: ChromShare[] = result.segments.map((sg) => {
      const inside = (dosageByChrom.get(sg.chrom) ?? [])
        .filter((m) => m.pos >= sg.startBp && m.pos <= sg.endBp)
      const lrrs = (cnByChrom.get(sg.chrom) ?? [])
        .filter((m) => m.pos >= sg.startBp && m.pos <= sg.endBp && m.log2R !== null)
        .map((m) => m.log2R as number)
      return {
        chrom: sg.chrom,
        label: `${sg.chrom}:${(sg.startBp / 1e6).toFixed(1)}-${(sg.endBp / 1e6).toFixed(1)}Mb`,
        share: inside.length ? median(inside.map((m) => m.patShare)) : NaN,
        informative: inside.length,
        biparental: true,
        lrrShift: lrrs.length ? median(lrrs) - genomeLrr : NaN,
      }
    })
    const byLabel = new Map(scanAlleleGains(segShares, chromShares).map((g) => [g.chrom, g]))
    for (const sg of result.segments) {
      const g = byLabel.get(`${sg.chrom}:${(sg.startBp / 1e6).toFixed(1)}-${(sg.endBp / 1e6).toFixed(1)}Mb`)
      if (g && (g.verdict === 'gain' || g.verdict === 'imbalance')) sg.alleleGain = g
    }
  }

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
    // meaningful where the array is reading a genome at all. Detecting structure in an array
    // that failed its own quality inference is not a conservative reading of weak data, it is
    // reading noise as biology. The event list stays empty, on the same condition and for the
    // same reason as the aneuploidy and segment scans above it.
    //
    // A DROP-IN ARRAY IS READING A GENOME, and its detectors run, each judged on what drop-in does to
    // it. Drop-in only ever adds heterozygous calls, so it can break a run of homozygosity and never
    // make one: a run or a depleted stretch found through it is there despite it. What drop-in does
    // corrupt is the B allele frequency, which triploidy reads, so that one is withheld under its
    // physical cause; and the complex-genome count reads copy number, so it goes where copy number goes.
    if (uninterpretable) {
      result.refusals!.push({ claim: 'findings', gate: stageGate, cause: stageCause, measured: stageMeasured })
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
      // AND THE SAME RULE THE COPY-NUMBER SCAN ABOVE APPLIES, for the same reason: a
      // chromosome already called whole is not scanned for a run inside it, since that is
      // the same event described twice. A chromosome at ONE copy is hemizygous end to end,
      // so it reads homozygous end to end and produces a full-length run; classified here
      // that run becomes an isodisomy, which asserts BOTH copies came from one parent on an
      // event that has one copy. The loss is the true row and it is reported through
      // `result.chroms`; a real isodisomy has two copies, leaves intensity at normal and is
      // never in `lost`, so it is untouched by this.
      //
      // `lost`, NOT `whole`. Filtering every aneuploid chromosome took the gains with it, and a
      // gain that is homozygous end to end is an isodisomic trisomy: a real and separately
      // meaningful finding that the gain row does not report.
      ...detectUpd(
        runsOfHomozygosity(selfMarkers, { chromEndBp: chromEnd })
          .filter((r) => !lost.has(r.chrom)),
        {
          zygosity: result.zygosity,
          backgroundHet: selfMarkers.length
            ? selfMarkers.filter((m) => m.het).length / selfMarkers.length : undefined,
          intensity: lrrAll.length > 0,
        },
      ),
    ]
    log('SCAN', `runs of homozygosity and ploidy over ${selfMarkers.length} called markers`)
    await breathe()
    if (dropInCell) {
      result.refusals!.push({ claim: 'triploidy', gate: 'stage.drop-in', cause: 'baf-displaced',
        measured: dropInMeasured })
    } else {
      const tri = detectTriploidy([...myBaf.values()])
      if (tri) findings.push(tri)
    }
    // A genome with too little undisturbed remainder cannot self-reference, which blocks
    // every origin call on the array rather than only on the affected chromosomes. It counts
    // copy-number events, so where copy number is withheld it is withheld with it, on record.
    if (copyWithheld) {
      result.refusals!.push(intensityWithheld
        ? { claim: 'complex', gate: 'normal-majority', cause: 'no-intensity-reference', measured: dropInMeasured }
        : { claim: 'complex', gate: 'copy-scan.chaotic-genome', cause: 'proxy',
          measured: { ...stageMeasured, heterozygositySpread: spread } })
    } else {
      const deviant = new Set([...whole, ...result.segments.map((sg) => sg.chrom)]).size
      const cx = detectComplex(deviant, cnByChrom.size)
      if (cx) findings.push(cx)
    }
    // A COPY-NEUTRAL CLAIM SAYS THERE ARE TWO COPIES, and where copy number has no reference the count
    // is unknown. A chromosome down to one copy reads homozygous end to end, which is a run; with
    // copy number read it is filtered as a loss above, and without it the run is classed as a
    // disomy of one parent. On GSE148488 that put a UPD or isodisomy on 10 chromosomes of 4 cells
    // whose copy number is withheld, every one of them published monosomic or deleted. So these
    // go where copy number goes, on record, under the same cause.
    const COPY_NEUTRAL = new Set(['cnn-loh', 'isodisomy', 'heterodisomy', 'segmental-upd'])
    const kept = copyWithheld ? findings.filter((f) => !COPY_NEUTRAL.has(f.cls)) : findings
    for (const f of findings) {
      if (kept.includes(f)) continue
      result.refusals!.push(intensityWithheld
        ? { claim: 'findings', chrom: f.chrom, gate: 'normal-majority', cause: 'no-intensity-reference',
          measured: dropInMeasured, withheld: f }
        : { claim: 'findings', chrom: f.chrom, gate: 'copy-scan.chaotic-genome', cause: 'proxy',
          measured: { ...stageMeasured, heterozygositySpread: spread }, withheld: f })
    }
    log('SCAN', `taxonomy: ${kept.length} finding`
      + `${kept.length === 1 ? '' : 's'} across `
      + `${new Set(kept.map((f) => f.chrom)).size} chromosome(s)`)
    await breathe()
    result.findings = kept
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
    const h = hetByChrom.get(chrom) ?? { offered: 0, informative: 0, het: 0 }
    const background = externalHetBackground(hetByChrom, chrom)
    const c = callHomologue(
      h.informative ? h.het / h.informative : NaN, h.informative, background,
    )
    return {
      where,
      kind,
      origin: c.verdict === 'other homologue'
        ? 'this parent, other homologue (meiotic)' : 'not determinable',
      // The extra-copy sentence is about a gain; a loss here carries the homologue reading alone.
      why: mat || event === 'loss' ? c.why
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
  // NOT GATED ON THE PARENT COUNT. A second parent is strictly more evidence: it lowers every
  // dosage floor it touches (an esc-single whole-chromosome loss goes from 0.348 to 0.232) and
  // takes nothing away. Scoping this block to one-parent runs meant loading the oocyte donor
  // emitted no dosage call and no Mendelian call at all, so the better-equipped run was the
  // silent one, and `parents: mat ? 2 : 1` below could never be anything but 1.
  // AND SEGMENTS COUNT AS SOMETHING TO SCORE. They were reaching the Mendelian channel and
  // nothing else, so on a sample where that channel is off, a detected segment left the run with
  // no origin output of any kind: measured on a paternal pronucleus carrying one constructed
  // 40 Mb loss, 0 dosage rows and 0 Mendelian rows for an event the scan had just reported.
  if (whole.size || (result.findings?.length ?? 0) > 0 || result.segments.length) {
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
    // THE ARRAY'S OWN NULL, built once. A unit is a whole chromosome: its median log2R. The scale
    // is a median absolute deviation rather than a standard deviation because a plain SD over
    // chromosomes is inflated 2.34x by a SINGLE affected chromosome, so one real trisomy destroys
    // the null meant to detect it. The MAD holds to three or four.
    const chromNull = nullScale(
      [...cnByChrom].filter(([c]) => isAutosome(c)).map(([, ms]) => {
        const v = ms.map((m) => m.log2R).filter((x): x is number => x !== null)
        return median(v)
      }).filter((x) => Number.isFinite(x)),
    )
    log('SCAN', `intensity null from ${chromNull.units} autosomes, `
      + `centre ${chromNull.centre.toFixed(4)}, scale ${chromNull.scale.toFixed(4)}`)

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
      const { region, background, inL, untRows } =
        gatherInterval(scanIndex, iv, { regionOnly: unreachable })
      const mean = (xs: number[]) => (xs.length
        ? xs.reduce((a, x) => a + x, 0) / xs.length : NaN)
      // `outL` and a difference of means are no longer used for the error: the region is
      // compared to the array's own units instead. Keeping the old local would invite it back.
      // MEASURED AGAINST THE ARRAY'S OWN UNITS, NOT sd/sqrt(n). The iid form assumes markers are
      // independent; on amplified material they are long-range dependent, so it understates the
      // error by 3.6x at 50 markers rising to 16.5x at 20,000, and it flagged 89.7 percent of
      // event-free chromosomes. 55 percent of those false positives read as gains, and a gain
      // inverts the sign map that loss and copy-neutral share, so the parent came out backwards.
      // See intensityNull.ts.
      const intensityZ = calibratedZ(mean(inL), chromNull)

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
        // Or too few calls to form the reference, which is the array's property and not a finding.
        noSelfReference: !!result.findings?.some((f) => f.cls === 'complex')
          || genomeCallRate < COMPLEX_CALL_RATE,
        // INTENSITY DETECTS OR IT SAYS NOTHING. Below the calibrated threshold the SIGN of a shift
        // is close to a coin flip: on 264 event-free chromosomes 59 percent of the old statistic's
        // false positives read POSITIVE, a positive shift resolves the class as a gain, and a gain
        // inverts the sign map that loss and copy-neutral share, so the parent came out backwards.
        // Passing a sub-threshold z contributes that coin flip to the class. Passing nothing
        // contributes nothing, which is the honest value of the evidence.
        intensityZ: intensityDetects(intensityZ, wholeChromosome) ? intensityZ : undefined,
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
          // THE FIGURE IS HISTORY AND HAS TO READ AS HISTORY. Since 5.19.0 no band-F row from
          // either channel carries a parent, so "this band recovers the parent 0.51 to 0.56 of
          // the time" has no referent in the shipped tool: there is no such parent to be right
          // or wrong about. It is quoted as the evidence that RETIRED that caller, which is
          // what makes the refusal below something other than caution.
          why: `${c.why}. Graded F: this interval carried no usable marker of its own, `
            + 'so nothing here bears on which parent it was. The parent is NOT named. That '
            + 'withhold is not caution: when an earlier build did name a parent in this band, '
            + 'an injection series on real arrays scored it right 0.51 to 0.56 of the time, '
            + 'which is chance, because every call reaching this band has an unresolved '
            + 'copy-number class and a gain inverts the sign that loss and copy-neutral share',
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
            + 'information at all. When an earlier build did name a parent in this band, an '
            + 'injection series on real arrays scored it right 0.51 to 0.56 of the time, '
            + 'which is chance. No parent is named here',
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
      result.dosageCalls.push(withheldIfSlotWrong(scoreInterval(`chr${chrom}`, { chrom }, true)))
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
      result.dosageCalls.push(withheldIfSlotWrong(call))
      scoredSoFar += 1
      // One yield every couple of findings, so the lines above appear as they are decided
      // rather than all at once when the loop ends.
      if (scoredSoFar % 2 === 0) await breathe()
    }

    // THE SEGMENTS GO THROUGH IT TOO, and for the same reason the findings do: a reader comparing
    // a segment row against a whole-chromosome row must be comparing one posterior and one set of
    // bands. Where this material cannot answer, the scorer says which channel refused and why, and
    // an F is not a result; silence was the one outcome that told the reader nothing.
    //
    // A segment on a chromosome already called whole is skipped, the same rule the copy-number
    // scan applies: that is the same event described twice, and the whole-chromosome row above
    // already carries it.
    const segs = result.segments.filter((sg) => !whole.has(sg.chrom))
    if (segs.length) {
      log('SCAN', `scoring parental origin on ${segs.length} segment`
        + `${segs.length === 1 ? '' : 's'}`)
    }
    for (const sg of segs) {
      const co = segmentCoords(sg)
      result.dosageCalls.push(withheldIfSlotWrong(scoreInterval(
        `chr${sg.chrom} ${(co.start / 1e6).toFixed(1)}-${(co.end / 1e6).toFixed(1)}Mb`,
        { chrom: sg.chrom, startBp: co.start, endBp: co.end },
        false,
        sg.kind === 'copy-gain' ? 'gain' : 'loss',
        sg.kind,
      )))
      await breathe()
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
  // Nor on a genome whose ploidy was refused: `unknown` may be the one-parent case this guard
  // exists for. See knownBiparental.
  // AND NOT ON A GENOME THAT CARRIES ONE PARENT THROUGHOUT, measured the same way as the allele
  // count: at least two autosomes must read both parents at the markers where the loaded parent
  // is homozygous. The zygosity gate alone passes a diploid androgenetic genome, which a complete
  // mole is, and a parent's own array in the sample slot; on either, every chromosome reads the
  // loaded parent's copy present and the other's absent, which is true and says nothing about any
  // interval. Applied only where the reference is bulk, because against an amplified reference the
  // per-chromosome reading does not separate (audit section E2) and the gate would refuse at
  // random.
  const canAttribute = knownBiparental(result.zygosity) && (!bulkReference || bothParentsSomewhere)
  // Also not gated on the parent count: this channel reads whether the LOADED parent's allele is
  // present, which a second parent neither supplies nor obstructs.
  // AND THE MENDELIAN CHANNEL DOES NOT RUN AT ALL WITH A WRONG PARENTAL SLOT. Unlike dosage, every
  // row it produces IS a parental attribution: there is no location-only part of it to keep. A row
  // reading "the loaded parent's copy is absent" is only meaningful once "the loaded parent" is
  // the person the report says it is.
  const mendelEvents = canAttribute && !parentSlotWrong ? [
    ...(result.chroms ?? []).filter((c) => c.aneuploidy || c.lostParent)
      .map((c) => ({ chrom: c.chrom, start: 0, end: Number.MAX_SAFE_INTEGER,
        label: `chr${c.chrom}` })),
    // A segment a likelihood scan found is named by that scan and is not asked again here:
    // nameAbsentParent counts opposite homozygotes at every marker, and on a stretch at one copy
    // the miscalls among them pass its bound.
    ...result.segments.filter((sg) => !scanFound.has(sg)).map((sg) => {
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
    //
    // AND WITH BOTH ARRAYS LOADED, EACH IS ASKED ONLY THE QUESTION IT CAN ANSWER. This channel
    // is right 0.8539 of the time about the loaded parent's OWN copy and 0.2890 about the other
    // parent's, on the same arrays and the same events; the difference is structural, since the
    // second reading rests on missing heterozygosity and missing heterozygosity is what dropout
    // produces. Loading only the first parent therefore answered about half of all losses by the
    // unreliable direction. See bothParentsOrigin.ts.
    result.oneParent = mendelEvents.map((sg) => {
      const co = { start: sg.start, end: sg.end }
      const pairs: [string, string][] = []
      const otherPairs: [string, string][] = []
      for (const [probe, gt] of myGt) {
        const q = markerPos.get(probe)
        if (!q || q.chrom !== sg.chrom || q.pos < co.start || q.pos > co.end) continue
        const pg = pat.gt.get(probe)
        if (pg) pairs.push([pg, gt])
        if (mat) {
          const mg = mat.gt.get(probe)
          if (mg) otherPairs.push([mg, gt])
        }
      }
      // Dropout comes from the inferred stage. A failed array carries no usable figure
      // rather than a flattering one, so it gets the most conservative dropout measured
      // on any stage instead of NaN, which would silently void every likelihood below.
      const ado = Number.isFinite(result.stage!.dropout) ? result.stage!.dropout : 0.308
      const c = callOneParentOrigin(pairs as never, ado)
      if (!mat) {
        // ARM (b). With one array the parent used to be left unnamed, because the only statistic
        // consulted was heterozygosity that is missing, and missing heterozygosity is what dropout
        // produces. Two statistics computed INSIDE this region settle it without the second parent:
        // whether one contribution is present or two, from markers where the loaded parent is
        // heterozygous, and whether the copy that remains is the loaded parent's, from the
        // opposite-homozygote rate at markers where that parent is homozygous. See nameAbsentParent.
        let het = 0; let markers = 0; let opp = 0; let oppMarkers = 0
        for (const [pg, gt] of pairs) {
          if (pg === 'AB') { markers += 1; if (gt === 'AB') het += 1 } else if (gt !== 'AB') {
            oppMarkers += 1; if (gt !== pg) opp += 1
          }
        }
        const absent = nameAbsentParent(het, markers, opp, oppMarkers, ado)
        const named = absent === 'loaded' ? soloRole : absent === 'other' ? otherRole : undefined
        return {
          where: sg.label,
          verdict: c.verdict,
          posterior: c.posterior,
          band: c.band,
          markers: c.markers,
          exclusive: c.exclusive,
          why: named
            ? `${c.why} With one array loaded the absent copy is named from the allele classes `
              + `this region carries: ${markers.toLocaleString()} markers where the `
              + `${soloRole} array is heterozygous read ${pct(het / Math.max(1, markers), 1)} `
              + `heterozygous, which is one contribution, and `
              + `${pct(opp / Math.max(1, oppMarkers), 2)} of ${oppMarkers.toLocaleString()} `
              + `markers where it is homozygous carry the opposite homozygote, which is `
              + `${absent === 'loaded' ? 'too high for its own copy to be present'
                : 'too low for its own copy to be absent'}. So the ${named} copy is the absent one.`
            : c.why,
          parent: named,
          oneArrayNamed: named ? true : undefined,
        }
      }
      const o = callOneParentOrigin(otherPairs as never, ado)
      // The rule takes the paternal call first, so the two are ordered by role rather than by
      // which array happened to be loaded into which slot.
      const bothCall = soloRole === 'paternal'
        ? callBothParentsOrigin(c, o) : callBothParentsOrigin(o, c)
      // Markers and exclusives stay those of the side the answer came from, so the evidence a
      // reader checks by hand is the evidence the verdict rests on.
      const side = bothCall.parent === soloRole ? c : o
      return {
        where: sg.label,
        verdict: bothCall.verdict,
        posterior: bothCall.posterior,
        band: bothCall.band,
        markers: side.markers,
        exclusive: side.exclusive,
        why: bothCall.why,
        twoParents: true,
        parent: bothCall.parent,
        corroborated: bothCall.corroborated,
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

  // THE SAMPLE-LEVEL FLAG, computed last because it summarises everything above it. Reusing the
  // bounds the taxonomy already measured rather than choosing new ones: see integrity.ts.
  {
    const autosomes = (result.chroms ?? []).filter((c) => isAutosome(c.chrom))
    const affected = new Set<string>()
    for (const c of autosomes) if (c.aneuploidy) affected.add(c.chrom)
    for (const sg of (result.segments ?? [])) affected.add(sg.chrom)
    for (const f of (result.findings ?? [])) if (f.chrom !== 'genome') affected.add(f.chrom)
    result.integrity = assessIntegrity({
      stage: result.stage?.stage ?? 'unknown',
      callRate: genomeCallRate,
      affectedAutosomes: affected.size,
      totalAutosomes: autosomes.length,
      events: (result.segments?.length ?? 0) + (result.findings?.length ?? 0)
        + autosomes.filter((c) => c.aneuploidy).length,
      parentAbsent: autosomes.filter((c) => c.lostParent && !c.aneuploidy).length,
      gainByAlleles: autosomes.filter((c) => c.alleleGain?.verdict === 'gain' && c.aneuploidy !== 'gain').length,
      allelesStand: (uninterpretable && alleleClaims.size > 0) || intensityWithheld,
      normalFraction: normal,
    })
    // A LABEL IS A REFUSAL TOO. These levels remove nothing, but tell the reader nothing on the
    // sample is reportable, which is the same to the person reading it. Recorded as soft.
    if (result.integrity.level !== 'intact') {
      const m = result.integrity.measures
      result.refusals!.push({ claim: 'sample', gate: `integrity.${result.integrity.level}`,
        cause: 'proxy', soft: true,
        measured: { callRate: m.callRate, deviantFraction: m.deviantFraction,
          affectedAutosomes: m.affectedAutosomes, totalAutosomes: m.totalAutosomes } })
    }
    log(result.integrity.level === 'uninterpretable' ? 'WARN' : 'DONE',
      `${sampleName}: ${result.integrity.headline}`)
  }
  // After the integrity call, so a sample not reportable logs what was seen, not a change.
  for (const [level, line] of changeLogLines(result, sampleName, soloRole)) log(level, line)

  return result
}
