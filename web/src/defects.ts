/**
 * What the defect display is told, separate from how it is drawn.
 *
 * Kept out of the .tsx for two reasons: node cannot strip JSX, so a component file cannot carry a
 * self-check, and this is the layer where getting it wrong matters. Showing a confident parent for
 * a region that never had one is the most damaging thing this UI could do, so the mapping that
 * decides it is testable on its own.
 */
import type { Segment } from './segments.ts'
import type { GainAnnotation } from './parentage.ts'
import { segmentCoords } from './segments.ts'
import { locus } from './stage.ts'
import type { StageCall } from './stage.ts'
import { ORIGIN_UNREACHABLE, taxonomyFor, type Finding } from './abnormalities.ts'

/** What is known about one defect, gathered from whichever channels produced it. */
export interface Defect {
  chrom: string
  startBp: number
  endBp: number
  /**
   * What the change IS.
   *
   * The first four are the original three plus the whole-chromosome case; everything after them
   * comes from the taxonomy, and they are the same field on purpose. A reader compares rows, so a
   * copy-neutral event and a deletion must arrive in one list with one shape rather than in two
   * displays that cannot be read against each other.
   */
  kind: 'copy-loss' | 'copy-gain' | 'parental-absence' | 'whole-chromosome'
    | 'cnn-loh' | 'isodisomy' | 'segmental-upd' | 'triploidy' | 'haploidy' | 'complex'
    | 'monosomy' | 'trisomy' | 'segmental-deletion' | 'segmental-duplication' | 'gamete-de-novo'
  /** The region in genome-browser form, chr6:39,302-294,904. */
  locus: string
  /** Refined interval text where the breakpoint was localised, e.g. '+/- 61 kb'. */
  interval?: string
  /** Origin as the tool determined it, with the reason. */
  origin: 'paternal' | 'maternal' | 'unclear'
  why: string
  /** How the origin was reached, which decides how much weight it carries. */
  basis?: 'two-parent' | 'one-parent' | 'sibling' | 'homologue' | 'dosage'
  /** Which channel produced it. Genotypes where they exist, dosage where they have collapsed. */
  channel?: 'genotype' | 'dosage'
  /** Self-referenced centroid shift and its z, the dosage channel's own evidence. */
  shift?: number
  z?: number
  /** Mosaic fraction the shift implies, under the class the posterior found most probable. */
  impliedF?: number
  informative?: number
  posterior?: number
  /**
   * Calibrated probability that the named parent is the right one, and the band it lands in.
   *
   * A probability rather than a likelihood ratio or a bare label, because a reader compares rows:
   * an LR of 12 means different things at 100 markers and at 800, a confidence of 0.92 does not.
   * Every band carries its number including the weakest, since each is calibrated within itself.
   * What changes across bands is the words beside the number.
   */
  confidence?: number
  /** 'A' very confident, 'B' confident, 'C' weak direction only, 'D' weak not for reporting. */
  band?: string
  /**
   * Why the confidence is what it is, because the causes have different remedies: a low fraction
   * wants more cells of the same embryo, poor amplification wants a re-amplification, few markers
   * cannot be helped at a fixed platform, and an unresolved class wants the second parent.
   */
  limitedBy?: string
  /**
   * Whether this change came from the gamete or arose after fertilisation.
   *
   * Present only when two or more units of the same embryo were arrayed. Genotype alone cannot
   * separate the two at any material quality, so on a single-unit run this is absent and the
   * report says what would supply it rather than leaving a blank.
   */
  mechanism?: 'meiotic' | 'post-zygotic' | 'unresolved'
  mechanismWhy?: string
  /** Markers carrying an allele the loaded parent does not have. This is the Mendelian evidence
   *  itself rather than a summary of it: a parent who is AA has no B to give, and dropout removes
   *  alleles without inventing one, so a non-zero count cannot come from amplification loss. */
  exclusive?: number
  /** The material this was read from, carried so a chip can state it without a second lookup. */
  stage?: string
  /** Dropout the call was parameterised with, and how that figure was arrived at. */
  dropout?: number
  dropoutBasis?: string
}

/** Build the display list from a run's segments and whatever origin annotations exist for them. */
export function defectsFrom(
  segments: readonly Segment[],
  gains: readonly GainAnnotation[],
  losses: readonly GainAnnotation[],
  /** Single-parent calls, keyed by the same `where` string. The loaded parent's role decides
   *  whether 'known-parent-lost' reads paternal or maternal. */
  oneParent: readonly {
    where: string, verdict: string, posterior: number, markers: number, exclusive: number,
    why: string,
    /** Capped at B: this channel's failure mode is dropout, which is the same event as the
     *  observation, so no amount of evidence earns the top band. */
    band?: string,
  }[] = [],
  loadedParent: 'paternal' | 'maternal' = 'paternal',
  stage?: StageCall,
  /** Whole-chromosome calls from allele dosage, keyed `chrN`. Fill a gap the genotype channel
   *  cannot fill rather than competing with it: they never cover the same event. */
  dosageCalls: readonly {
    where: string, verdict: string, shift: number, z: number, impliedF: number,
    window: number, why: string,
    /** Calibrated confidence in the named parent, and the band it lands in. */
    confidence?: number, band?: string, limitedBy?: string,
  }[] = [],
): Defect[] {
  const byOne = new Map(oneParent.map((o) => [o.where, o]))
  const byDosage = new Map(dosageCalls.map((d) => [d.where, d]))
  const other = loadedParent === 'paternal' ? 'maternal' : 'paternal'
  const byWhere = new Map<string, GainAnnotation>()
  for (const a of [...gains, ...losses]) byWhere.set(a.where, a)
  return segments.map((sg) => {
    const co = segmentCoords(sg)
    const where = `chr${sg.chrom} ${(co.start / 1e6).toFixed(1)}-${(co.end / 1e6).toFixed(1)}Mb`
    const ann = byWhere.get(where)
    let origin: Defect['origin'] = ann?.origin?.startsWith('paternal') ? 'paternal'
      : ann?.origin?.startsWith('maternal') ? 'maternal' : 'unclear'
    // A single-parent call names a parent where two-parent dosage could not, so it fills the gap
    // rather than overriding: two parents remain the stronger evidence where both were loaded.
    const one = byOne.get(where)
    if (origin === 'unclear' && one) {
      if (one.verdict === 'known-parent-lost') origin = loadedParent
      else if (one.verdict === 'other-parent-lost') origin = other
    }
    // Dosage fills what genotypes could not reach. A whole chromosome is detected by its genotype
    // call rate collapsing, so on those events the genotype channel has no evidence left and this
    // is the only channel that can speak. It never overrides a genotype answer.
    const dose = byDosage.get(`chr${sg.chrom}`)
    // 'loaded-parent' and 'other-parent' replace the old 'known-parent-lost' / 'other-parent-lost'
    // on THIS channel. The old names presumed a loss, and a loss is exactly what the dosage channel
    // cannot presume: a gain inverts the sign, so the same shift names the opposite parent. The
    // verdict now comes from a posterior that marginalises the class rather than assuming it.
    const usedDosage = origin === 'unclear' && !!dose
      && (dose.verdict === 'loaded-parent' || dose.verdict === 'other-parent')
    if (usedDosage && dose) {
      origin = dose.verdict === 'loaded-parent' ? loadedParent : other
    }
    return {
      chrom: sg.chrom,
      startBp: co.start,
      endBp: co.end,
      locus: locus(sg.chrom, co.start, co.end),
      kind: sg.kind === 'copy-gain' ? 'copy-gain'
        : sg.kind === 'copy-loss' ? 'copy-loss' : 'parental-absence',
      interval: sg.refined ? co.interval : undefined,
      origin,
      why: ann?.why ?? (usedDosage ? dose?.why : undefined) ?? one?.why
        ?? 'No parental genotype was loaded for this run, so allele dosage cannot be oriented and '
          + 'the parent is not determined. The position and its interval are unaffected.',
      basis: ann ? 'two-parent'
        : usedDosage ? 'dosage'
          : (one && origin !== 'unclear') ? 'one-parent' : undefined,
      channel: usedDosage ? 'dosage' : (one ? 'genotype' : undefined),
      // EVERY CHANNEL SUPPLIES ITS OWN CONFIDENCE, taken from whichever one actually produced the
      // origin above. Before this, only dosage did, so a two-parent genotype call, the strongest
      // evidence here, printed a bare parent name while a weak dosage call printed 0.62 beside it.
      // The formatting told a reader the opposite of the truth.
      confidence: usedDosage ? dose?.confidence
        : (ann && origin !== 'unclear') ? ann.confidence
          : (one && origin !== 'unclear') ? one.posterior : undefined,
      band: usedDosage ? dose?.band
        : (ann && origin !== 'unclear') ? ann.band
          : (one && origin !== 'unclear') ? one.band : undefined,
      limitedBy: usedDosage ? dose?.limitedBy : undefined,
      shift: dose && Number.isFinite(dose.shift) ? dose.shift : undefined,
      z: dose && Number.isFinite(dose.z) ? dose.z : undefined,
      impliedF: dose && Number.isFinite(dose.impliedF) ? dose.impliedF : undefined,
      informative: one?.markers,
      posterior: one && Number.isFinite(one.posterior) ? one.posterior : undefined,
      exclusive: one?.exclusive,
      stage: stage?.stage,
      dropout: stage && Number.isFinite(stage.dropout) ? stage.dropout : undefined,
      dropoutBasis: stage?.basis === 'none' ? undefined : stage?.basis,
    } as Defect
  })
}

// ---------------------------------------------------------------------------------------------
// Display helpers. Here rather than in the .tsx because node cannot strip JSX, and this is the
// layer where getting it wrong matters: a parent shown without its confidence reads as certain.

/**
 * The headline: whose copy, how sure, or that it could not be said.
 *
 * THE CONFIDENCE TRAVELS WITH THE PARENT AND IS NOT A SEPARATE CHIP. A parent named without a
 * number reads as certain, and on amplified material most of them are not: the measured accuracy
 * of a top-band call is 0.995 while the weakest band is 0.62. Those two must not look alike at a
 * glance, so the band label sits in the headline where the parent is, and the chips below carry
 * the evidence rather than the verdict.
 */
export const BAND_WORD: Record<string, string> = {
  A: 'very confident', B: 'confident', C: 'weak, direction only', D: 'weak, not for reporting',
}
/**
 * What each class is called in a headline, in words rather than in the code's vocabulary.
 *
 * A reader of a clinical report should not have to know that 'cnn-loh' means the chromosome is
 * still there in two copies but both came from one parent. Every class gets a phrase.
 */
export const KIND_WORD: Record<string, string> = {
  'copy-gain': 'extra copy',
  'copy-loss': 'copy lost',
  'parental-absence': "parent's alleles absent",
  'whole-chromosome': 'whole chromosome affected',
  monosomy: 'one copy only',
  trisomy: 'three copies',
  'segmental-deletion': 'segment deleted',
  'segmental-duplication': 'segment duplicated',
  'cnn-loh': 'both copies from one parent, copy number normal',
  isodisomy: 'both copies are the same parental homologue',
  'segmental-upd': 'segment from one parent only',
  triploidy: 'a third chromosome set',
  haploidy: 'one chromosome set only',
  complex: 'genome too disturbed to reference against',
  'gamete-de-novo': 'segment changed, gamete or post-fertilisation',
}

export const headline = (d: Defect): string => {
  const what = KIND_WORD[d.kind] ?? (d.kind === 'copy-gain' ? 'extra copy' : 'copy lost')
  const where = d.chrom === 'genome' ? 'genome' : `chr${d.chrom}`
  if (d.origin === 'unclear') {
    // A class that can never carry a parent must not read like one that merely failed to.
    return `${where} ${what}, ${originBlockedByClass(d.kind)
      ? 'no parental origin exists for this class'
      : 'origin not determined'}`
  }
  const conf = d.confidence !== undefined && Number.isFinite(d.confidence)
    ? ` · ${d.confidence.toFixed(3)}${d.band ? ` ${BAND_WORD[d.band] ?? d.band}` : ''}`
    : ''
  // A gamete-borne change and one that arose after fertilisation are different findings for a
  // patient, so the timing sits in the headline beside the parent rather than in a chip below.
  const when = d.mechanism === 'meiotic' ? ' · from the gamete'
    : d.mechanism === 'post-zygotic' ? ' · after fertilisation' : ''
  return `${d.origin.toUpperCase()} ${what} · ${where}${conf}${when}`
}

/** Bands C and D are dimmed, so a weak number cannot be mistaken for a strong one at a glance. */
export const bandColour = (d: Defect): string => {
  if (d.origin === 'unclear') return 'var(--om-text-dim)'
  return d.band === 'C' || d.band === 'D' ? 'var(--om-text-dim)' : 'var(--om-defect)'
}

// ---------------------------------------------------------------------------------------------
// The taxonomy's findings, entering the same list as everything else.

/**
 * Turn a taxonomy finding into a defect, so it reaches every display the older events already do.
 *
 * THE POINT IS THAT THERE IS NO SECOND LIST. A copy-neutral event, an isodisomy and a deletion are
 * different measurements of the same kind of thing, and a reader compares them against each other.
 * Giving the new classes their own panel would have been easier and would have reproduced the
 * defect this release exists to fix, where the strongest evidence carried the least visible
 * confidence because it happened to be displayed somewhere else.
 *
 * ORIGIN IS BLOCKED FOR SOME CLASSES BY THE CLASS ITSELF, not by weak evidence, and the two must
 * not look alike. A triploidy has no parental origin available at any quality because band
 * structure cannot say whose the extra set is; a complex genome has none because self-reference has
 * failed. Both arrive here with the reason attached rather than as an empty field.
 */
export function findingToDefect(
  f: Finding,
  stage?: StageCall,
  /** The dosage call scored over this finding's own interval, where the class allows one. */
  scored?: {
    verdict: string, shift: number, z: number, impliedF: number, why: string,
    confidence?: number, band?: string, limitedBy?: string,
  },
  loadedParent: 'paternal' | 'maternal' = 'paternal',
): Defect {
  const blocked = f.originBlocked ?? (ORIGIN_UNREACHABLE.has(f.cls)
    ? taxonomyFor(f.cls)?.origin : undefined)
  const other = loadedParent === 'paternal' ? 'maternal' : 'paternal'
  // A blocked class is never scored, so a scored call on one would be a caller error rather than
  // evidence, and is ignored here rather than trusted.
  const use = blocked ? undefined : scored
  const named = use && (use.verdict === 'loaded-parent' || use.verdict === 'other-parent')
  return {
    chrom: f.chrom,
    startBp: f.startBp,
    endBp: f.endBp,
    locus: f.chrom === 'genome' ? 'whole genome' : locus(f.chrom, f.startBp, f.endBp),
    kind: f.cls as Defect['kind'],
    origin: named ? (use.verdict === 'loaded-parent' ? loadedParent : other) : 'unclear',
    why: [f.evidence, blocked, use?.why, f.flag].filter(Boolean).join('. '),
    // Absent where nothing scored an origin, so a missing basis never implies a channel spoke and
    // declined. Where one did, this is the SAME dosage channel the older events use.
    basis: named ? 'dosage' : undefined,
    channel: use ? 'dosage' : undefined,
    confidence: named ? use.confidence : undefined,
    band: named ? use.band : undefined,
    limitedBy: named ? use.limitedBy : undefined,
    shift: use && Number.isFinite(use.shift) ? use.shift : undefined,
    z: use && Number.isFinite(use.z) ? use.z : undefined,
    impliedF: use && Number.isFinite(use.impliedF) ? use.impliedF : undefined,
    stage: stage?.stage,
  }
}

/** True where the class itself, rather than the evidence, is what stops a parent being named. */
export const originBlockedByClass = (kind: string): boolean =>
  ORIGIN_UNREACHABLE.has(kind as never)

/**
 * Attach per-event timing to defects, matched by overlap rather than by identical edges.
 *
 * MATCHED BY OVERLAP FOR THE SAME REASON THE UNIFORMITY TEST USES IT. Two biopsies of one embryo do
 * not place a breakpoint identically, and the defect list and the uniformity list are built from
 * different passes over the same events, so requiring identical coordinates would drop the match on
 * exactly the events that have one.
 *
 * Absent uniformity leaves every defect untouched, which is the single-unit case: the report then
 * says what a second array would supply rather than showing a blank field.
 */
export function withMechanism(
  defects: Defect[],
  uniformity?: readonly {
    chrom: string, startBp: number, endBp: number,
    mechanism: 'meiotic' | 'post-zygotic' | 'unresolved', why: string,
  }[],
): Defect[] {
  if (!uniformity?.length) return defects
  return defects.map((d) => {
    const hit = uniformity.find((u) => u.chrom === d.chrom
      && u.startBp <= d.endBp && d.startBp <= u.endBp)
    return hit ? { ...d, mechanism: hit.mechanism, mechanismWhy: hit.why } : d
  })
}
