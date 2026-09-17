/**
 * How many X chromosomes, and whether one of them is a male karyotype or a monosomy.
 *
 * Heterozygosity cannot separate those two. A 46,XY male and a 45,X female both carry one X and
 * neither is heterozygous outside the pseudoautosomal region, so a sex call taken from the
 * heterozygosity ratio alone reads every 45,X as male. Two further measurements do separate them,
 * and neither is optional:
 *
 *   X copy number   chrX's own median intensity against the sample's autosomes. Measured on the
 *                   CytoScan miscarriage series, which states a karyotype per sample: one X sits
 *                   at -0.49 to -0.63 log2 (20 stated 45,X and 27 normal males, 10th to 90th
 *                   centile), two X at -0.05 to +0.01 (33 normal females). Nothing measured in
 *                   between, which is why the autosomal magnitude floor separates them.
 *   chrY            present or absent. THIS is what tells the two one-X constitutions apart, and
 *                   a panel that cannot measure chrY at all cannot tell them apart at all.
 *
 * A panel may carry chrY as intensity with no genotypes: the CytoScan 750K calls no chrY marker on
 * a male or a female alike, so the call-rate half of sexing.ts reads absence on every array. Where
 * that is so, and only on bulk DNA, chrY is read against the array's OWN single-copy level rather
 * than an absolute cut. A male's Y sits beside his X (median +0.09 log2 over 26 normal males); an
 * absent Y sits 1.8 to 2.3 below it (20 stated 45,X and 33 normal females). An absolute threshold
 * taken from another platform lands inside the male cluster here, which called four normal males a
 * monosomy; the array's own X does not move between platforms.
 *
 * ON DIPLOID MATERIAL ONLY. A haploid product carries one of everything, so one X is what it is
 * supposed to have; asking this question of a pronucleus or a sperm would report normal
 * segregation as a monosomy. The copy count is relative to the sample's own autosomes, so a
 * triploid 69,XXX reads as two X here. Triploidy is a whole-genome state and is reported by its
 * own channel.
 */
import { intensityDetects } from './intensityNull.ts'
import { COPY_SHIFT_FLOOR } from './parentage.ts'
import { sexCall, type SexTally } from './sexing.ts'

/**
 * Fewest chrY markers for intensity to answer on its own, on a panel that genotypes none of them.
 *
 * A file carrying a handful of chrY probes is a panel too thin to ask. The platforms this path
 * exists for carry thousands: the CytoScan 750K writes 4,135.
 */
export const Y_CN_MIN_PROBES = 200

/**
 * How far below the array's own chrX level a chrY may sit and still be a chrY.
 *
 * Measured on the CytoScan miscarriage series, chrY minus chrX: +0.09 median on 26 normal males,
 * -1.78 on 20 stated 45,X, -2.25 on 33 normal females. The cut sits in the middle of that gap with
 * about a log2 of margin on each side.
 */
export const Y_BESIDE_X = -0.9

/** What the sex chromosomes were found to be. */
export type Constitution = 'XX' | 'XY' | '45,X' | '47,XXY' | '47,XXX' | 'unresolved'

export interface SexChromosomeCall {
  /** X copies relative to the sample's own autosomes. Null where nothing could be measured. */
  xCopies: 1 | 2 | 3 | null
  yBearing: boolean | null
  constitution: Constitution
  /** True only where the constitution named is itself an abnormality. */
  aneuploid: boolean
  why: string
}

export interface SexChromosomeInput {
  /** chrX median log2 ratio against the sample's autosomal centre. */
  xShift: number
  /** That shift over the array's own autosomal spread, where one could be built. */
  xZ?: number
  /** The sample's own chrY tally, read by sexing.ts. */
  y: SexTally
  /** Whether this array's intensity can be read at face value: bulk DNA at a bulk call rate. */
  bulk?: boolean
  /** Whether this material is expected to carry two of each chromosome. */
  diploid: boolean
}

/**
 * Whether this sample carries a chrY, including on a panel that genotypes none of it.
 *
 * Thousands of chrY markers and not one call is a property of the PLATFORM, not of the sample: an
 * absent chromosome still produces calls, so a panel that genotypes chrY at all produces some on a
 * female too. Where that is what the file shows, the call rate says nothing and only intensity is
 * left, which is readable on bulk DNA and not on an amplified single cell. Null is returned rather
 * than a guess in the case that remains.
 */
export function yPresent(o: SexChromosomeInput): boolean | null {
  const call = sexCall(o.y)
  const panelGenotypesY = o.y.yCalled > 0 || o.y.yTotal < Y_CN_MIN_PROBES
  if (panelGenotypesY) return call.yBearing
  if (o.bulk === true && Number.isFinite(call.lrrShift) && Number.isFinite(o.xShift)) {
    return call.lrrShift - o.xShift >= Y_BESIDE_X
  }
  return null
}

export function callSexChromosomes(input: SexChromosomeInput): SexChromosomeCall {
  const o = { ...input, yBearing: yPresent(input) }
  const none = (why: string): SexChromosomeCall => ({
    xCopies: null, yBearing: o.yBearing, constitution: 'unresolved', aneuploid: false, why,
  })
  if (!o.diploid) {
    return none('this material is not expected to carry two of each chromosome, so one X is '
      + 'ordinary rather than evidence of anything')
  }
  if (!Number.isFinite(o.xShift)) {
    return none('chromosome X carries no intensity on this file, so its copy number was not '
      + 'measured')
  }
  // The same two gates the autosomes are held to: the shift is bigger than this array's own
  // noise, and it is big enough to be a chromosome.
  const detects = (o.xZ === undefined || intensityDetects(o.xZ, true))
    && Math.abs(o.xShift) >= COPY_SHIFT_FLOOR
  const xCopies = !detects ? 2 : o.xShift < 0 ? 1 : 3
  const level = `chrX sits ${o.xShift.toFixed(2)} log2 from this sample's autosomes`

  if (xCopies === 1) {
    if (o.yBearing === true) {
      return {
        xCopies,
        yBearing: o.yBearing,
        constitution: 'XY',
        aneuploid: false,
        why: `${level}, which is one copy, and chrY is present. One X and a Y is an ordinary male `
          + 'karyotype',
      }
    }
    if (o.yBearing === false) {
      return {
        xCopies,
        yBearing: o.yBearing,
        constitution: '45,X',
        aneuploid: true,
        why: `${level}, which is one copy, and chrY is measurable on this panel and absent. One X `
          + 'and no Y is a monosomy, not a male karyotype',
      }
    }
    return {
      xCopies,
      yBearing: null,
      constitution: 'unresolved',
      aneuploid: false,
      why: `${level}, which is one copy. Whether that is a male karyotype or a monosomy X turns `
        + 'entirely on chrY, and this file carries too few chrY markers to say. The two are not '
        + 'separable here',
    }
  }
  if (xCopies === 3) {
    return {
      xCopies,
      yBearing: o.yBearing,
      constitution: o.yBearing === true ? '47,XXY' : '47,XXX',
      aneuploid: true,
      why: `${level}, which is more than two copies${o.yBearing === true
        ? ', with chrY present' : ''}`,
    }
  }
  return {
    xCopies,
    yBearing: o.yBearing,
    constitution: o.yBearing === true ? '47,XXY' : 'XX',
    aneuploid: o.yBearing === true,
    why: o.yBearing === true
      ? `${level}, which is two copies, and chrY is present as well`
      : `${level}, which is two copies`,
  }
}
