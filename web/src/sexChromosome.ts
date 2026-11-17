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
import { CALL_COLLAPSE, COPY_SHIFT_FLOOR } from './parentage.ts'
import type { Sex } from './ingest.ts'
import type { SexMaterial } from './sexMaterial.ts'
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
export type Constitution = 'XX' | 'XY' | '45,X' | '47,XXY' | '47,XXX' | 'one X' | 'XX or one X'
  | 'unresolved'

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
  /** Why it is not, from sexMaterial.ts; a haploid cell gets its own reading. */
  material?: SexMaterial
  /** Why chrX has no level, where it has none: what the file lacks. */
  xMissing?: string
  /** Two X only where chrX is inside the array's own noise; set with no parent loaded. */
  twoInsideNoise?: boolean
  /** The heterozygosity sex (ingest.ts callSex) and its ratio, from the sample's own calls. */
  hetSex?: Sex
  hetRatio?: number | null
  /** chrX call rate over the genome's median: under CALL_COLLAPSE it has no template. */
  xCallFraction?: number
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

/**
 * With no parent, on material shown to hold two copies, heterozygous calls on chrX are a second
 * reading of the X count: real on two genomes, which is the only material this is asked of. The
 * bands are callSex's. A reading is printed only where nothing measured contradicts it.
 */
export function callSexChromosomes(input: SexChromosomeInput): SexChromosomeCall {
  const c = byIntensity(input)
  if (!input.hetSex || !input.twoInsideNoise || !input.diploid) return c
  const het = input.hetSex === 'female' ? 2 : input.hetSex === 'male' ? 1 : null
  if (het === null) return c
  const ratio = input.hetRatio == null ? '' : ` ${input.hetRatio.toFixed(2)}`
  const hetWhy = `heterozygous calls on chrX read ${het === 2 ? 'two X' : 'one X'} (chrX to autosome `
    + `heterozygosity ratio${ratio}, in the ${input.hetSex} band)`
  const y = c.yBearing
  const no = (why: string): SexChromosomeCall => ({
    xCopies: null, yBearing: y, constitution: 'unresolved', aneuploid: false, why,
  })
  if (c.xCopies !== null && c.xCopies !== het) {
    return no(`the readings disagree: ${hetWhy}; by intensity, ${c.why}`)
  }
  if (het === 2 && y === true) return no(`the readings disagree: ${hetWhy}, and chrY is present`)
  if (c.constitution !== 'unresolved') return { ...c, why: `by intensity, ${c.why}; and ${hetWhy}` }
  if (het === 2) {
    return { xCopies: 2, yBearing: y, constitution: 'XX', aneuploid: false,
      why: `${hetWhy}; intensity alone did not separate one X from two (${c.why})` }
  }
  // Two copy steps under is not one X, however few heterozygous calls X carries.
  if (y === true && !(input.xShift <= -2 * COPY_SHIFT_FLOOR)) {
    return { xCopies: 1, yBearing: y, constitution: 'XY', aneuploid: false,
      why: `${hetWhy}, and chrY is present by its calls and intensity` }
  }
  return no(`${hetWhy}, and ${y === false ? 'chrY is genotyped on this panel and not present'
    : 'this file cannot say whether chrY is present'}: a male karyotype and a monosomy X are not `
    + 'separable from heterozygosity')
}

function byIntensity(input: SexChromosomeInput): SexChromosomeCall {
  const o = { ...input, yBearing: yPresent(input) }
  const none = (why: string): SexChromosomeCall => ({
    xCopies: null, yBearing: o.yBearing, constitution: 'unresolved', aneuploid: false, why,
  })
  const xWhy = o.xMissing ?? 'chromosome X carries no intensity on this file'
  // The same two gates the autosomes are held to: the shift is bigger than this array's own
  // noise, and it is big enough to be a chromosome.
  const detects = (o.xZ === undefined || intensityDetects(o.xZ, true))
    && Math.abs(o.xShift) >= COPY_SHIFT_FLOOR
  const level = `chrX sits ${o.xShift.toFixed(2)} log2 from this sample's autosomes`
  if (o.twoInsideNoise && o.material?.staged) return whatXAndYShow(o, level, detects)
  if (!o.diploid && o.material?.haploid) {
    // One genome: one X at the autosomes' level, or none and a Y. Neither is a sex.
    const x = !Number.isFinite(o.xShift) ? `${xWhy}, so where chrX sits was not measured`
      : !detects ? `${level}, which in one genome is one X`
        : o.xShift < 0 ? `${level}, under them by more than a copy step, which in one genome is no X`
          : `${level}, over them by more than a copy step, which is more than one X in one genome`
    const y = o.yBearing === true ? 'chrY is present by its calls and intensity'
      : o.yBearing === false ? 'chrY is genotyped on this panel and not present'
        : 'this file cannot say whether chrY is present'
    return none(`${o.material.why}, so one sex chromosome is its normal complement and says `
      + `nothing of a sex. ${x}; ${y}`)
  }
  if (!o.diploid) {
    return none(o.material?.why ?? 'this material is not expected to carry two of each chromosome, '
      + 'so one X is ordinary rather than evidence of anything')
  }
  if (!Number.isFinite(o.xShift)) {
    return none(`${xWhy}, so its copy number was not measured`)
  }
  // Beyond the noise and short of the copy step is neither level: on amplified single cells one X
  // reads -0.3 and the copy step is 0.40, so "not a copy step" is not "two copies".
  if (o.twoInsideNoise && !detects) {
    const zOut = o.xZ !== undefined && intensityDetects(o.xZ, true)
    if (zOut || Math.abs(o.xShift) >= COPY_SHIFT_FLOOR) {
      return none(`${level}, ${zOut ? "beyond this array's own noise but short of the "
        + `${COPY_SHIFT_FLOOR.toFixed(2)} log2 copy step` : `past the ${COPY_SHIFT_FLOOR.toFixed(2)} log2 `
        + "copy step but inside this array's own noise"}, so one X or two is not separable here`)
    }
    if (o.yBearing === true) {
      return none(`${level}, not a one-copy step, and chrY is present. A one-X level compressed by `
        + 'amplification and 47,XXY are not separable without a parent')
    }
  }
  const xCopies = !detects ? 2 : o.xShift < 0 ? 1 : 3

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
  // TWO X WITH A Y IS A CLAIM, NOT A DEFAULT. Two X here only means chrX was not shown a copy step
  // under the autosomes, which on amplified material an XY cell often is not. 47,XXY stands where
  // bulk DNA holds chrX at its autosomes' level and heterozygous calls on X do not read one X.
  if (!o.twoInsideNoise && o.yBearing === true) {
    const twoSteps = o.xShift <= -2 * COPY_SHIFT_FLOOR
    if (o.hetSex === 'male' && !twoSteps) {
      const ratio = o.hetRatio == null ? '' : ` ${o.hetRatio.toFixed(2)}`
      return {
        xCopies: 1, yBearing: true, constitution: 'XY', aneuploid: false,
        why: `heterozygous calls on chrX read one X (chrX to autosome heterozygosity ratio${ratio}, `
          + `in the male band), and chrY is present by its calls and intensity; ${level}, not a `
          + 'copy step, which on this material does not separate one X from two',
      }
    }
    const atLevel = Math.abs(o.xShift) < COPY_SHIFT_FLOOR
      && !(o.xZ !== undefined && intensityDetects(o.xZ, true))
    if (!(o.bulk === true && atLevel)) {
      const called = !(o.xCallFraction !== undefined && o.xCallFraction < CALL_COLLAPSE)
      return {
        xCopies: null, yBearing: true, constitution: 'unresolved', aneuploid: false,
        why: twoSteps && called
          ? `${level}: chromosome X sits further below the autosomes than one copy, and its calls `
            + 'do not show it absent either, so this file does not settle what the cell holds'
          : `${level}, not a copy step, and chrY is present. On ${o.bulk === true
            ? 'bulk DNA with chrX away from the level of the autosomes' : 'amplified material'} this `
            + 'X level cannot separate one X from two, so XY and 47,XXY are not separable here',
      }
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

/**
 * With no parent, what chrX and chrY themselves show, with the ploidy added only where it is read.
 * One genome holds chrX at its autosomes' level or holds none, so chrX present a copy step under
 * the autosomes with a chrY is two genomes, XY, whatever the heterozygosity says. chrX at the
 * autosomes' level with no chrY is XX in two genomes and one X in one, which intensity cannot
 * separate. No aneuploid constitution is printed with no parent.
 */
function whatXAndYShow(
  o: SexChromosomeInput & { yBearing: boolean | null }, level: string, detects: boolean,
): SexChromosomeCall {
  const m = o.material!
  const none = (why: string): SexChromosomeCall => ({
    xCopies: null, yBearing: o.yBearing, constitution: 'unresolved', aneuploid: false, why,
  })
  const say = (constitution: Constitution, xCopies: 1 | 2 | null, why: string): SexChromosomeCall => ({
    xCopies, yBearing: o.yBearing, constitution, aneuploid: false, why,
  })
  if (!Number.isFinite(o.xShift)) {
    return none(`${o.xMissing ?? 'chromosome X carries no intensity on this file'}, so its level was `
      + 'not measured')
  }
  const zOut = o.xZ !== undefined && intensityDetects(o.xZ, true)
  const step = COPY_SHIFT_FLOOR.toFixed(2)
  const yText = o.yBearing === true ? 'chrY is present by its calls and intensity'
    : o.yBearing === false ? 'chrY is genotyped on this panel and not present'
      : 'this file cannot say whether chrY is present'
  if (!detects && (zOut || Math.abs(o.xShift) >= COPY_SHIFT_FLOOR)) {
    return none(`${level}, ${zOut ? `beyond this array's own noise but short of the ${step} log2 copy step`
      : `past the ${step} log2 copy step but inside this array's own noise`}, so one X or two is not `
      + 'separable here')
  }
  if (detects && o.xShift > 0) {
    return none(`${level}, over the autosomes by more than a copy step. An extra X is not printed with `
      + 'no parent loaded')
  }
  if (detects) {
    const cf = o.xCallFraction
    if (cf !== undefined && Number.isFinite(cf) && cf < CALL_COLLAPSE) {
      return none(`${level}, and chrX calls at ${(cf * 100).toFixed(0)}% of the genome's median rate, `
        + `under the ${(CALL_COLLAPSE * 100).toFixed(0)}% a chromosome that is there reaches: no X `
        + `template. ${yText}`)
    }
    if (o.xShift <= -2 * COPY_SHIFT_FLOOR) {
      return none(`${level}: chromosome X sits further below the autosomes than one copy, and its calls `
        + 'do not show it absent either, so this file does not settle what the cell holds')
    }
    if (o.yBearing === true) {
      return say('XY', 1, `${level}, a copy step under the autosomes, and ${yText}. One genome cannot `
        + 'hold chrX a step under its autosomes, so this is two genomes: XY')
    }
    return none(`${level}, a copy step under the autosomes, and ${yText}. With no parent loaded a `
      + 'missing X is not printed: a male karyotype and a monosomy X are not separated here')
  }
  if (o.yBearing === true) {
    return none(`${level}, not a one-copy step, and chrY is present. A one-X level compressed by `
      + 'amplification and 47,XXY are not separable without a parent')
  }
  if (m.haploid) {
    return say('one X', 1, `the normal complement of one genome: this cell is staged haploid and `
      + `${level}, the level of its autosomes`)
  }
  if (m.diploid) {
    return say('XX', 2, `${level}, which is two copies, and two copies of the autosomes are shown`)
  }
  return say('XX or one X', null, 'X at the level of the autosomes, no Y read: XX if this cell holds '
    + `two genomes, one X if it holds one, and this file does not settle which (${level}; ${m.why})`)
}
