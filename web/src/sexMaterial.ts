/**
 * Whether this array shows two copies of its autosomes, which is what the level of chrX is read
 * against.
 *
 * chrX at the autosomes' level is two X in a diploid cell and one X in a haploid one: every
 * chromosome of one genome is at one copy. So intensity alone cannot tell a pronucleus from a
 * diploid XX cell, and a heterozygosity sex reads both as female, since on one genome every
 * heterozygous call is noise on X and the autosomes alike. What separates them is how many
 * genomes the array carries.
 *
 * With a parent loaded the stage reads that from the parent's alleles, and this keeps the stage's
 * answer. With no parent it has only heterozygosity, which drop-in inflates on one genome without
 * reaching the diploid ceiling, so a pronucleus can sit on a diploid rung. There the heterozygosity
 * must clear HAPLOID_MAX_HET with the stage's own drop-in removed: haploid pronuclei read 0.2% to
 * 1.4% so, and 10.5% is the line a diploid reaches at any stage. A diploid the stage read from its
 * BAF band, at or under that line, keeps that reading: the band is not heterozygosity.
 */
import type { Gate, Sex } from './ingest.ts'
import { HAPLOID_MAX_HET, type StageCall } from './stage.ts'

export interface SexMaterial {
  /** The array shows two copies of its autosomes, so chrX can be read as one X or two. */
  diploid: boolean
  haploid: boolean
  /** The stage resolved: not failed, not unknown. */
  staged: boolean
  /** Why not, where not. */
  why: string
}

const pct = (x: number) => `${(x * 100).toFixed(1)}%`

export function sexMaterial(o: {
  stage?: StageCall
  /** Heterozygous fraction of called markers, as read. */
  hetRate: number
  /** The founder-panel drop-in the stage was given, NaN where none was measured. */
  dropIn: number
  parentLoaded: boolean
}): SexMaterial {
  const s = o.stage
  if (!s || s.stage === 'failed' || s.stage === 'unknown') {
    return {
      diploid: false, haploid: false, staged: false,
      why: `the stage did not resolve (${s?.stage ?? 'not read'}), so how many genomes this array `
        + 'carries is not known, and chrX at the level of the autosomes could be one X in one '
        + 'genome or two X in two',
    }
  }
  if (s.stage === 'haploid' || s.ploidy === 'haploid') {
    return { diploid: false, haploid: true, staged: true, why: 'this cell is staged haploid: one genome' }
  }
  if (o.parentLoaded) return { diploid: true, haploid: false, staged: true, why: '' }
  const h = o.hetRate
  const d = o.dropIn
  const real = Number.isFinite(d) && d > 0 && d < 1 ? Math.max(0, (h - d) / (1 - d)) : h
  if (real > HAPLOID_MAX_HET) return { diploid: true, haploid: false, staged: true, why: '' }
  if (s.ploidy === 'diploid' && h <= HAPLOID_MAX_HET) {
    return { diploid: true, haploid: false, staged: true, why: '' }
  }
  return {
    diploid: false, haploid: false, staged: true,
    why: `${pct(h)} heterozygous as read${Number.isFinite(d) ? `, ${pct(real)} with this array's `
      + `${pct(d)} drop-in removed` : ''}, at or under the ${pct(HAPLOID_MAX_HET)} a diploid `
      + 'reaches at any stage. With no parent loaded nothing else shows a second genome, so two '
      + 'copies of the autosomes are not shown, and chrX at their level could be one X in one '
      + 'genome or two X in two',
  }
}

/**
 * The heterozygosity sex (ingest.ts callSex) as the page prints it for a sample read with no
 * parent: only where two copies are shown. Elsewhere it would read female on a pronucleus. With
 * chrY present, heterozygous calls on X are not a woman's (sexChromosome.ts: 47,XXY or a mix of
 * cells), with a parent or without.
 */
export function shownSex(sex: Sex, m: SexMaterial | undefined, yBearing?: boolean | null): string {
  if (sex === 'female' && yBearing === true && (!m || m.diploid)) {
    return 'not read: chrY present with a heterozygous X'
  }
  if (!m || m.diploid) return sex
  return m.haploid ? 'not read: staged haploid, one sex chromosome is its normal complement'
    : 'not read: two copies of the autosomes not shown'
}

/** The quality gates with the heterozygosity sex gate stated as not read, where two copies are not shown. */
export function sexGates(gates: Gate[], m: SexMaterial | undefined): Gate[] {
  if (!m || m.diploid) return gates
  return gates.map((x) => (x.name === 'sex call'
    ? { ...x, verdict: 'report_only' as const, detail: `Not read: ${m.why}.` } : x))
}
