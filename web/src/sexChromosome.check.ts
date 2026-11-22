// Self-check for the sex chromosome call. Run: node src/sexChromosome.check.ts
//
// The numbers are measured, not assumed: chrX median log2 against the sample's own autosomes over
// the CytoScan miscarriage series, whose karyotype is stated per sample. One X reads about -0.55
// whether the sample is a male or a monosomy, which is the whole reason chrY has to be consulted.
import assert from 'node:assert/strict'
import { callSexChromosomes } from './sexChromosome.ts'
import type { SexTally } from './sexing.ts'

/** A chrY tally on a panel that genotypes chrY, which is what Axiom arrays do. */
const genotypedY = (present: boolean): SexTally => ({
  yCalled: present ? 48 : 0,
  yTotal: 51,
  autoCalled: 190_000,
  autoTotal: 200_000,
  yLrr: Array(51).fill(present ? 0.30 : -1.05),
  autoLrr: Array(2000).fill(0),
  seen: 200_000,
})

/** A panel that carries chrY as intensity only: the CytoScan 750K genotypes none of its 4,135. */
const intensityOnlyY = (lrr: number): SexTally => ({
  yCalled: 0,
  yTotal: 4135,
  autoCalled: 190_000,
  autoTotal: 200_000,
  yLrr: Array(4135).fill(lrr),
  autoLrr: Array(2000).fill(0),
  seen: 200_000,
})

/** Too few chrY markers to ask at all. */
const thinY: SexTally = {
  yCalled: 0, yTotal: 4, autoCalled: 190_000, autoTotal: 200_000,
  yLrr: [0, 0, 0, 0], autoLrr: Array(2000).fill(0), seen: 200_000,
}

const diploid = { diploid: true }

// --- 1. ONE X IS NOT A DIAGNOSIS ON ITS OWN --------------------------------------------------
{
  // GSM6322010, stated 45,X: -0.570. A normal male of the same series: -0.551. The same number.
  const monosomy = callSexChromosomes({ ...diploid, xShift: -0.570, y: genotypedY(false) })
  const male = callSexChromosomes({ ...diploid, xShift: -0.551, y: genotypedY(true) })
  assert.equal(monosomy.constitution, '45,X')
  assert.equal(monosomy.aneuploid, true)
  assert.equal(male.constitution, 'XY')
  assert.equal(male.aneuploid, false)
  assert.equal(monosomy.xCopies, male.xCopies, 'the intensity is the same; only chrY differs')
}

// --- 2. A PANEL THAT CANNOT SEE chrY REFUSES -------------------------------------------------
{
  // The CytoScan 750K genotypes no chrY marker, so sexing.ts returns null rather than false.
  const c = callSexChromosomes({ ...diploid, xShift: -0.570, y: thinY })
  assert.equal(c.constitution, 'unresolved')
  assert.equal(c.xCopies, 1, 'the copy number is still measured and still reported')
  assert.equal(c.aneuploid, false)
  assert.ok(c.why.includes('chrY'), 'the reason names the measurement that is missing')
}

// --- 3. TWO X -------------------------------------------------------------------------------
{
  // Normal females of the same series run -0.05 to +0.01.
  assert.equal(callSexChromosomes({ ...diploid, xShift: -0.010, y: genotypedY(false) }).constitution,
    'XX')
  assert.equal(callSexChromosomes({ ...diploid, xShift: -0.054, y: thinY }).constitution,
    'XX', 'two X needs no chrY to be said')
  // A triploid 69,XXX carries three of everything, so chrX against its own autosomes is two.
  assert.equal(callSexChromosomes({ ...diploid, xShift: -0.040, y: thinY }).constitution,
    'XX')
  // Both sex chromosomes present is an abnormality in a way two X alone is not.
  // Bulk DNA, chrX at the autosomes' level, heterozygous on X as two X are, and a chrY.
  const kline = callSexChromosomes({
    ...diploid, bulk: true, hetSex: 'female', hetRatio: 1.0, xShift: -0.010, y: genotypedY(true),
  })
  assert.equal(kline.constitution, '47,XXY')
  assert.equal(kline.aneuploid, true)
}

// --- 4. HAPLOID MATERIAL IS NOT ASKED --------------------------------------------------------
{
  // A maternal pronucleus carries one X and that is what it is supposed to carry.
  const c = callSexChromosomes({ diploid: false, xShift: -0.570, y: genotypedY(false) })
  assert.equal(c.constitution, 'unresolved')
  assert.equal(c.aneuploid, false)
  assert.equal(c.xCopies, null)
}

// --- 5. NOISE IS NOT A COPY CHANGE -----------------------------------------------------------
{
  // Below the magnitude floor, and a z that says the shift is inside this array's own spread.
  assert.equal(callSexChromosomes({ ...diploid, xShift: -0.30, y: thinY }).xCopies, 2)
  assert.equal(
    callSexChromosomes({ ...diploid, xShift: -0.60, xZ: -1.2, y: genotypedY(false) }).constitution,
    'XX', 'a big shift on a noisy array is not a copy change')
  assert.equal(
    callSexChromosomes({ ...diploid, xShift: -0.60, xZ: -22, y: genotypedY(false) }).constitution,
    '45,X')
  assert.equal(callSexChromosomes({ ...diploid, xShift: NaN, y: genotypedY(true) }).constitution,
    'unresolved')
}

// --- 6. A PANEL THAT GENOTYPES NO chrY -------------------------------------------------------
{
  // Measured on the CytoScan miscarriage series: a normal male reads chrY at -0.46 with chrX at
  // -0.55, and a stated 45,X reads chrY at -2.35 with the same chrX. The absolute cut taken from
  // a panel that genotypes chrY sits inside the male cluster here, so the array's own chrX is the
  // reference instead.
  const male = callSexChromosomes({
    ...diploid, bulk: true, xShift: -0.551, y: intensityOnlyY(-0.464),
  })
  assert.equal(male.constitution, 'XY', 'a chrY beside the chrX is a chrY, whatever the genotyper did')
  const monosomy = callSexChromosomes({
    ...diploid, bulk: true, xShift: -0.570, y: intensityOnlyY(-2.348),
  })
  assert.equal(monosomy.constitution, '45,X')
  // Amplified material is not admitted: its intensity is what the two-signal rule exists for.
  assert.equal(callSexChromosomes({
    ...diploid, xShift: -0.551, y: intensityOnlyY(-0.464),
  }).constitution, 'unresolved')
  // And a female on the same platform is not a male just because chrX is where it should be.
  assert.equal(callSexChromosomes({
    ...diploid, bulk: true, xShift: -0.010, y: intensityOnlyY(-2.263),
  }).constitution, 'XX')
}

// --- 7. TWO X WITH A Y ON AMPLIFIED MATERIAL IS NOT A DEFAULT -----------------------------------
{
  // A single cell with a parent: chrX -0.30, inside the noise, a chrY, and no heterozygous calls on X.
  const xy = callSexChromosomes({ ...diploid, hetSex: 'male', hetRatio: 0.05, xShift: -0.30, xZ: -1.8,
    y: genotypedY(true) })
  assert.equal(xy.constitution, 'XY', xy.why)
  // The same level with heterozygosity between the bands: not 47,XXY, and the reason says why.
  const weak = callSexChromosomes({ ...diploid, hetSex: 'ambiguous', hetRatio: 0.3, xShift: 0.17, xZ: 0.5,
    y: genotypedY(true) })
  assert.equal(weak.constitution, 'unresolved')
  assert.ok(weak.why.includes('cannot separate one X from two'), weak.why)
  // Bulk DNA at the autosomes' level whose heterozygous calls read one X: XY, not 47,XXY.
  assert.equal(callSexChromosomes({ ...diploid, bulk: true, hetSex: 'male', hetRatio: 0.05, xShift: -0.12,
    xZ: -2.1, y: genotypedY(true) }).constitution, 'XY')
  // Two copy steps down and still called: not read.
  assert.equal(callSexChromosomes({ ...diploid, hetSex: 'ambiguous', xShift: -1.77, xZ: -4.5,
    xCallFraction: 0.9, y: genotypedY(true) }).constitution, 'unresolved')
}

// --- 8. ON A CHIP THAT GENOTYPES chrY, A WOMAN'S Y MARKERS READ NO Y -------------------------------
{
  // A woman's bulk array on this chip: 0 of 812 chrY markers called, chrY intensity near chrX, half
  // of the chrY markers at a pole where her own homozygous calls sit there 97% of the time.
  const woman: SexTally = { ...intensityOnlyY(-0.12), yTotal: 812, yLrr: Array(812).fill(-0.12),
    yBaf: 812, yBafPole: 406, homBaf: 18_000, homBafPole: 17_460 }
  const c = callSexChromosomes({ ...diploid, bulk: true, chipGenotypesY: true, hetSex: 'female', hetRatio: 1.0,
    xShift: 0.19, xZ: 1.0, y: woman })
  assert.equal(c.yBearing, false)
  assert.equal(c.constitution, 'XX', c.why)
  // The same file with nothing said of the chip still takes the intensity-only route, which is
  // what a chip that never genotypes chrY needs (section 6).
  assert.equal(callSexChromosomes({ ...diploid, bulk: true, xShift: -0.551, y: intensityOnlyY(-0.464) })
    .constitution, 'XY')
}

console.log('sexChromosome.check.ts OK: one X is a male or a monosomy, and chrY is what decides')
