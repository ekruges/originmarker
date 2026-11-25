// Self-check for the sex reading with no parent. Run: node src/sexMaterial.check.ts
//
// chrX at the autosomes' level is two X in two genomes and one X in one, so the sex is read only
// where two copies are shown; a pronucleus whose heterozygosity is drop-in must not read XX.
import assert from 'node:assert/strict'
import { classify, emptyTally } from './parentage.ts'
import { callSexChromosomes, sexRow, type SexChromosomeCall } from './sexChromosome.ts'
import { sexGates, sexMaterial, shownSex } from './sexMaterial.ts'
import type { SexTally } from './sexing.ts'
import type { StageCall } from './stage.ts'

const stage = (s: Partial<StageCall>) => ({ stage: 'single-cell', ploidy: 'diploid', ...s }) as StageCall
const noY: SexTally = {
  yCalled: 0, yTotal: 51, autoCalled: 190_000, autoTotal: 200_000,
  yLrr: Array(51).fill(-1.05), autoLrr: Array(2000).fill(0), seen: 200_000,
}

// --- 1. TWO COPIES SHOWN, OR NOT ---------------------------------------------------------------
{
  const np = { parentLoaded: false }
  // A pronucleus on a diploid rung: 15% heterozygous, 14% of it drop-in. One genome.
  const pn = sexMaterial({ ...np, stage: stage({}), hetRate: 0.15, dropIn: 0.14 })
  assert.equal(pn.diploid, false)
  assert.equal(pn.haploid, false)
  assert.ok(pn.why.includes('drop-in removed'), 'the reason quotes the numbers it rests on')
  // A diploid cell: 16% heterozygous, 2% drop-in.
  assert.equal(sexMaterial({ ...np, stage: stage({}), hetRate: 0.16, dropIn: 0.02 }).diploid, true)
  // Ploidy unknown under heavy drop-in, real heterozygosity well over one genome's.
  assert.equal(sexMaterial({ ...np, stage: stage({ stage: 'trophectoderm', ploidy: 'unknown' }),
    hetRate: 0.30, dropIn: 0.15 }).diploid, true)
  // Read diploid on the BAF band under the haploid line: the band is not heterozygosity.
  assert.equal(sexMaterial({ ...np, stage: stage({ stage: 'blastomere' }), hetRate: 0.08, dropIn: 0.01 })
    .diploid, true)
  assert.equal(sexMaterial({ ...np, stage: stage({ stage: 'haploid', ploidy: 'haploid' }),
    hetRate: 0.01, dropIn: 0.01 }).haploid, true)
  for (const s of ['failed', 'unknown'] as const) {
    const m = sexMaterial({ ...np, stage: stage({ stage: s, ploidy: 'unknown' }), hetRate: 0.16, dropIn: 0.02 })
    assert.equal(m.diploid, false)
    assert.ok(m.why.includes('did not resolve'))
  }
  // With a parent the stage's answer stands: it read the genome count from the alleles.
  assert.equal(sexMaterial({ parentLoaded: true, stage: stage({}), hetRate: 0.15, dropIn: 0.14 })
    .diploid, true)
}

// --- 2. NO SEX PRINTED WHERE TWO COPIES ARE NOT SHOWN --------------------------------------------
{
  const pn = sexMaterial({ parentLoaded: false, stage: stage({}), hetRate: 0.15, dropIn: 0.14 })
  const c = callSexChromosomes({ diploid: pn.diploid, material: pn, xShift: 0.02, y: noY })
  assert.equal(c.constitution, 'unresolved', 'chrX at the autosomes is not XX on one genome')
  assert.equal(c.why, pn.why)
  assert.equal(shownSex('female', pn), 'not read: two copies of the autosomes not shown')
  assert.equal(shownSex('female', undefined), 'female', 'with a parent the page is unchanged')
  const g = sexGates([{ name: 'sex call', verdict: 'usable', value: 1, detail: 'female band' }], pn)
  assert.equal(g[0].verdict, 'report_only')
  assert.ok(!g[0].detail.includes('female'))

  const hap = sexMaterial({ parentLoaded: false, stage: stage({ stage: 'haploid', ploidy: 'haploid' }),
    hetRate: 0.01, dropIn: NaN })
  const one = callSexChromosomes({ diploid: false, material: hap, xShift: 0.02, y: noY })
  assert.equal(one.constitution, 'unresolved')
  assert.ok(one.why.includes('in one genome is one X'), one.why)
  const none = callSexChromosomes({ diploid: false, material: hap, xShift: -1.1, xZ: -20, y: noY })
  assert.ok(none.why.includes('in one genome is no X'), none.why)
  assert.ok(shownSex('female', hap).startsWith('not read'))
}

// --- 2b. BEYOND THE NOISE AND SHORT OF A COPY STEP IS NOT TWO X ---------------------------------
{
  const yes: SexTally = { ...noY, yCalled: 48, yLrr: Array(51).fill(0.3) }
  // A normal male single cell with no parent: chrX -0.31 at z -6.7 and a chrY. Not 47,XXY.
  const c = callSexChromosomes({ diploid: true, twoInsideNoise: true, xShift: -0.31, xZ: -6.7, y: yes })
  assert.equal(c.constitution, 'unresolved', c.why)
  // Inside both gates with a chrY: not 47,XXY with no parent.
  assert.equal(callSexChromosomes({ diploid: true, twoInsideNoise: true, xShift: -0.32, xZ: -0.8, y: yes })
    .constitution, 'unresolved')
  // Past the copy step on a noisy array is not two X.
  assert.equal(callSexChromosomes({ diploid: true, twoInsideNoise: true, xShift: -0.58, xZ: -3.7, y: noY })
    .constitution, 'unresolved')
  assert.equal(callSexChromosomes({ diploid: true, twoInsideNoise: true, xShift: 0.08, xZ: 3.4, y: noY })
    .constitution, 'XX', 'inside the noise two X stands')
  // With a parent loaded the 6.5.0 reading stands, but not two X by default with a chrY on
  // amplified material.
  assert.equal(callSexChromosomes({ diploid: true, xShift: -0.31, xZ: -6.7, y: yes }).constitution, 'unresolved')
}

// --- 2c. HETEROZYGOSITY ON X, WITH NO PARENT, ON TWO COPIES ONLY ----------------------------------
{
  const yes: SexTally = { ...noY, yCalled: 48, yLrr: Array(51).fill(0.3) }
  const np = { diploid: true, twoInsideNoise: true }
  // Intensity between levels, heterozygous X: XX, and the reason says which evidence.
  const xx = callSexChromosomes({ ...np, hetSex: 'female', hetRatio: 1.1, xShift: 0.32, xZ: 6, y: noY })
  assert.equal(xx.constitution, 'XX')
  assert.ok(xx.why.startsWith('heterozygous calls on chrX read two X'), xx.why)
  // Heterozygous X and a chrY: 47,XXY or a mix of cells, which one array does not separate.
  const bad = callSexChromosomes({ ...np, hetSex: 'female', hetRatio: 1.1, xShift: -0.32, xZ: -0.8, y: yes })
  assert.equal(bad.constitution, '47,XXY or a mix of cells')
  assert.ok(bad.why.startsWith('chrY is present by its calls and intensity; chrX reads two different X'), bad.why)
  // No heterozygous X and a chrY: XY.
  assert.equal(callSexChromosomes({ ...np, hetSex: 'male', hetRatio: 0.05, xShift: -0.32, xZ: -0.8, y: yes })
    .constitution, 'XY')
  // Intensity one copy, heterozygosity two: disagree.
  assert.equal(callSexChromosomes({ ...np, hetSex: 'female', hetRatio: 1, xShift: -0.6, xZ: -20, y: noY })
    .constitution, 'unresolved')
  // Not shown two copies, or a parent loaded: heterozygosity is not consulted.
  assert.equal(callSexChromosomes({ diploid: false, twoInsideNoise: true, hetSex: 'female', xShift: 0, y: noY })
    .constitution, 'unresolved')
  // With a parent, a chrY beside heterozygous calls on X is 47,XXY or a mix of cells (sexing.check.ts 6).
  assert.equal(callSexChromosomes({ diploid: true, hetSex: 'female', xShift: -0.31, xZ: -6.7, y: yes })
    .constitution, '47,XXY or a mix of cells')
}

// --- 2d. WHAT X AND Y SHOW, WITH NO PARENT -------------------------------------------------------
{
  const yes: SexTally = { ...noY, yCalled: 48, yLrr: Array(51).fill(0.3) }
  const shown = sexMaterial({ parentLoaded: false, stage: stage({}), hetRate: 0.16, dropIn: 0.02 })
  const notShown = sexMaterial({ parentLoaded: false, stage: stage({}), hetRate: 0.15, dropIn: 0.14 })
  const hap = sexMaterial({ parentLoaded: false, stage: stage({ stage: 'haploid', ploidy: 'haploid' }),
    hetRate: 0.01, dropIn: NaN })
  const np = (m: typeof shown, xShift: number, xZ: number, y: SexTally, xCallFraction = 1) =>
    callSexChromosomes({ diploid: m.diploid, material: m, twoInsideNoise: true, xShift, xZ, y, xCallFraction })
  // X present a step under with a Y is two genomes whatever the heterozygosity says.
  assert.equal(np(notShown, -0.45, -12, yes).constitution, 'XY')
  // Two copy steps under and still called: not XY, by either route.
  const far = np(notShown, -1.69, -9.9, yes, 0.74)
  assert.equal(far.constitution, 'unresolved')
  assert.ok(far.why.includes('does not settle what the cell holds'), far.why)
  assert.equal(callSexChromosomes({ diploid: true, material: shown, twoInsideNoise: true, hetSex: 'male',
    hetRatio: 0.05, xShift: -1.94, xZ: -19, y: yes, xCallFraction: 0.96 }).constitution, 'unresolved')
  // No X template and a Y: not XY.
  assert.equal(np(hap, -1.9, -40, yes, 0.3).constitution, 'unresolved')
  // At the autosomes' level, no Y: one X on one genome, XX on two, the conditional line otherwise.
  assert.equal(np(hap, 0.08, 2, noY).constitution, 'one X')
  assert.equal(np(shown, 0.08, 2, noY).constitution, 'XX')
  const cond = np(notShown, 0.08, 2, noY)
  assert.equal(cond.constitution, 'XX or one X')
  assert.ok(cond.why.startsWith('X at the level of the autosomes, no Y read'), cond.why)
  // A step under with no Y read, and over the autosomes: nothing aneuploid with no parent.
  assert.equal(np(shown, -0.6, -20, noY).constitution, 'unresolved')
  assert.equal(np(shown, 0.6, 20, noY).constitution, 'unresolved')
  // Stage unresolved: not read.
  const failed = sexMaterial({ parentLoaded: false, stage: stage({ stage: 'failed', ploidy: 'unknown' }),
    hetRate: 0.16, dropIn: 0.02 })
  assert.equal(np(failed, 0.08, 2, noY).constitution, 'unresolved')
}

// --- 3. THE REASON FOR A MISSING LEVEL IS WHAT THE FILE LACKS ------------------------------------
{
  const c = callSexChromosomes({ diploid: true, xShift: NaN, y: noY,
    xMissing: 'this file has fewer than 200 chromosome X markers outside the pseudoautosomal region' })
  assert.ok(c.why.startsWith('this file has fewer than 200'), c.why)
}

// --- 4. A chrX ROW WITH NO PARENT ----------------------------------------------------------------
{
  // Call rate and intensity on 22 autosomes and chrX, no parental marker anywhere.
  const t = emptyTally()
  t.noParent = true
  for (const c of [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X']) {
    t.callByChrom.set(c, [270, 300])
    t.lrrByChrom.set(c, Array.from({ length: 300 }, (_, i) => (c === 'X' ? 0.3 : 0) + (i % 7) * 0.01))
  }
  const r = classify(t, NaN)
  const x = r.chroms.find((c) => c.chrom === 'X')
  assert.ok(x, 'chrX has a row with no parent')
  assert.ok(Math.abs(x.lrrShift - 0.3) < 1e-9, 'and its level against the autosomes')
  assert.equal(x.aneuploidy, undefined, 'the copy call stays autosomal')
  // With a parent loaded the chrX row comes from its markers only, as in 6.5.0.
  t.noParent = false
  assert.equal(classify(t, NaN).chroms.find((c) => c.chrom === 'X'), undefined)
}

// --- 5. THE SEX-CHROMOSOME ROW OF THE CHANGE LIST -----------------------------------------------
{
  const yes: SexTally = { ...noY, yCalled: 48, yLrr: Array(51).fill(0.3) }
  const shown = sexMaterial({ parentLoaded: false, stage: stage({}), hetRate: 0.16, dropIn: 0.02 })
  const notShown = sexMaterial({ parentLoaded: false, stage: stage({}), hetRate: 0.15, dropIn: 0.14 })
  const hap = sexMaterial({ parentLoaded: false, stage: stage({ stage: 'haploid', ploidy: 'haploid' }),
    hetRate: 0.01, dropIn: NaN })
  const read = (m: typeof shown, y: SexTally, xCallFraction: number) => callSexChromosomes({
    diploid: m.diploid, material: m, twoInsideNoise: true, xShift: -1.9, xZ: -20, y, xCallFraction })
  // No X template: the reading says so, and the ploidy decides what it is.
  const gone = read(notShown, yes, 0.38)
  assert.equal(gone.xNone, true, gone.why)
  assert.equal(gone.constitution, 'unresolved')
  assert.equal(read(notShown, yes, 0.74).xNone, undefined, 'calls still made are not no X')
  assert.equal(sexRow(read(shown, yes, 0.38), shown)?.dir, 'loss')
  assert.equal(sexRow(read(shown, yes, 0.38), shown)?.says, 'no X')
  assert.equal(sexRow(gone, notShown)?.dir, 'none', 'one Y genome or two with no X: a finding')
  assert.equal(sexRow(read(hap, yes, 0.38), hap), null, 'a haploid Y genome has no X by nature')
  for (const m of [shown, notShown, hap]) {
    assert.equal(sexRow(read(m, noY, 0.4), m)?.says, 'no X and no Y')
  }
  // A constitution stated as a fact is a change; one naming two possibilities is a finding.
  const as = (constitution: SexChromosomeCall['constitution']): SexChromosomeCall =>
    ({ xCopies: null, yBearing: null, constitution, aneuploid: false, why: 'x' })
  assert.equal(sexRow(as('45,X'))?.dir, 'loss')
  assert.equal(sexRow(as('47,XXX'))?.dir, 'gain')
  assert.equal(sexRow(as('47,XXY'))?.dir, 'gain')
  assert.equal(sexRow(as('47,XXY or a mix of cells'))?.dir, 'none')
  for (const c of ['XX', 'XY', 'one X', 'XX or one X', 'unresolved'] as const) assert.equal(sexRow(as(c)), null, c)
}

// --- 6. THE SEX LINE, THE SEX WORD, THE BADGE AND THE X ROW SAY ONE THING -------------------------
{
  const { sexSays } = await import('./sexMaterial.ts')
  const yes: SexTally = { ...noY, yCalled: 48, yLrr: Array(51).fill(0.3) }
  const shown = sexMaterial({ parentLoaded: false, stage: stage({}), hetRate: 0.16, dropIn: 0.02 })
  const notShown = sexMaterial({ parentLoaded: false, stage: stage({}), hetRate: 0.15, dropIn: 0.14 })
  const hap = sexMaterial({ parentLoaded: false, stage: stage({ stage: 'haploid', ploidy: 'haploid' }), hetRate: 0.01, dropIn: NaN })
  const read = (m: typeof shown, y: SexTally, xShift = -1.9, xCallFraction = 0.4) => callSexChromosomes({
    diploid: m.diploid, material: m, twoInsideNoise: true, xShift, xZ: -20, y, xCallFraction })
  // No X and no Y, on two copies not shown and on a haploid cell alike: one statement everywhere.
  for (const m of [notShown, hap]) {
    const c = read(m, noY)
    const say = sexSays('ambiguous', m, c)
    assert.ok(say.line.startsWith('No X and no Y found: '), say.line)
    assert.equal(say.word, 'no sex named: no sex chromosome was found')
    assert.equal(say.badge, 'no X, no Y')
    assert.equal(sexRow(c, m)?.dir, 'loss')
  }
  assert.equal(sexSays('ambiguous', shown, read(shown, yes)).word, 'no sex named: no X was found')
  const unsure = sexSays('ambiguous', notShown, read(notShown, yes))
  assert.ok(unsure.line.startsWith('No X found; ') && unsure.word.startsWith('not read: no X was found') && unsure.badge === null)
  // An XY the reading prints on two copies not shown is the word it names; not reportable, no badge.
  const xy = read(notShown, yes, -0.45, 1)
  assert.equal(xy.constitution, 'XY')
  assert.equal(sexSays('ambiguous', notShown, xy).word, 'male')
  assert.equal(sexSays('ambiguous', notShown, xy).badge, 'XY')
  assert.equal(sexSays('ambiguous', notShown, xy, false).badge, null)
  // Every reading keeps line, word and badge in step with its row.
  for (const [m, y, x, cf] of [[shown, noY, -1.9, 0.4], [shown, yes, -1.9, 0.4], [notShown, yes, -1.9, 0.4],
    [hap, yes, -1.9, 0.4], [notShown, yes, -0.45, 1], [shown, noY, 0.02, 1]] as const) {
    const c = read(m, y, x, cf); const row = sexRow(c, m); const say = sexSays('ambiguous', m, c)
    if (row) assert.ok(say.line.startsWith(row.head) && say.word === row.word && say.badge === row.badge, say.line)
    if (say.badge && !row) assert.equal(say.badge, c.constitution)
  }
}

console.log('sexMaterial.check.ts OK: no sex where two copies are not shown, and a chrX row with no parent')
