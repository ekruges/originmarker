// Self-check for the operator's stage declaration. Run: node src/declaredStage.check.ts
//
// What is tested is the thing that decides whether a reader is warned: that a declaration is USED,
// that the inference still runs underneath it, and that a disagreement about PLOIDY is reported as
// such rather than lumped in with a disagreement between two diploid stages. Getting the last one
// wrong lets a haploid array be scored as biparental with no notice, which is a confident wrong
// parental call on material that has only one parent in it.
import assert from 'node:assert/strict'
import { reconcileStage, ploidyOf } from './declaredStage.ts'
import { stageDefaults, type Stage, type StageCall } from './stage.ts'

const inferred = (stage: Stage): StageCall => ({
  ...stageDefaults(stage), why: `the array reads ${stage}`,
})

// --- 1. no declaration changes nothing --------------------------------------------------------
{
  const a = reconcileStage(inferred('blastomere'), undefined, stageDefaults)
  assert.equal(a.used.stage, 'blastomere')
  assert.equal(a.agrees, true)
  assert.equal(a.ploidyConflict, false)
  assert.equal(a.notice, '', 'silence when there is nothing to reconcile')
  assert.equal(a.declared, undefined)
}

// --- 2. a declaration WINS, and the inference is still reported ---------------------------------
//
// The operator dissected the cell. The tool guessed from heterozygosity. Overriding silently would
// be worse than either, so both must survive into the result.
{
  const a = reconcileStage(inferred('haploid'), 'blastomere', stageDefaults)
  assert.equal(a.used.stage, 'blastomere', 'the declaration is what the run uses')
  assert.equal(a.inferred.stage, 'haploid', 'and the array reading survives beside it')
  assert.equal(a.agrees, false)
  assert.ok(a.notice.includes('blastomere') && a.notice.includes('haploid'),
    'the notice must name both, or a reader cannot tell what disagreed')
  // Dropout must come from the DECLARED stage: it parameterises every likelihood below, and the
  // operator has just said this is not the material the inference named.
  assert.equal(a.used.dropout, stageDefaults('blastomere').dropout)
  assert.notEqual(a.used.dropout, stageDefaults('haploid').dropout)
}

// --- 3. PLOIDY DISAGREEMENT IS ITS OWN CATEGORY -------------------------------------------------
//
// The case Joy named: the operator says diploid, the array says haploid. Two diploid stages
// disagreeing moves a detection floor. Disagreeing about ploidy moves whether parent of origin is
// a question that can be asked at all, so it must not read the same.
{
  const ploidyRow = reconcileStage(inferred('haploid'), 'blastomere', stageDefaults)
  assert.equal(ploidyRow.ploidyConflict, true)
  assert.ok(/PLOIDY/.test(ploidyRow.notice), 'a ploidy conflict must say so in the notice')
  // Both readings of the conflict are given, because which one is wrong changes what to do.
  assert.ok(/dropout/i.test(ploidyRow.notice), 'says what it means if the declaration is right')
  assert.ok(/not what the tube says|no comparison/i.test(ploidyRow.notice),
    'and what it means if the array is right')

  // Two diploid stages: a real disagreement, but not this one.
  const floorRow = reconcileStage(inferred('single-cell'), 'trophectoderm', stageDefaults)
  assert.equal(floorRow.agrees, false)
  assert.equal(floorRow.ploidyConflict, false)
  assert.ok(!/PLOIDY/.test(floorRow.notice))
  assert.ok(/detection floor/.test(floorRow.notice),
    'a same-ploidy disagreement must say what it actually changes')

  // And the mirror: operator says haploid, array reads diploid.
  const mirror = reconcileStage(inferred('bulk'), 'haploid', stageDefaults)
  assert.equal(mirror.ploidyConflict, true)
  assert.ok(/drop-in|second cell/i.test(mirror.notice),
    'the mirror case has a different explanation and must not reuse the other one')
}

// --- 4. every stage has a ploidy, and it is the right one ---------------------------------------
//
// An amplified single cell of a diploid line is DIPLOID however much dropout it carries. Calling it
// haploid here would silence the parental channels on exactly the material they work on.
{
  assert.equal(ploidyOf('bulk'), 'diploid')
  assert.equal(ploidyOf('trophectoderm'), 'diploid')
  assert.equal(ploidyOf('blastomere'), 'diploid')
  assert.equal(ploidyOf('single-cell'), 'diploid',
    'an amplified single cell of a diploid line is diploid, whatever its dropout')
  assert.equal(ploidyOf('haploid'), 'haploid')
  assert.equal(ploidyOf('failed'), 'unknown')
  assert.equal(ploidyOf('unknown'), 'unknown')
  // 'failed' and 'unknown' name the absence of a call, so they must never raise a ploidy conflict:
  // there is nothing to conflict with.
  assert.equal(reconcileStage(inferred('failed'), 'blastomere', stageDefaults).ploidyConflict, false)
  assert.equal(reconcileStage(inferred('bulk'), 'unknown', stageDefaults).ploidyConflict, false)
}

// --- 5. declaring what the array already says is agreement, not an override ---------------------
{
  const a = reconcileStage(inferred('trophectoderm'), 'trophectoderm', stageDefaults)
  assert.equal(a.agrees, true)
  assert.equal(a.ploidyConflict, false)
  assert.ok(a.notice.includes('the same'), 'confirmation is worth saying, quietly')
}

// --- 6. no em dashes in anything an operator reads ----------------------------------------------
{
  for (const [inf, dec] of [['haploid', 'blastomere'], ['bulk', 'haploid'],
    ['single-cell', 'trophectoderm'], ['trophectoderm', 'trophectoderm']] as [Stage, Stage][]) {
    const a = reconcileStage(inferred(inf), dec, stageDefaults)
    assert.ok(!a.notice.includes('—'), `em dash in the notice for ${inf} vs ${dec}`)
    assert.ok(!a.used.why.includes('—'))
  }
}

console.log('declaredStage.check.ts: the declaration wins, the inference survives beside it, '
  + 'and a ploidy disagreement reads differently from a floor disagreement')
