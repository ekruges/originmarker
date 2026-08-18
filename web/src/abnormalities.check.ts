// Self-check for the abnormality taxonomy. Run: node src/abnormalities.check.ts
//
// Three properties carry this module.
//
//   A CLASS THAT CANNOT BE ANSWERED MUST STILL APPEAR. The reason this file exists is that silence
//   read as absence: a report listing what was found said nothing about what was never checked.
//   An impossible class that quietly stops being listed is the original defect returning.
//
//   THE COPY-NEUTRAL DETECTOR MUST NOT FIRE ON A DELETION. Depleted heterozygosity is equally the
//   signature of a lost copy, and calling that copy-neutral would name the wrong event and then
//   invert its origin, since loss and gain read opposite signs.
//
//   NO RUN FOUND IS NOT NO DISOMY. Heterodisomy is invisible to a run-based detector and 32% of
//   confirmed disomy carries no run at all, so the negative must say so in its own words.
import assert from 'node:assert/strict'
import {
  TAXONOMY, taxonomyFor, detectLoh, detectUpd, detectTriploidy, detectComplex, callUniformity,
  unanswerable, ORIGIN_UNREACHABLE, LCSH_REPORT_MB, UPD_NO_STRETCH_RATE, LOH_DEPLETION,
  CLUSTERED_DROPOUT_DEPLETION, COMPLEX_DEVIANT_FRACTION, TRIPLOID_BANDS, DIPLOID_BAND,
  runsOfHomozygosity, groupUnits, overlaps, unitsCarrying,
  type WindowStat, type RunOfHomozygosity,
} from './abnormalities.ts'

// --- 1. EVERY CLASS IS ENUMERATED, INCLUDING THE ONES THAT CANNOT BE ANSWERED ---------------------
{
  const want = ['monosomy', 'trisomy', 'segmental-deletion', 'segmental-duplication', 'cnn-loh',
    'isodisomy', 'heterodisomy', 'segmental-upd', 'triploidy', 'haploidy', 'complex',
    'gamete-de-novo', 'reverse-segregation', 'tandem-vs-inserted']
  for (const c of want) {
    const t = taxonomyFor(c as never)
    assert.ok(t, `${c} is missing from the taxonomy`)
    assert.ok(t.by.length > 10, `${c} must say what it rests on`)
    assert.ok(t.origin.length > 10, `${c} must say what its origin rests on, or that none exists`)
  }
  // The impossible ones must carry a stated LIMIT, or a reader cannot tell a hard wall from a
  // floor that better data would clear.
  for (const c of ['heterodisomy', 'reverse-segregation', 'tandem-vs-inserted'] as const) {
    const t = taxonomyFor(c)!
    assert.equal(t.detectable, 'no')
    assert.ok(t.limit && t.limit.length > 20, `${c} must state why it is a hard limit`)
  }
  assert.equal(TAXONOMY.length, want.length, 'the taxonomy and this list must not drift apart')
}

// --- 2. SILENCE IS NEVER ABSENCE ------------------------------------------------------------------
//
// The unanswerable set must SHRINK as the user supplies more, which is what makes it actionable
// rather than a disclaimer: heterodisomy becomes reachable with both parents, and the timing of a
// segmental change with a second unit of the same embryo.
{
  const oneParentOneUnit = unanswerable(false, 1).map((t) => t.cls)
  assert.ok(oneParentOneUnit.includes('heterodisomy'))
  assert.ok(oneParentOneUnit.includes('gamete-de-novo'))
  assert.ok(oneParentOneUnit.includes('reverse-segregation'))

  const twoParents = unanswerable(true, 1).map((t) => t.cls)
  assert.ok(!twoParents.includes('heterodisomy'), 'a second parent must lift the heterodisomy wall')
  assert.ok(twoParents.includes('reverse-segregation'), 'but not the ones that are truly hard')

  const twoUnits = unanswerable(false, 2).map((t) => t.cls)
  assert.ok(!twoUnits.includes('gamete-de-novo'), 'a second unit must lift the timing wall')

  // The three hard limits survive everything.
  const best = unanswerable(true, 4).map((t) => t.cls)
  assert.deepEqual(best.sort(), ['reverse-segregation', 'tandem-vs-inserted'].sort(),
    'exactly the platform-level limits remain when the user has supplied everything they can')
}

// --- 3. COPY-NEUTRAL LOH DOES NOT FIRE ON A DELETION -----------------------------------------------
{
  const w = (het: number, logR?: number): WindowStat =>
    ({ chrom: '7', startBp: 1e6, endBp: 13e6, called: 1000, het, logR })
  // Background is taken from the windows themselves, so give it a normal population to sit in.
  const normal = Array.from({ length: 8 }, () => w(170))

  // Heterozygosity gone, copy number unchanged: copy-neutral.
  const flat = detectLoh([...normal, w(20, 0.01)])
  assert.equal(flat.length, 1, 'a depleted window with unchanged intensity is copy-neutral LOH')
  assert.equal(flat[0].cls, 'cnn-loh')
  assert.ok(flat[0].evidence.includes('copy number unchanged'))

  // Heterozygosity gone AND intensity moved: a deletion, which this detector must leave alone.
  const dropped = detectLoh([...normal, w(20, -0.9)])
  assert.equal(dropped.length, 0,
    'a depleted window whose intensity has moved is a DELETION. Calling it copy-neutral would name '
    + 'the wrong event and then invert its origin, since loss and gain read opposite signs')

  // No intensity at all: reported, but explicitly not established as copy-neutral.
  const blind = detectLoh([...normal, w(20)])
  assert.equal(blind.length, 1)
  assert.ok(blind[0].evidence.includes('loss-or-copy-neutral'),
    'without intensity the class must be left open rather than assumed')

  // The false-positive floor must travel with every call.
  assert.ok(flat[0].flag!.includes('40%'))
  assert.ok(LOH_DEPLETION > CLUSTERED_DROPOUT_DEPLETION,
    'the calling threshold must sit ABOVE the depletion that dropout alone already produces on 1% '
    + 'of event-free windows, or the detector calls one clean window in a hundred')

  // A globally low-heterozygosity sample must not be called end to end, because the background is
  // the array's own.
  assert.equal(detectLoh(Array.from({ length: 8 }, () => w(20))).length, 0)
}

// --- 4. NO RUN FOUND IS NOT NO DISOMY --------------------------------------------------------------
{
  const run = (chrom: string, mb: number, whole = false): RunOfHomozygosity =>
    ({ chrom, startBp: 0, endBp: mb * 1e6, markers: Math.round(mb * 300), wholeChromosome: whole })

  assert.equal(detectUpd([run('7', 5)]).length, 0, `${LCSH_REPORT_MB} Mb is the reporting floor`)

  const one = detectUpd([run('7', 40, true)])
  assert.equal(one.length, 1)
  assert.equal(one[0].cls, 'isodisomy', 'a whole-chromosome run is isodisomy')
  assert.ok(one[0].evidence.includes('Heterodisomy leaves heterozygosity normal'),
    'the negative must state what it cannot see')
  assert.ok(one[0].evidence.includes('32%'),
    'and the share of real disomy that carries no run at all must travel with it')
  assert.equal(one[0].flag, undefined, 'a single chromosome is not consanguinity')

  assert.equal(detectUpd([run('7', 30)])[0].cls, 'segmental-upd',
    'a run short of the whole chromosome is segmental')

  // Runs on many chromosomes at once are shared ancestry, and the calls are FLAGGED rather than
  // withdrawn: the runs are real, it is what they imply that changes.
  const many = detectUpd(['1', '4', '9', '15', '20'].map((c) => run(c, 30)))
  assert.equal(many.length, 5)
  for (const f of many) {
    assert.ok(f.flag && f.flag.includes('shared parental ancestry'),
      'runs across many chromosomes must flag consanguinity on every one of them')
  }
  assert.ok(UPD_NO_STRETCH_RATE > 0.3)
}

// --- 5. TRIPLOIDY NEEDS THE HALF BAND VACATED, NOT JUST THE THIRDS OCCUPIED ------------------------
{
  const at = (centre: number, n: number) => Array.from({ length: n }, () => centre)

  const tri = detectTriploidy([...at(TRIPLOID_BANDS[0], 400), ...at(TRIPLOID_BANDS[1], 400)])
  assert.ok(tri, 'mass at both thirds with the half band empty is a third set')
  assert.equal(tri!.cls, 'triploidy')
  // The ploidy needs no parent, but the PARENT of the extra set is a different question.
  assert.ok(tri!.originBlocked!.includes('Digynic and diandric'))
  assert.ok(ORIGIN_UNREACHABLE.has('triploidy'))

  // A diploid genome occupies the half band, and that occupancy alone must veto the call even if
  // the thirds are also populated by noise. This is the discriminating half of the statistic.
  const diploid = detectTriploidy([...at(DIPLOID_BAND, 500), ...at(TRIPLOID_BANDS[0], 400)])
  assert.equal(diploid, null, 'an occupied half band excludes a third set')

  assert.equal(detectTriploidy(at(TRIPLOID_BANDS[0], 50)), null, 'too few markers to attempt it')
}

// --- 6. A COMPLEX GENOME BLOCKS EVERY ORIGIN, AND IS NOT A QUALITY FAILURE -------------------------
{
  assert.equal(detectComplex(2, 22, 0.95), null, 'two deviant autosomes is not a chaotic genome')

  const chaotic = detectComplex(12, 22, 0.95)
  assert.ok(chaotic)
  assert.ok(chaotic!.originBlocked!.includes('no undisturbed part left'))
  assert.ok(chaotic!.originBlocked!.includes('rather than a fault of the array'),
    'the difference between a disturbed genome and a poor array is what a reader acts on')
  assert.ok(ORIGIN_UNREACHABLE.has('complex'))

  // The same conclusion arrives from a collapsed call rate, for the same reason.
  assert.ok(detectComplex(1, 22, 0.5))
  assert.ok(COMPLEX_DEVIANT_FRACTION < 0.5)
}

// --- 7. TIMING NEEDS A SECOND UNIT, AND SAYS SO WHEN IT HAS ONLY ONE -------------------------------
{
  const alone = callUniformity(1, 1)
  assert.equal(alone.mechanism, 'unresolved')
  assert.ok(alone.why.includes('second array of the same embryo'),
    'the remedy must be stated, because it is one the user may already be able to supply')

  assert.equal(callUniformity(3, 3).mechanism, 'meiotic')
  assert.equal(callUniformity(1, 3).mechanism, 'post-zygotic')
  assert.ok(callUniformity(3, 3).why.includes('64 of'))
}

// --- 8. THE RUN FINDER TOLERATES DROP-IN, BECAUSE REAL RUNS ARE NOT CLEAN ------------------------
//
// Allele drop-in puts a heterozygous call at 4.35% of truly homozygous markers on amplified
// material, measured over 113 same-genome pairs. A run-finder that broke on the first heterozygote
// would return nothing on exactly the material this tool exists for.
{
  const mk = (n: number, hetAt: number[] = []) =>
    Array.from({ length: n }, (_, i) => ({ chrom: '7', pos: i * 100_000, het: hetAt.includes(i) }))

  const clean = runsOfHomozygosity(mk(400))
  assert.equal(clean.length, 1, 'a clean stretch is one run')
  assert.equal(clean[0].markers, 400)

  // Scattered drop-in must not shatter the run.
  const noisy = runsOfHomozygosity(mk(400, [50, 137, 240, 361]))
  assert.equal(noisy.length, 1, `drop-in must not break a real run, got ${noisy.length} pieces`)
  assert.equal(noisy[0].markers, 400)

  // A genuinely heterozygous stretch must still break it.
  const broken = runsOfHomozygosity(mk(400, Array.from({ length: 120 }, (_, i) => 200 + i)))
  assert.ok(broken.every((r) => r.markers < 400),
    'a truly heterozygous region must end the run rather than being tolerated')

  // Runs under the marker floor are not returned at all.
  assert.equal(runsOfHomozygosity(mk(40)).length, 0)

  // And the whole-chromosome flag is what separates isodisomy from a long segmental run.
  const whole = runsOfHomozygosity(mk(400), { chromEndBp: new Map([['7', 40_000_000]]) })
  assert.equal(whole[0].wholeChromosome, true, '400 markers spanning 39.9 Mb of a 40 Mb chromosome')
  const part = runsOfHomozygosity(mk(400), { chromEndBp: new Map([['7', 159_000_000]]) })
  assert.equal(part[0].wholeChromosome, false, 'the same run on a full-length chromosome is segmental')
}

// --- 9. UNITS OF ONE EMBRYO ARE MEASURED, NOT DECLARED --------------------------------------------
//
// Two biopsies of one embryo are the same genome and concord like replicate arrays, 95.8% against
// 54.9% for a parent-offspring pair. Asking the user to declare the grouping would have been
// easier and worse: a mislabelled group yields a confidently wrong MECHANISM, and the separation
// here is 41 points wide.
{
  // Concordance matrix: A and B are one embryo, C is a sibling, D is unrelated.
  const conc: Record<string, Record<string, number>> = {
    A: { A: 1, B: 0.958, C: 0.549, D: 0.31 },
    B: { A: 0.958, B: 1, C: 0.551, D: 0.30 },
    C: { A: 0.549, B: 0.551, C: 1, D: 0.29 },
    D: { A: 0.31, B: 0.30, C: 0.29, D: 1 },
  }
  const g = groupUnits(['A', 'B', 'C', 'D'], (a, b) => conc[a][b])
  assert.equal(g.length, 3, `three embryos, got ${g.length}`)
  assert.deepEqual(g[0], ['A', 'B'], 'the two units of one embryo group together')
  assert.deepEqual(g[1], ['C'])
  assert.deepEqual(g[2], ['D'])

  // A chain of near-misses must not merge two embryos through an intermediate, which is why
  // membership requires agreement with EVERY member rather than any one of them.
  const chain: Record<string, Record<string, number>> = {
    X: { X: 1, Y: 0.95, Z: 0.60 },
    Y: { X: 0.95, Y: 1, Z: 0.95 },
    Z: { X: 0.60, Y: 0.95, Z: 1 },
  }
  const c = groupUnits(['X', 'Y', 'Z'], (a, b) => chain[a][b])
  assert.ok(c.some((grp) => grp.includes('Z') && !grp.includes('X')),
    'X and Z disagree, so they must not end up in one embryo via Y')

  // An unmeasurable pair (too few shared markers) must not group.
  assert.equal(groupUnits(['P', 'Q'], () => NaN).length, 2)
}

// --- 10. UNIFORMITY USES OVERLAP, NOT IDENTICAL EDGES ----------------------------------------------
//
// Two biopsies of one embryo do not place a breakpoint identically. Requiring identical edges would
// report every genuinely uniform event as non-uniform, which inverts the mechanism call.
{
  const ev = (chrom: string, a: number, b: number) => ({ chrom, startBp: a, endBp: b })
  assert.ok(overlaps(ev('7', 1e6, 20e6), ev('7', 19e6, 40e6)))
  assert.ok(!overlaps(ev('7', 1e6, 10e6), ev('7', 11e6, 20e6)))
  assert.ok(!overlaps(ev('7', 1e6, 20e6), ev('8', 1e6, 20e6)), 'different chromosomes never overlap')

  const target = ev('7', 5e6, 25e6)
  // Present in all three units, with edges that differ by megabases: uniform.
  assert.equal(unitsCarrying(target, [
    [ev('7', 4.8e6, 25.2e6)], [ev('7', 5.4e6, 24.1e6)], [ev('7', 6e6, 30e6)],
  ]), 3)
  assert.equal(callUniformity(3, 3).mechanism, 'meiotic')

  // Present in one of three: confined to a lineage.
  assert.equal(unitsCarrying(target, [[ev('7', 5e6, 25e6)], [], [ev('9', 1e6, 9e6)]]), 1)
  assert.equal(callUniformity(1, 3).mechanism, 'post-zygotic')
}

console.log('abnormalities.check.ts: all assertions passed, including every class enumerated with '
  + 'its limits, copy-neutral LOH declining to fire on a deletion, and no-run-found stating that '
  + 'it is not no-disomy')
