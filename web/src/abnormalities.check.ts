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
  runsOfHomozygosity, groupUnits, overlaps, unitsCarrying, mergeLoh, LOH_SEGMENT_MARKERS,
  callTriploidyOrigin,
  type Finding,
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
  // A WINDOW BIG ENOUGH TO BE MEASURED. At 1,000 markers and a background of 0.17 this fixture
  // expected 170 heterozygotes, and the false-positive sweep on real diploid arrays puts the
  // zero-crossing at 1,800 to 2,100 markers for exactly that background. The window is now 2,400,
  // which expects 408, just past the measured floor. The rates are unchanged, so the test still
  // asks what it always asked.
  const w = (het: number, logR?: number): WindowStat =>
    ({ chrom: '7', startBp: 1e6, endBp: 13e6, called: 2800, het, logR })
  // Background is taken from the windows themselves, so give it a normal population to sit in.
  const normal = Array.from({ length: 8 }, () => w(476))

  // Heterozygosity gone, copy number unchanged: copy-neutral.
  const flat = detectLoh([...normal, w(56, 0.01)])
  assert.equal(flat.length, 1, 'a depleted window with unchanged intensity is copy-neutral LOH')
  assert.equal(flat[0].cls, 'cnn-loh')
  assert.ok(flat[0].evidence.includes('copy number unchanged'))

  // Heterozygosity gone AND intensity moved: a deletion, which this detector must leave alone.
  const dropped = detectLoh([...normal, w(56, -0.9)])
  assert.equal(dropped.length, 0,
    'a depleted window whose intensity has moved is a DELETION. Calling it copy-neutral would name '
    + 'the wrong event and then invert its origin, since loss and gain read opposite signs')

  // A WHOLE-CHROMOSOME WINDOW MUST SAY SO, because it decides which detection floor applies and
  // the two differ by everything. A 12 Mb-scale interval on amplified material has NO floor at any
  // mosaic fraction with one parent, while a whole chromosome on the same material has 0.399.
  // Every copy-neutral finding is produced from one window per chromosome, so reporting them as
  // segments sent all of them to the segment floor and came back not-evaluable: the class was
  // detected and then refused an origin for a reason that did not apply to it.
  const whole = detectLoh([...normal.map((x) => ({ ...x, wholeChromosome: true })),
    { ...w(56, 0.01), wholeChromosome: true }])
  assert.equal(whole.length, 1)
  assert.equal(whole[0].wholeChromosome, true,
    'a window spanning the chromosome must be reported as a whole chromosome, or it is scored '
    + 'against a floor that does not exist for it')
  assert.equal(flat[0].wholeChromosome, false, 'and a sub-chromosomal one must not claim to be')

  // No intensity at all: reported, but explicitly not established as copy-neutral.
  const blind = detectLoh([...normal, w(56)])
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
  assert.equal(detectLoh(Array.from({ length: 8 }, () => w(56))).length, 0)
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
  // The PLOIDY needs no parent; WHOSE the extra set is does, and having one it is answerable. See
  // section 13: this class was on the unreachable list on reasoning that was true of band
  // structure alone.
  assert.ok(!ORIGIN_UNREACHABLE.has('triploidy'),
    'with a genotyped parent the extra set is nameable, so triploidy is not origin-unreachable')

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
  // The DOSAGE channel refuses here and must: it is self-referenced and there is nothing left to
  // reference against. The obligate-het channel does not, because it is Mendelian. Refusing every
  // channel on the strength of one channel's limit is what section 13 corrects.
  assert.ok(!ORIGIN_UNREACHABLE.has('complex'),
    'a complex genome keeps the Mendelian channel even though it loses the self-referenced one')

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

// --- 12. BOTH SCALES ARE SCANNED, AND THE REDUNDANCY IS COLLAPSED --------------------------------
//
// Copy-neutral LOH is scanned at whole-chromosome scale AND in sliding windows. The chromosome pass
// is the one with a detection floor on amplified material, since a 12 Mb-scale interval has none at
// any mosaic fraction with one parent; the windows are what give a partial event its own extent
// instead of reporting the whole chromosome. Scanning both leaves duplicates, so the widest
// interval covering a position wins and a segment that is only its chromosome restated is dropped.
{
  const f = (startBp: number, endBp: number, whole = false): Finding => ({
    cls: 'cnn-loh', chrom: '7', startBp, endBp, wholeChromosome: whole, evidence: 'x',
  })

  // A whole chromosome plus the windows inside it collapses to the chromosome alone.
  const nested = mergeLoh([f(0, 159e6, true), f(10e6, 30e6), f(20e6, 40e6)])
  assert.equal(nested.length, 1, `nested windows must collapse, got ${nested.length}`)
  assert.equal(nested[0].wholeChromosome, true,
    'and the survivor is the WIDER one, because it is the one that can be scored: a precise extent '
    + 'that comes back not-evaluable is worth less than a coarse one carrying a parent')

  // Half-overlapping windows over one real event become one interval at its full extent.
  const joined = mergeLoh([f(10e6, 30e6), f(20e6, 40e6), f(30e6, 50e6)])
  assert.equal(joined.length, 1, 'overlapping windows are one event seen several times')
  assert.equal(joined[0].startBp, 10e6)
  assert.equal(joined[0].endBp, 50e6, 'joined to its full extent rather than truncated')

  // Genuinely separate events on one chromosome stay separate.
  const apart = mergeLoh([f(10e6, 20e6), f(90e6, 110e6)])
  assert.equal(apart.length, 2, 'a gap between events must not be bridged')

  // And chromosomes never merge into each other.
  assert.equal(mergeLoh([f(10e6, 20e6), { ...f(10e6, 20e6), chrom: '8' }]).length, 2)
  assert.ok(LOH_SEGMENT_MARKERS >= 200, 'a window needs enough markers for a rate to mean anything')
}

// --- 13. THE TWO CLASSES DECLARED UNREACHABLE, AND WERE NOT ---------------------------------------
//
// Both were ruled out on reasoning that was true of ONE channel and applied to every channel. A
// triploid's extra set is named from allele fraction at the loaded parent's homozygous markers; a
// complex genome refuses the dosage channel and keeps the Mendelian one.
{
  const at = (share: number, n: number, g: 'AA' | 'BB' = 'AA') =>
    Array.from({ length: n }, () => [g, g === 'AA' ? 1 - share : share] as const)

  // Two of three copies are the loaded parent's: the extra set is theirs.
  const theirs = callTriploidyOrigin(at(2 / 3, 600))
  assert.equal(theirs.origin, 'extra-set-loaded-parent')
  assert.ok(theirs.why.includes('two thirds'))

  // One of three: the extra set is the other parent's.
  const others = callTriploidyOrigin(at(1 / 3, 600))
  assert.equal(others.origin, 'extra-set-other-parent')

  // BOTH HOMOZYGOTES MUST ORIENT THE SAME WAY, or the statistic acquires a directional null of the
  // kind this project has already had to remove once.
  assert.equal(callTriploidyOrigin(at(2 / 3, 600, 'BB')).origin, 'extra-set-loaded-parent')
  assert.equal(callTriploidyOrigin(at(1 / 3, 600, 'BB')).origin, 'extra-set-other-parent')

  // Occupancy, not per-marker calling: wide scatter must not defeat it, because per-marker band
  // assignment needs a BAF spread under 0.053 and fails on every amplified class.
  const noisy = Array.from({ length: 600 }, (_, i) =>
    ['AA', 1 - (2 / 3 + 0.09 * Math.sin(i * 2.399963))] as const)
  assert.equal(callTriploidyOrigin(noisy).origin, 'extra-set-loaded-parent',
    'per-marker dispersion must not defeat an occupancy statistic')

  // And it declines rather than guessing when there is not enough to count.
  assert.equal(callTriploidyOrigin(at(2 / 3, 50)).origin, 'unresolved')

  // NEITHER CLASS IS ON THE UNREACHABLE LIST ANY MORE, and the three real limits still are.
  assert.ok(!ORIGIN_UNREACHABLE.has('triploidy'))
  assert.ok(!ORIGIN_UNREACHABLE.has('complex'))
  for (const c of ['heterodisomy', 'reverse-segregation', 'tandem-vs-inserted'] as const) {
    assert.ok(ORIGIN_UNREACHABLE.has(c), `${c} is a real limit and must stay`)
  }
  // The taxonomy text must say HOW, not merely that it is possible.
  assert.ok(taxonomyFor('triploidy')!.origin.includes('two thirds'))
  assert.ok(taxonomyFor('complex')!.origin.includes('Mendelian rather than self-referenced'))
}

console.log('abnormalities.check.ts: all assertions passed, including every class enumerated with '
  + 'its limits, copy-neutral LOH declining to fire on a deletion, and no-run-found stating that '
  + 'it is not no-disomy')

// --- COPY-NEUTRAL LOH MUST NOT FIRE ON A GENOME THAT HAS NO HETEROZYGOSITY --------------------
//
// The detector measures depletion RELATIVE to the array's own mean heterozygosity. As that mean
// approaches zero the ratio is dominated by counting noise and any threshold is cleared by nothing
// at all. Measured on this code before the guards: an event-free genome of 800 windows at a true
// heterozygosity of 0.001 returned 39 copy-neutral findings. An external review's simulation of the
// same function returned more still, and a real run on uniparental material returned 183.
{
  let seed = 20260819
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  const eventFree = (het: number) => {
    const w = []
    for (let i = 0; i < 800; i += 1) {
      let h = 0
      for (let k = 0; k < 2000; k += 1) if (rnd() < het) h += 1
      w.push({ chrom: String((i % 22) + 1), startBp: i * 1e6, endBp: i * 1e6 + 1e6,
        called: 2000, het: h, logR: 0, wholeChromosome: false })
    }
    return w
  }

  for (const het of [0.001, 0.002, 0.005, 0.010, 0.030]) {
    const found = detectLoh(eventFree(het))
    assert.equal(found.length, 0,
      `an event-free genome at heterozygosity ${het} returned ${found.length} copy-neutral `
      + 'findings, which is the artefact these guards exist to stop')
  }

  // The decisive guard where the caller knows the zygosity: a genome with one parental
  // contribution is homozygous by construction, so the event is not defined in it at any numbers.
  assert.equal(detectLoh(eventFree(0.30), { zygosity: 'uniparental_homozygous' }).length, 0,
    'a uniparental genome must produce no copy-neutral findings whatever its windows look like')

  // AND A REAL EVENT MUST STILL BE FOUND, or the guards have simply disabled the detector.
  const planted = eventFree(0.30)
  for (let i = 100; i < 160; i += 1) planted[i].het = Math.round(planted[i].called * 0.30 * 0.2)
  assert.equal(detectLoh(planted).length, 60,
    'a genuine 80% depletion on a clean genome must still be detected on every planted window')
  console.log('  copy-neutral: 0 false windows on event-free genomes at every heterozygosity '
    + 'tested, uniparental excluded outright, 60 of 60 planted events still found')
}
