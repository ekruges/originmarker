// WHAT A RUN MUST NOT REPORT. Run: node --experimental-strip-types src/scoreSample.check.ts
//
// Three ways this file has said one thing and reported another, each pinned here on a synthetic
// array so a reader can see the event that was built and the row that came out:
//
//   1. AN ARRAY THE TOOL REJECTED still shipped calls. The guard emptied the taxonomy findings
//      alone, so whole-chromosome aneuploidy and copy-number segments walked past it, and the
//      Mendelian channel then named a parent for each of them. Measured on GSM4472448, a real
//      trophectoderm array carrying a constructed one-copy chromosome, with the stage forced to
//      failed: 1 aneuploidy, 3 segments and 4 origin calls naming a parent, three of those at
//      posterior 0.793 to 0.794 with the array already declared unreadable.
//   2. A ONE-COPY CHROMOSOME reported a loss and an isodisomy side by side. Isodisomy asserts
//      BOTH copies came from one parent, on an event that has one copy.
//   3. A DETECTED SEGMENT got no origin output of any kind wherever the Mendelian channel was
//      off, which is every uniparental sample.
//
// The fixture is built here rather than read, so this runs with no data to supply. It is the
// shape of the arrays this tool is for and not a substitute for them: the measurements behind
// each of the three sit in the audit, on real material.
import assert from 'node:assert/strict'
import { headerMap, parseRow, accumulate, accumulateBaf, emptyBafSums, finishProfile } from './ingest.ts'
import { emptyParent, collectParentRow, finishParent, emptyCollected, collectRow, scoreSample } from './scoreSample.ts'
import { mendelParent } from './defects.ts'

/** A deterministic stream, so a failure is the code changing and never the draw changing. */
function rng(seed: number): () => number {
  let x = seed >>> 0
  return () => {
    x ^= x << 13; x >>>= 0
    x ^= x >> 17
    x ^= x << 5; x >>>= 0
    return x / 0x100000000
  }
}

const HEAD = 'probeset_id\tchr\tposition\tlog2R\tbaf\tcopy_number\tgenotype\tBestProbeset'
const CHROMS = Array.from({ length: 22 }, (_, i) => String(i + 1))
/** The segment chromosome carries more, because a segment is not reported under 2,400 markers,
 *  and it must stay a SEGMENT: a chromosome whose call rate collapses below 0.60 of the array's
 *  median is a whole-chromosome event and is not scanned for a segment inside it. */
const PER_CHROM = (c: string) => (c === '9' ? 9_000 : 2_400)
const PARENT_HET = 0.170
/** How often the other parent transmits the same allele. Two people share the common allele at
 *  most markers of a panel like this; independent draws put the sample at 47% heterozygous, which
 *  no genome reaches and which the stage inference rejects outright. */
const SHARED = 0.86

const ONE_COPY = '5'
const ISODISOMY = '7'
const SEGMENT = '9'
const ISO_TRISOMY = '11'
const MAT_LOST = '13'

/**
 * One bulk parent and one biparental sample carrying three events.
 *
 *   chr5  ONE COPY. The loaded parent's copy is gone, so every call is the allele that parent does
 *         not have, which is what hemizygosity reads as, and the intensity is down by one copy.
 *   chr7  TRUE ISODISOMY. Both copies from the loaded parent: heterozygosity is gone and the
 *         intensity is UNTOUCHED, because copy number is still two.
 *   chr9  A SEGMENT. 3,000 consecutive markers not called, with the intensity agreeing.
 *   chr13 THE MATERNAL COPY GONE, the mirror of chr5. One copy, and every call is an allele the
 *         MOTHER cannot have supplied. This is the direction test: chr5 and chr13 differ only in
 *         which parent's copy was removed, so a rule that inverted the two would name both wrong
 *         while every other assertion in this file still passed.
 *   chr11 AN ISODISOMIC TRISOMY. Three copies AND homozygous end to end. The two facts are
 *         separate: the gain is one, and both copies coming from one parent is the other, which is
 *         the trisomy-rescue and imprinting mechanism. Suppressing runs of homozygosity on every
 *         aneuploid chromosome rather than only on the one-copy ones deleted the second.
 */
function synthesise(): { parent: string; mother: string; sample: string } {
  const r = rng(20260826)
  const parent = [HEAD]
  const mother = [HEAD]
  const sample = [HEAD]
  for (const chrom of CHROMS) {
    const n = PER_CHROM(chrom)
    // PER-CHROMOSOME INTENSITY WOBBLE, because the whole-chromosome test is read against the
    // spread of this array's own chromosome medians. A fixture whose chromosomes all sit at
    // exactly zero has a null of width zero, and then any offset at all reads as an event: the
    // segment chromosome was reported as a whole-chromosome loss purely because its median was
    // dragged a hundredth of a log2 unit by the segment inside it.
    const chromOffset = r() * 0.06 - 0.03
    for (let i = 0; i < n; i += 1) {
      const id = `AX-${chrom}-${i}`
      const pos = 100_000 + i * 90_000
      const pg = r() < PARENT_HET ? 1 : (r() < 0.5 ? 0 : 2)
      parent.push([id, chrom, pos, (r() * 0.1 - 0.05).toFixed(4),
        pg === 1 ? (0.5 + (r() - 0.5) * 0.06).toFixed(4) : (pg === 0 ? 0.02 : 0.98).toFixed(4),
        '2.0', String(pg), '1'].join('\t'))

      const fromLoaded = pg === 1 ? (r() < 0.5 ? 0 : 1) : (pg === 0 ? 0 : 1)
      const fromOther = r() < SHARED ? fromLoaded : 1 - fromLoaded
      // THE MOTHER, built to be consistent with what she transmitted rather than drawn beside it.
      // Heterozygous at the same rate as the father, and otherwise homozygous for the allele the
      // sample received from her, so no marker asks her to transmit an allele she does not carry.
      const mg = r() < PARENT_HET ? 1 : 2 * fromOther
      mother.push([id, chrom, pos, (r() * 0.1 - 0.05).toFixed(4),
        mg === 1 ? (0.5 + (r() - 0.5) * 0.06).toFixed(4) : (mg === 0 ? 0.02 : 0.98).toFixed(4),
        '2.0', String(mg), '1'].join('\t'))
      let sg = fromLoaded + fromOther
      let lrr = chromOffset + r() * 0.1 - 0.05
      let called = r() < 0.98
      if (chrom === ONE_COPY) {
        sg = 2 * fromOther
        lrr -= 0.9
      } else if (chrom === ISODISOMY) {
        sg = 2 * fromLoaded
      } else if (chrom === SEGMENT && i >= 3_000 && i < 5_600) {
        called = false
        lrr -= 1.5
      } else if (chrom === ISO_TRISOMY) {
        sg = 2 * fromLoaded
        lrr += Math.log2(3 / 2)
      } else if (chrom === MAT_LOST) {
        sg = 2 * fromLoaded
        lrr -= 0.9
      }
      sample.push([id, chrom, pos, lrr.toFixed(4),
        sg === 1 ? 0.5 + (r() - 0.5) * 0.08 : sg === 0 ? 0.02 + r() * 0.02 : 0.96 + r() * 0.02,
        '2.0', called ? String(sg) : '-1', '1'].join('\t'))
    }
  }
  return {
    parent: parent.join('\n'), mother: mother.join('\n'), sample: sample.join('\n'),
  }
}

function rowsOf(text: string) {
  const lines = text.split('\n')
  const map = headerMap(lines[0])
  assert.ok(map, 'the fixture header must parse through the shipped ingest')
  const out = []
  for (let i = 1; i < lines.length; i += 1) {
    const row = parseRow(lines[i], map)
    if (row) out.push(row)
  }
  return out
}

const { parent, mother, sample } = synthesise()
const parentRows = rowsOf(parent)
const motherRows = rowsOf(mother)
const sampleRows = rowsOf(sample)

function indexOf(rows: ReturnType<typeof rowsOf>) {
  const pacc = emptyParent()
  const byChrom = new Map(); const bafSums = emptyBafSums(); let first = ''
  for (const row of rows) {
    if (!first) first = row.probesetId
    accumulate(row, byChrom); accumulateBaf(row, bafSums); collectParentRow(row, pacc)
  }
  return finishParent(pacc, finishProfile('p', byChrom, bafSums, first).build.build)
}

async function run(declaredStage?: 'failed', withMother = false) {
  const pat = indexOf(parentRows)
  const mat = withMother ? indexOf(motherRows) : null
  const acc = emptyCollected(pat, mat)
  const byChrom = new Map(); const bafSums = emptyBafSums(); let first = ''
  for (const row of sampleRows) {
    if (!first) first = row.probesetId
    accumulate(row, byChrom); accumulateBaf(row, bafSums); collectRow(row, pat, mat, acc)
  }
  const profile = finishProfile('sample', byChrom, bafSums, first)
  return scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: 'sample',
    log: () => {}, declaredStage,
  })
}

// ---------------------------------------------------------------- the array as it reads
{
  const r = await run()
  assert.notEqual(r.stage?.stage, 'failed',
    'the fixture must resolve to a stage, or it is exercising the rejection path below rather '
    + `than the detection path. Got ${r.stage?.why}`)

  const loss = (r.chroms ?? []).find((c) => c.chrom === ONE_COPY)
  assert.equal(loss?.aneuploidy, 'loss',
    `chr${ONE_COPY} was built at one copy, so it must be reported as a whole-chromosome loss. `
    + `Got ${loss?.aneuploidy}. Nothing below means anything if the event itself is not detected.`)

  const disomyHere = (r.findings ?? []).filter((f) => f.chrom === ONE_COPY
    && (f.cls === 'isodisomy' || f.cls === 'segmental-upd'))
  assert.deepEqual(disomyHere, [],
    `chr${ONE_COPY} has ONE copy and is reported as a loss, so a run of homozygosity on it is `
    + 'that same event described twice. Classified it becomes isodisomy, which asserts BOTH '
    + `copies came from one parent. Got ${disomyHere.map((f) => f.cls).join(', ')}.`)

  // AND A GAIN THAT IS HOMOZYGOUS END TO END KEEPS BOTH FACTS. The filter above exists for the
  // one-copy artefact; applied to gains it deleted an isodisomic trisomy, and the gain row it left
  // behind says nothing about both copies coming from one parent.
  const triGain = (r.chroms ?? []).find((c) => c.chrom === ISO_TRISOMY)
  assert.equal(triGain?.aneuploidy, 'gain',
    `chr${ISO_TRISOMY} was built at three copies and must be reported as a gain. `
    + `Got ${triGain?.aneuploidy}.`)
  const triUpd = (r.findings ?? []).filter((f) => f.chrom === ISO_TRISOMY
    && (f.cls === 'isodisomy' || f.cls === 'segmental-upd'))
  assert.ok(triUpd.length > 0,
    `chr${ISO_TRISOMY} has THREE copies and is homozygous end to end, which is an isodisomic `
    + 'trisomy: a separate fact from the gain and the mechanism behind trisomy rescue and '
    + `imprinting disorders. Got ${(r.findings ?? []).filter((f) => f.chrom === ISO_TRISOMY)
      .map((f) => f.cls).join(', ') || 'nothing'}.`)

  const real = (r.findings ?? []).filter((f) => f.chrom === ISODISOMY && f.cls === 'isodisomy')
  assert.equal(real.length, 1,
    `chr${ISODISOMY} carries TWO copies, both from the loaded parent, at normal intensity: that `
    + 'is isodisomy and it must still be reported as isodisomy. Got '
    + `${(r.findings ?? []).filter((f) => f.chrom === ISODISOMY).map((f) => f.cls).join(', ')}. `
    + 'Suppressing this one would be deleting a class rather than fixing a label.')

  const seg = (r.segments ?? []).filter((s) => s.chrom === SEGMENT)
  assert.ok(seg.length >= 1,
    `chr${SEGMENT} carries a 2,600-marker copy-loss segment and the scan must find it. Got none.`)
  const segRows = (r.dosageCalls ?? []).filter((c) => c.where.startsWith(`chr${SEGMENT} `))
  assert.ok(segRows.length >= 1,
    'a detected segment must leave the run with an origin row: a call, or a refusal that names '
    + 'the channel and the reason. It reached the Mendelian channel and nothing else, so on every '
    + 'sample where that channel is off it got silence. Got no row.')
}

// ---------------------------------------------------------------- the same array, rejected
//
// `integrity.ts` calls a failed array uninterpretable and writes the limit this asserts: "Nothing
// on this sample is reportable. Do not read its changes, its parental calls or its counts."
// Position does not survive here either, though it survives the other two uninterpretable
// verdicts: on a low call rate or a genome with no undisturbed remainder the array IS reading a
// genome and only the parent is withheld, where a failed array is the case with no genome being
// read, so there is no event to have a position.
{
  const r = await run('failed')
  assert.equal(r.integrity?.level, 'uninterpretable',
    'a failed array is uninterpretable, and this file must suppress on the same condition rather '
    + `than hold a second opinion beside it. Got ${r.integrity?.level}.`)
  const aneu = (r.chroms ?? []).filter((c) => c.aneuploidy).map((c) => c.chrom)
  assert.deepEqual(aneu, [],
    `an array that did not resolve to a stage must report no aneuploidy: it was found by the `
    + `measurement the tool has just declared broken. Got chr${aneu.join(', chr')}.`)
  assert.deepEqual(r.segments, [], 'and no segments, for the same reason')
  assert.deepEqual(r.findings, [], 'and no taxonomy findings, for the same reason')
  assert.deepEqual(r.dosageCalls ?? [], [], 'and no dosage origin rows')
  assert.deepEqual(r.oneParent ?? [], [],
    'and no Mendelian origin rows. This is the one that shipped a named parent at posterior '
    + '0.79 on an array the tool had already said it could not read.')
  assert.deepEqual(r.gains ?? [], [], 'and no gain rows')
  assert.deepEqual(r.losses ?? [], [], 'and no loss rows')
}


// ---------------------------------------------------------------- BOTH PARENTS, AND WHICH IS WHICH
//
// The Mendelian channel is right 0.8539 of the time about a loaded parent's OWN copy being absent
// and 0.2890 about the other parent's, measured on material where the answer is known from
// dissection. With one array loaded, roughly half of all losses were therefore answered by the
// direction that is wrong about seven times in ten. With both, each array is asked only the
// question it can answer.
//
// chr5 and chr13 are the same event with the parents swapped, which is what makes this a direction
// test rather than a detection test: a rule that inverted the two would name both wrong while
// every other assertion in this file still passed.
{
  const one = await run()
  const two = await run(undefined, true)

  // DETECTION MUST NOT MOVE. Loading a second array changes who is named, never what was found.
  const found = (r: Awaited<ReturnType<typeof run>>) =>
    (r.chroms ?? []).filter((c) => c.aneuploidy).map((c) => `${c.chrom}:${c.aneuploidy}`).sort()
  assert.deepEqual(found(two), found(one),
    'the second parental array is evidence about origin, not about copy number, and must not add '
    + `or remove an event. One parent: ${found(one).join(' ')}. Two: ${found(two).join(' ')}.`)

  const rowFor = (r: Awaited<ReturnType<typeof run>>, chrom: string) =>
    (r.oneParent ?? []).find((o) => o.where === `chr${chrom}`)

  // THE DIRECTIONS, both ways round.
  const patGone = rowFor(two, ONE_COPY)
  assert.ok(patGone, `chr${ONE_COPY} must carry a Mendelian row with both parents loaded`)
  assert.equal(patGone!.twoParents, true, 'and it must be marked as a two-parent row')
  assert.equal(patGone!.parent, 'paternal',
    `chr${ONE_COPY} was built with the FATHER's copy removed: every call is an allele he cannot `
    + `have supplied. Got ${patGone!.parent}, verdict ${patGone!.verdict}. If this reads maternal `
    + 'the rule is inverted and every event on every real run names the wrong parent.')

  const matGone = rowFor(two, MAT_LOST)
  assert.ok(matGone, `chr${MAT_LOST} must carry a Mendelian row`)
  assert.equal(matGone!.parent, 'maternal',
    `chr${MAT_LOST} was built with the MOTHER's copy removed. Got ${matGone!.parent}, verdict `
    + `${matGone!.verdict}.`)

  // AND THE SINGLE-PARENT RUN IS THE THING THIS REPLACES. With only the father loaded, chr13 is
  // the weak direction: his copy is present, so the channel has to infer the mother's absence from
  // heterozygosity that is not there, which is exactly what dropout produces. The row must not
  // silently carry a parent field it did not earn.
  const matGoneSolo = rowFor(one, MAT_LOST)
  assert.ok(matGoneSolo, 'the single-parent run must still produce the row')
  assert.notEqual(matGoneSolo!.twoParents, true,
    'a single-parent row must not claim to be a two-parent one')
  assert.equal(matGoneSolo!.parent, undefined,
    'and it must not carry an explicit parent, since with one array the parent is derived from '
    + 'the verdict and the role of the array that was loaded')

  // A GAIN'S ROW IS TRUE AND ITS ORIGIN IS STILL WITHHELD, and the two are not in tension.
  // chr11 is an isodisomic trisomy: three copies from the father, none from the mother. The
  // Mendelian channel is right that the maternal copy is absent. But with three copies present and
  // none of them hers, the EXTRA one is HIS, so passing that answer through as the origin of a
  // copy-gain would print the wrong parent at the confidence of a right one.
  const triRow = rowFor(two, ISO_TRISOMY)
  assert.equal(triRow?.parent, 'maternal',
    'the row states which copy is ABSENT, and on an isodisomic trisomy that is the maternal one')
  assert.equal(mendelParent(triRow!, 'paternal', true), null,
    'but the ORIGIN OF A GAIN is not the parent whose copy is absent: it is the other one, and '
    + 'this tool does not measure it. Naming the absent parent here is the sign inversion that '
    + 'has produced every backwards call in this codebase.')
  assert.equal(mendelParent(triRow!, 'paternal', false), 'maternal',
    'and the same row read as a LOSS does name that parent, so the gate is the class and not the '
    + 'row')

  // NEITHER PARENT IS NAMED WHERE BOTH COPIES ARE THERE. Every remaining autosome in this fixture
  // is intact, so anything named on one is a false call.
  const built = new Set([`chr${ONE_COPY}`, `chr${MAT_LOST}`, `chr${ISO_TRISOMY}`])
  const intact = (two.oneParent ?? []).filter((o) => !built.has(o.where) && !o.where.includes('Mb'))
  const falseCalls = intact.filter((o) => o.parent)
  assert.deepEqual(falseCalls.map((o) => `${o.where}=${o.parent}`), [],
    'a parent named on a chromosome carrying both copies is a false call')

  // AND THE ROW SAYS WHY, in words an operator can act on.
  assert.ok(/0\.8539/.test(patGone!.why),
    'the row has to carry the accuracy of the direction it used, or a reader cannot weigh it')
  assert.ok(!patGone!.why.includes('\u2014'), 'em dash')
}

console.log('scoreSample.check.ts: a rejected array reports nothing, a one-copy chromosome is not '
  + 'an isodisomy, a real isodisomy still is, a segment leaves with an origin row, and with both '
  + 'parents loaded each loss is named from the array that can see it')
