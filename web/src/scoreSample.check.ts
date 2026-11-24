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
import { relatednessAssessable, REL } from './relatedness.ts'
import { headerMap, parseRow, accumulate, accumulateBaf, emptyBafSums, finishProfile } from './ingest.ts'
import { emptyParent, collectParentRow, finishParent, emptyCollected, collectRow, scoreSample } from './scoreSample.ts'
import { mendelParent, runAlerts } from './defects.ts'
import { buildReportPdf } from './syngamyPdf.ts'
import { readFileSync } from 'node:fs'
import { PANEL_HET_FLOOR } from './parentage.ts'

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
 *  median is a whole-chromosome event and is not scanned for a segment inside it. The two
 *  homozygous chromosomes carry enough markers to expect LOH_MIN_EXPECTED_HET heterozygotes at
 *  this sample's own heterozygosity, which a real array's chromosome clears many times over. */
const PER_CHROM = (c: string) => (c === '9' ? 9_000 : c === '7' || c === '11' ? 4_000 : 2_400)
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
/** Biparental chromosomes carrying an extra copy: the paternal share moves to 0.64 or 0.36, which
 *  is two thirds shrunk the way whole-genome amplification shrinks it, and the intensity rises. */
const PAT_GAIN = '3'
const MAT_GAIN = '4'
/** The same shift with the intensity BELOW the cell's level, which an extra copy cannot do. */
const IMBALANCE = '8'
/** The same shift with the intensity only just above the cell's level: no extra copy is read. */
const WEAK = '10'

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
      // A GAIN IS A SHIFT OF THE ALLELE RATIO AT THE MARKERS WHERE THE PARENTS ARE OPPOSITE
      // HOMOZYGOTES, which is the only place a cell with both parents shows whose copy is extra.
      let baf = sg === 1 ? 0.5 + (r() - 0.5) * 0.08 : sg === 0 ? 0.02 + r() * 0.02 : 0.96 + r() * 0.02
      const opposite = pg !== 1 && mg !== 1 && pg !== mg
      if (opposite && (chrom === PAT_GAIN || chrom === MAT_GAIN || chrom === IMBALANCE || chrom === WEAK)) {
        // The paternal allele is A where the father is AA, so a paternal share of 0.64 is a B
        // allele frequency of 0.36 there and 0.64 where he is BB.
        const paternalShare = chrom === MAT_GAIN ? 0.36 : 0.64
        baf = (pg === 0 ? 1 - paternalShare : paternalShare) + (r() - 0.5) * 0.12
        sg = 1
      }
      if (chrom === PAT_GAIN || chrom === MAT_GAIN) lrr += 0.3
      else if (chrom === IMBALANCE) lrr -= 0.5
      else if (chrom === WEAK) lrr += 0.05
      sample.push([id, chrom, pos, lrr.toFixed(4), baf,
        '2.0', called ? String(sg) : '-1', '1'].join('\t'))
    }
  }
  // CHROMOSOME Y, so the parental-slot check has something to read. Without it `sexCall` returns
  // null and refuses to judge, which is the right behaviour on a panel with no Y coverage and is
  // exactly why the rest of this fixture is unaffected by that check. The father calls his Y at
  // full intensity; the mother has the same probes, all no-call, at the intensity of an absent
  // chromosome. That is what the two really look like.
  const TAB = String.fromCharCode(9)
  for (let i = 0; i < 400; i += 1) {
    const id = `Y_${i}`
    const pos = String(1000 + i * 10_000)
    parent.push([id, 'Y', pos, '0.0200', '0.9800', '1.0', 'AA', '1'].join(TAB))
    mother.push([id, 'Y', pos, '-3.5000', '', '0.0', 'NC', '1'].join(TAB))
    // The sample is a female embryo here: it took no Y. Only the PARENT slots are under test.
    sample.push([id, 'Y', pos, '-3.5000', '', '0.0', 'NC', '1'].join(TAB))
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

function indexOf(rows: ReturnType<typeof rowsOf>, inferred?: Parameters<typeof finishParent>[2]) {
  const pacc = emptyParent()
  const byChrom = new Map(); const bafSums = emptyBafSums(); let first = ''
  for (const row of rows) {
    if (!first) first = row.probesetId
    accumulate(row, byChrom); accumulateBaf(row, bafSums); collectParentRow(row, pacc)
  }
  return finishParent(pacc, finishProfile('p', byChrom, bafSums, first).build.build, inferred)
}

/** Any array in the parental slot, under any declared role. The slot check is about this pair. */
async function runWith(
  patRows: ReturnType<typeof rowsOf>,
  matRows: ReturnType<typeof rowsOf> | null,
  soloRole: 'paternal' | 'maternal',
  declaredStage?: 'failed',
) {
  const pat = indexOf(patRows)
  const mat = matRows ? indexOf(matRows) : null
  const acc = emptyCollected(pat, mat)
  const byChrom = new Map(); const bafSums = emptyBafSums(); let first = ''
  for (const row of sampleRows) {
    if (!first) first = row.probesetId
    accumulate(row, byChrom); accumulateBaf(row, bafSums); collectRow(row, pat, mat, acc)
  }
  const profile = finishProfile('sample', byChrom, bafSums, first)
  return scoreSample({
    acc, profile, pat, mat, soloRole, sibs: [], sampleName: 'sample',
    log: () => {}, declaredStage,
  })
}

/** The same run on another set of sample rows, both parents loaded or not, with a panel or not. */
async function runRows(
  rows: ReturnType<typeof rowsOf>, withMother: boolean,
  panel?: ReadonlyMap<string, 'AA' | 'BB'>, inferred?: Parameters<typeof finishParent>[2],
) {
  const pat = indexOf(parentRows, inferred)
  const mat = withMother ? indexOf(motherRows) : null
  const acc = emptyCollected(pat, mat, panel)
  const byChrom = new Map(); const bafSums = emptyBafSums(); let first = ''
  for (const row of rows) {
    if (!first) first = row.probesetId
    accumulate(row, byChrom); accumulateBaf(row, bafSums); collectRow(row, pat, mat, acc)
  }
  const profile = finishProfile('sample', byChrom, bafSums, first)
  return scoreSample({
    acc, profile, pat, mat, soloRole: 'paternal', sibs: [], sampleName: 'sample', log: () => {},
  })
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
  // WITH BOTH PARENTS LOADED, which is where the allele readings run at all. They are computed
  // before the stage verdict is acted on, so each has to be cleared with the rest.
  const both = await run('failed', true)
  assert.deepEqual((both.chroms ?? []).filter((c) => c.mixed || c.alleleGain || c.lostParent).map((c) => c.chrom), [],
    'a refused array carries no allele reading of any kind: absent parent, mixed signal or extra copy')
  assert.deepEqual(r.losses ?? [], [], 'and no loss rows')

  // A REFUSAL FOR DROP-IN IS NOT A REFUSAL OF THE ALLELES. The same genome with 30% of its
  // homozygous calls turned heterozygous: drop-in on that scale is what the cells the heterozygosity
  // ceiling refuses in GSE148488 carry, and it pushes this sample over the ceiling. It adds
  // heterozygous calls; it does not add an absent parent's allele as a homozygote, so the parent
  // whose copy is gone is still exactly as readable, and must still be named, correctly, on every
  // chromosome the clean run names it on. The intensity claims stay withheld: that is a separate
  // question this does not answer.
  {
    const dr = rng(20261110)
    const droppedIn = rowsOf(sample).map((row) => (
      /^\d+$/.test(row.chrom) && (row.genotype === 'AA' || row.genotype === 'BB') && dr() < 0.30
        ? { ...row, genotype: 'AB' as const } : row))
    const refused = await runRows(droppedIn as never, true)
    assert.equal(refused.stage?.stage, 'failed',
      `30% drop-in must take this sample over the heterozygosity ceiling, or this tests nothing. Got `
      + `${refused.stage?.stage}: ${refused.stage?.why}`)
    assert.equal(refused.stage?.refusedBy, 'heterozygosity-ceiling')
    const named = Object.fromEntries((refused.chroms ?? []).filter((c) => c.lostParent)
      .map((c) => [c.chrom, c.lostParent]))
    assert.deepEqual(named,
      { [ONE_COPY]: 'paternal', [ISODISOMY]: 'maternal', [ISO_TRISOMY]: 'maternal', [MAT_LOST]: 'maternal' },
      `the alleles name the same absent parents through the drop-in. Got ${JSON.stringify(named)}`)
    assert.deepEqual((refused.chroms ?? []).filter((c) => c.aneuploidy).map((c) => c.chrom), [],
      'while the intensity claims stay withheld on the refused array')
    assert.ok(!(refused.refusals ?? []).some((x) => x.claim === 'lostParent'),
      'and a claim that stands is not recorded as refused')
    // AND THE LABEL AGREES WITH THE RESULT. A banner telling the reader not to read parental calls,
    // above parental calls the result is showing, is the contradiction that sends a reader back to ask
    // which of the two to believe.
    assert.equal(refused.integrity?.level, 'disturbed')
    assert.ok(/which parent is read from the alleles/.test(refused.integrity?.headline ?? ''),
      `the label says what stands. Got: ${refused.integrity?.headline}`)
    assert.ok(!/Nothing on this sample is reportable/.test(refused.integrity?.limits ?? ''),
      'and does not tell the reader to ignore the calls it is showing')

    // WITH A FOUNDER PANEL THE SAME ARRAY IS NOT REFUSED AT ALL. The panel stands in for parents that
    // always agree; here it is the markers where this family's two parents share a homozygote, which
    // is what an unrelated panel approximates. The stage removes the drop-in and finds one genome.
    const panel = new Map<string, 'AA' | 'BB'>()
    {
      const mg = new Map(motherRows.map((r) => [r.probesetId, r.genotype]))
      for (const r of parentRows) {
        if ((r.genotype === 'AA' || r.genotype === 'BB') && mg.get(r.probesetId) === r.genotype
          && /^\d+$/.test(r.chrom)) panel.set(r.probesetId, r.genotype)
      }
    }
    const staged = await runRows(droppedIn as never, true, panel)
    assert.notEqual(staged.stage?.stage, 'failed',
      `with the panel, drop-in is removed and one genome is found. Got ${staged.stage?.why}`)
    assert.ok(Number.isFinite(staged.stage?.dropInCorrected), 'and the drop-in removed is stated')
    assert.deepEqual(Object.fromEntries((staged.chroms ?? []).filter((c) => c.lostParent)
      .map((c) => [c.chrom, c.lostParent])), named, 'the same absent parents are named')
    assert.deepEqual((staged.chroms ?? []).filter((c) => c.aneuploidy).map((c) => c.chrom), [],
      'copy number is still not read: no normal majority has been shown')
    // The claims that read a genome's copy number and need no parent. The extra-copy reading also
    // passes an intensity gate, but it names a parent and exists only with two loaded, like lostParent,
    // so it is a parent-of-origin claim and is not part of this comparison.
    const intensityRefusals = (r: Awaited<ReturnType<typeof runRows>>) => (r.refusals ?? [])
      .filter((x) => x.cause === 'no-intensity-reference' && (x.claim === 'aneuploidy' || x.claim === 'segment'))
      .map((x) => `${x.claim}:${x.chrom}`).sort()
    assert.ok(intensityRefusals(staged).includes(`aneuploidy:${ONE_COPY}`),
      `and the loss it would have reported is recorded under its physical cause. Got `
      + `${intensityRefusals(staged).join(' ')}`)

    // 1:1. NOTHING ABOUT THE STAGE OR THE INTENSITY MAY DEPEND ON HOW MANY PARENTS ARE LOADED. The
    // laboratory's material has one parent, so this is the property that makes the result usable there.
    const oneParent = await runRows(droppedIn as never, false, panel)
    assert.deepEqual(
      { stage: oneParent.stage?.stage, ploidy: oneParent.stage?.ploidy, d: oneParent.stage?.dropInCorrected },
      { stage: staged.stage?.stage, ploidy: staged.stage?.ploidy, d: staged.stage?.dropInCorrected },
      'the stage is identical with one parent and with two')
    // And with that parent reconstructed, whose banner states a heterozygosity no measured array
    // gave: a rebuild's estimate moves with its products and must not move the ceiling.
    const rebuilt = await runRows(droppedIn as never, false, panel, {
      hRetained: 0.30, spuriousAbsence: 0.005, contamination: 0.01, products: 8, side: 'paternal' })
    assert.deepEqual(
      { stage: rebuilt.stage?.stage, ploidy: rebuilt.stage?.ploidy, d: rebuilt.stage?.dropInCorrected },
      { stage: staged.stage?.stage, ploidy: staged.stage?.ploidy, d: staged.stage?.dropInCorrected },
      'a reconstructed parent\'s heterozygosity scaled the stage ceiling')
    const wholeRefusals = (r: Awaited<ReturnType<typeof runRows>>) =>
      intensityRefusals(r).filter((x) => x.startsWith('aneuploidy:'))
    assert.deepEqual(wholeRefusals(oneParent), wholeRefusals(staged),
      'and so is every whole-chromosome copy-number claim withheld')
    // And every segment refusal. The copy-number scan once skipped a chromosome the alleles named as
    // missing a parent, a name that exists only with two parents loaded, so segments differed between
    // the runs. It scans every chromosome now and reports a restated whole loss once.
    const segRefusals = (r: Awaited<ReturnType<typeof runRows>>) =>
      intensityRefusals(r).filter((x) => x.startsWith('segment:'))
    assert.deepEqual(segRefusals(oneParent), segRefusals(staged),
      'a segment refusal differs between one parent and two: the copy-number scan depends on the parents')
  }

  // AND EVERY ONE OF THOSE IS RECORDED, with what it would have said. The census scores refusals
  // against published truth from these records, so a cleared claim with no record is invisible.
  const refs = r.refusals ?? []
  assert.ok(refs.length > 0, 'a refused array records what it refused')
  for (const x of refs) {
    assert.ok(Object.keys(x.measured).length > 0, `a refusal of ${x.claim} carries no measurement`)
  }
  const lossRef = refs.find((x) => x.claim === 'aneuploidy' && x.chrom === ONE_COPY)
  assert.equal(lossRef?.withheld, 'loss',
    `chr${ONE_COPY} is reported as a loss when the array is read, so refusing the array must record `
    + `that loss as withheld. Got ${JSON.stringify(lossRef)}.`)
  assert.equal(lossRef?.gate, 'stage.declared', 'and name the gate that withheld it')
  assert.ok(refs.some((x) => x.claim === 'sample' && x.soft),
    'the uninterpretable label is recorded as a soft refusal of the whole sample')
  const bothRefs = both.refusals ?? []
  assert.ok(bothRefs.some((x) => x.claim === 'lostParent' || x.claim === 'alleleGain' || x.claim === 'mixed'),
    'with both parents loaded the allele readings are computed, cleared, and recorded as withheld')
  // AND A SAMPLE THAT IS READ CARRIES NO STAGE REFUSAL.
  const read = await run()
  assert.equal((read.refusals ?? []).filter((x) => x.gate.startsWith('stage.')).length, 0,
    'an array that resolves to a stage must carry no stage refusal')
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

  // A LOSS NEVER CARRIES THE ONE-PARENT GAIN SENTENCE, which says whose the extra copy is.
  assert.ok((one.losses ?? []).length > 0, 'the one-parent run has loss rows to read')
  for (const l of one.losses ?? []) {
    assert.ok(!/the extra copy is that parent/.test(l.why ?? ''), `a loss row names whose extra copy it is: ${l.where}: ${l.why}`)
  }

  // DETECTION MUST NOT MOVE. Loading a second array changes who is named, never what was found.
  const found = (r: Awaited<ReturnType<typeof run>>) =>
    (r.chroms ?? []).filter((c) => c.aneuploidy).map((c) => `${c.chrom}:${c.aneuploidy}`).sort()
  assert.deepEqual(found(two), found(one),
    'the second parental array is evidence about origin, not about copy number, and must not add '
    + `or remove an event. One parent: ${found(one).join(' ')}. Two: ${found(two).join(' ')}.`)

  // WHAT THE SECOND ARRAY DOES ADD IS THE ABSENT PARENT, FROM THE ALLELES ALONE. At markers where
  // the parents are opposite homozygotes, a chromosome reading one parent's allele and nothing
  // else is missing the other parent's copy. chr5 lost the loaded parent's copy, chr13 the
  // mother's, chr7 carries two copies from the father and none from the mother, and chr11 three
  // from the father and none from the mother: all four are missing a parent, and the alleles say
  // which, without saying how many copies remain. That is why chr7 keeps its isodisomy finding and
  // gains no loss: the copy number is the intensity's to decide, and on chr7 the intensity reads
  // two. On chr11 it reads three, and the alleles then say whose all three are.
  const lostParentOf = (r: Awaited<ReturnType<typeof run>>) =>
    Object.fromEntries((r.chroms ?? []).filter((c) => c.lostParent).map((c) => [c.chrom, c.lostParent]))
  assert.deepEqual(lostParentOf(one), {},
    'with one parent loaded the other allele has no owner, so no parent is named from the alleles')
  assert.deepEqual(lostParentOf(two),
    { [ONE_COPY]: 'paternal', [ISODISOMY]: 'maternal', [ISO_TRISOMY]: 'maternal', [MAT_LOST]: 'maternal' },
    `the alleles name the absent parent on every chromosome that is missing one. Got ${JSON.stringify(lostParentOf(two))}`)
  assert.equal((two.chroms ?? []).find((c) => c.chrom === ISODISOMY)?.aneuploidy, undefined,
    'and an absent parent at normal intensity is not a loss: it is the isodisomy the findings carry')
  assert.ok((two.oneParent ?? []).some((o) => o.where === `chr${ISODISOMY}`),
    'a chromosome the alleles flagged reaches the Mendelian channel like any detected event')

  // AN EXTRA COPY, AND WHOSE, FROM THE ALLELES. chr3 carries an extra paternal copy and chr4 an
  // extra maternal one, both with the intensity up, so the direction test is the same as for the
  // losses above: a rule that inverted the two would name both wrong while every other assertion
  // here still passed. chr8 shifts the same way with the intensity DOWN, which an extra copy
  // cannot do, and must name nobody.
  const ag = (r: Awaited<ReturnType<typeof run>>, c: string) => (r.chroms ?? []).find((x) => x.chrom === c)?.alleleGain
  assert.equal(ag(one, PAT_GAIN), undefined, 'with one parent loaded nothing is read from the shares')
  assert.equal(ag(two, PAT_GAIN)?.verdict, 'gain', `chr${PAT_GAIN} carries an extra paternal copy. Got ${JSON.stringify(ag(two, PAT_GAIN))}`)
  assert.equal(ag(two, PAT_GAIN)?.parent, 'paternal', 'and the extra copy is the father\'s')
  assert.equal(ag(two, MAT_GAIN)?.verdict, 'gain', `chr${MAT_GAIN} carries an extra maternal copy. Got ${JSON.stringify(ag(two, MAT_GAIN))}`)
  assert.equal(ag(two, MAT_GAIN)?.parent, 'maternal', 'and the extra copy is the mother\'s')
  assert.equal(ag(two, IMBALANCE)?.verdict, 'imbalance', 'a share shift with the intensity below the cell is an imbalance')
  assert.equal(ag(two, IMBALANCE)?.parent, undefined, 'and it names no parent')
  assert.equal(ag(two, WEAK)?.verdict, 'imbalance', 'a shift with the intensity under the floor is not an extra copy')
  assert.equal(ag(two, WEAK)?.parent, undefined, 'and names no parent either')
  const called = (two.chroms ?? []).filter((c) => c.alleleGain && ['gain', 'imbalance'].includes(c.alleleGain.verdict))
    .map((c) => c.chrom).sort()
  assert.deepEqual(called, [PAT_GAIN, MAT_GAIN, IMBALANCE, WEAK].sort(),
    `no other chromosome of this fixture is an allele gain or imbalance. Got chr${called.join(', chr')}`)

  // AND NOT ON A GENOME THAT CARRIES ONE PARENT THROUGHOUT. The father's own array as the sample,
  // both parents loaded: every autosome reads his allele alone, and no chromosome of it is missing
  // a parent in any sense the row could carry. Found on the bundled example, where this put eight
  // rows on the sperm donor's replicate array.
  const himself = await runRows(parentRows, true)
  assert.deepEqual(lostParentOf(himself), {},
    `the father's own genome names no absent parent chromosome by chromosome, got ${JSON.stringify(lostParentOf(himself))}`)

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
  // AND IT NOW NAMES HER, from two statistics that need no second array: whether this region
  // carries one contribution or two, and whether the copy that remains is the loaded father's.
  // See nameAbsentParent. The direction has to be right, not merely present: if this reads
  // paternal then every one-array run names the wrong parent on the commoner of the two cases.
  assert.equal(matGoneSolo!.parent, 'maternal',
    `chr${MAT_LOST} was built with the MOTHER's copy removed and only the father is loaded. Got `
    + `${matGoneSolo!.parent}. Measured on 110 hand-called losses of the public corpus this names `
    + '94 correctly and 1 wrongly; a flip here is that measurement inverted.')
  assert.equal(matGoneSolo!.oneArrayNamed, true,
    'and the row must say the name came from one array, since a reader cannot otherwise tell it '
    + 'from a Mendelian call made against both parents')
  assert.ok(/opposite homozygote/.test(matGoneSolo!.why ?? ''),
    'the reason must state the evidence, so the call can be checked by hand')

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

  // NEITHER PARENT IS NAMED WHERE BOTH PARENTS' COPIES ARE THERE. Every remaining autosome in this
  // fixture is intact, so anything named on one is a false call. chr7 is built: two copies, both
  // the father's, so its row rightly says the maternal copy is absent.
  const built = new Set([`chr${ONE_COPY}`, `chr${MAT_LOST}`, `chr${ISO_TRISOMY}`, `chr${ISODISOMY}`])
  const intact = (two.oneParent ?? []).filter((o) => !built.has(o.where) && !o.where.includes('Mb'))
  const falseCalls = intact.filter((o) => o.parent)
  assert.deepEqual(falseCalls.map((o) => `${o.where}=${o.parent}`), [],
    'a parent named on a chromosome carrying both copies is a false call')

  // AND THE ROW SAYS WHY, in words an operator can act on.
  assert.ok(/0\.8539/.test(patGone!.why),
    'the row has to carry the accuracy of the direction it used, or a reader cannot weigh it')
  assert.ok(!patGone!.why.includes('\u2014'), 'em dash')
}


// ------------------------------------------------------- THE WRONG ARRAY IN THE PARENTAL SLOT
//
// An operator opens a folder holding a sperm donor and four egg donors, picks the wrong one, and
// marks it paternal. Every genotype in that file is real and every calculation downstream is
// arithmetically correct. What is wrong is the LABEL, and the label is the answer: every event
// comes out named for the parent it did not come from.
//
// Measured on the 12 donor arrays of GSE148488, roles from the GEO titles rather than from this
// tool: 4 sperm donors read a Y call ratio of 1.029 to 1.039, and all 8 egg donors read exactly
// 0.000. Twelve of twelve, with nothing in between.
{
  // The MOTHER's array in the paternal slot. She is a real, complete, high-quality array.
  const swapped = await runWith(motherRows, null, 'paternal')
  assert.equal(swapped.parentSanity?.conflict, 'paternal-array-has-no-y',
    `an array with no Y in the paternal slot must be caught. Got `
    + `${JSON.stringify(swapped.parentSanity)}`)

  // AND NOTHING MAY BE NAMED. This is the whole point: the events survive, the names do not.
  const named = [
    ...(swapped.oneParent ?? []).filter((o) => o.parent
      || o.verdict === 'known-parent-lost' || o.verdict === 'other-parent-lost')
      .map((o) => `mendel ${o.where}`),
    ...(swapped.dosageCalls ?? []).filter((d) => (d as { parent?: string }).parent)
      .map((d) => `dosage ${d.where}`),
  ]
  assert.deepEqual(named, [],
    `no parent may be named when the parental slot is wrong. Got: ${named.join(', ')}`)

  // AND THE EVENTS SURVIVE. Refusing the whole run would throw away real findings: they were
  // found by channels that never read the parental array's sex.
  const events = (swapped.chroms ?? []).filter((c) => c.aneuploidy).length
  assert.ok(events > 0,
    'the copy-number events must still be reported: a wrong parental slot corrupts the naming, '
    + `not the detection. Found ${events}.`)

  // AND IT REACHES A SURFACE. Every one of these checks was computed and displayed nowhere.
  const alerts = runAlerts(swapped)
  assert.ok(alerts.length > 0, 'the conflict must appear in the shared alert list')
  assert.ok(alerts.some((a) => /CARRIES NO CHROMOSOME Y/.test(a.headline)),
    'the sex conflict must be one of them, wherever it sits in the order: '
    + alerts.map((a) => a.headline).join(' | '))

  // THE CONTROL. The father in the paternal slot is not a conflict, and still names parents.
  const fine = await runWith(parentRows, null, 'paternal')
  assert.equal(fine.parentSanity?.conflict, undefined,
    `the real father in the paternal slot must NOT trip the check. `
    + `Got ${JSON.stringify(fine.parentSanity)}`)
  const fineNamed = (fine.oneParent ?? []).filter((o) =>
    o.verdict === 'known-parent-lost' || o.verdict === 'other-parent-lost').length
  assert.ok(fineNamed > 0,
    'and the correct run must still name parents, or the guard has simply broken the tool')

  // THE MIRROR. The father declared MATERNAL is the same mistake the other way round.
  const mirrored = await runWith(parentRows, null, 'maternal')
  assert.equal(mirrored.parentSanity?.conflict, 'maternal-array-has-y')

  // AND IT REACHES THE REPORT, which is the artefact that leaves the lab. A warning that exists
  // only in a JSON field is not a warning; the page, the report and the command line each render
  // the same list, and this pins the one that gets printed and filed.
  const pdf = await buildReportPdf({
    files: [{
      name: 'sample', size: 1, sha256: 'x', role: 'sample', markers: 1000, result: swapped,
    }],
    donorHeterozygosity: 0.17, startedAt: null, generatedAt: '2026-01-01 00:00:00Z',
    tool: 'OriginMarker', reportId: 'CHECK', fromExamples: true,
  } as never)
  const txt = Buffer.from(await pdf.arrayBuffer()).toString('latin1')
  assert.ok(txt.startsWith('%PDF-') && txt.includes('%%EOF'),
    'the report must still build on a run whose parental slot is wrong, rather than throwing')
  assert.ok(/CHROMOSOME Y/.test(txt),
    'and the conflict has to be IN it: a report that prints the findings without the reason they '
    + 'carry no parent is the failure this whole check is about')
}

console.log('scoreSample.check.ts: a rejected array reports nothing, a one-copy chromosome is not '
  + 'an isodisomy, a real isodisomy still is, a segment leaves with an origin row, and with both '
  + 'parents loaded each loss is named from the array that can see it')

// ---------------------------------------------------------------- the verdict the material carries
//
// THE MEASUREMENT THAT FORCED THIS, over 106 children and 318 pairs in
// audit/relationship-separability.ts. On blastomeres the shipped verdict read `unrelated` for the
// child's own confirmed father 9 times out of 9, and `unrelated` for a stranger 9 times out of 9.
// The same word for opposite facts is worse than no word: it reads as a finding either way. Only
// bulk separates across both independent measurements of this, so only bulk carries a verdict.
{
  const r = await run()
  const stage = r.stage?.stage
  assert.ok(r.relationship, 'the relationship block must still be produced')
  assert.ok(Number.isFinite(r.relationship!.oppositeHomRate),
    'the RATE is a measurement and is reported whatever the material')
  if (relatednessAssessable(stage)) {
    assert.notEqual(r.relationship!.verdict, REL.notAssessable,
      `${stage} separates a parent from a stranger, so it must carry a verdict`)
  } else {
    assert.equal(r.relationship!.verdict, REL.notAssessable,
      `${stage} does not separate a true parent from a stranger, so a verdict of `
      + `"${r.relationship!.verdict}" states something the array cannot support`)
  }
  // And the split is on the material rather than on the number, so a stranger and a true parent on
  // the same material are treated the same way. Withholding only the awkward answers would be
  // worse than withholding none.
  assert.equal(relatednessAssessable('trophectoderm'), false,
    'trophectoderm overlapped in both measurements and must not carry a verdict')
  assert.equal(relatednessAssessable('blastomere'), false,
    'blastomere had a 0.0021 pooled margin, which is two samplings not overlapping yet')
  assert.equal(relatednessAssessable('bulk'), true,
    'bulk separated by 0.0389 and is where OPPOSITE_HOM_MAX was measured')
}

// ------------------------------------------------- no parent is named on a genome of unknown ploidy
//
// THE BUG THIS PINS. The Mendelian channel was gated with `!zygosity.startsWith('uniparental')`,
// which reads as "biparental" and is not: it is also true for `unknown`, which is the tool
// REFUSING to say how many parental contributions are present. So a genome whose ploidy could not
// be established still had parents named on it, by a channel whose own comment says it names the
// WRONG parent on a one-parent genome. Four rows on this fixture, at full confidence.
//
// THE FIXTURE IS THE REAL SITUATION, not a contrived one. Removing the B-allele frequencies is
// what a genotype-only export looks like, and it is what the converted second platform effectively
// is: the band cannot be read, the genotype fallback lands between its boundaries, and zygosity
// comes back `unknown` while the array still resolves to a stage and still carries its events.
{
  /** The same sample with no B-allele frequency column. Nothing else is touched. */
  const noBaf = (() => {
    const ls = sample.split('\n')
    const out = [ls[0]]
    for (let i = 1; i < ls.length; i += 1) {
      if (!ls[i]) continue
      const f = ls[i].split('\t')
      f[4] = ''
      out.push(f.join('\t'))
    }
    return out.join('\n')
  })()

  const pat = indexOf(parentRows)
  const acc = emptyCollected(pat, null)
  const byChrom = new Map(); const bafSums = emptyBafSums(); let first = ''
  for (const row of rowsOf(noBaf)) {
    if (!first) first = row.probesetId
    accumulate(row, byChrom); accumulateBaf(row, bafSums); collectRow(row, pat, null, acc)
  }
  const profile = finishProfile('nobaf', byChrom, bafSums, first)
  const out = await scoreSample({
    acc, profile, pat, mat: null, soloRole: 'paternal', sibs: [], sampleName: 'nobaf',
    log: () => {},
  })

  // The fixture must reproduce the situation, or the assertion below proves nothing.
  assert.equal(out.zygosity, 'unknown',
    `this fixture exists to exercise a REFUSED zygosity, got ${out.zygosity}`)
  assert.notEqual(out.stage?.stage, 'failed',
    'the array must still resolve to a stage, or the run is refused for an unrelated reason')
  assert.ok((out.segments?.length ?? 0) + (out.chroms ?? []).filter((c) => c.aneuploidy).length > 0,
    'the array must still carry events, or there is nothing for the channel to name a parent on')

  // THE REGRESSION. Detection is unaffected by ploidy and is kept; ATTRIBUTION is not.
  assert.equal((out.oneParent ?? []).length, 0,
    'no parent may be named on a genome whose ploidy the tool declined to call: '
    + `got ${(out.oneParent ?? []).length} rows at zygosity ${out.zygosity}`)

  // THE CONTROL. The same fixture WITH its B-allele frequencies reads diploid and still names
  // parents, so the guard above is withholding on ploidy rather than switching the channel off.
  const ok = await run()
  assert.equal(ok.zygosity, 'diploid',
    `the control must be a known two-parent genome, got ${ok.zygosity}`)
  assert.ok((ok.oneParent ?? []).length > 0,
    'a known biparental genome must still have its parents named, or this guard has silenced the '
    + 'channel outright rather than restricting it')
}


// --- A RECONSTRUCTED PARENT REACHES THE SCORER WITH ITS OWN PARAMETERS, OR NOT AT ALL -----------
//
// THE PATTERN THIS GUARDS, which is now three defects deep. In each case the CAPABILITY existed,
// a test covered the FUNCTION, and the CALLER skipped it, so every test passed while the shipped
// path did the wrong thing:
//
//   1. `om reconstruct --out` hand-rolled an array header instead of calling inferredArrayText,
//      writing a reconstructed person's genotype with no banner and no mark (cli/om.check.ts).
//   2. `scoreSample` called classify() without `spuriousAbsence`, so every sample scored against
//      a reconstruction was judged by the ceiling calibrated for a measured array, even though
//      parentage.check.ts section 16 covers that term and passes.
//   3. `finishParent` derived parental heterozygosity from a reconstruction's own genotypes,
//      which are homozygous at every marker, so it read 0.0% against a banner stating 17.332%.
//      PANEL_HET_FLOOR then withheld the second-parent axis for the whole run. The guard was
//      right; it was being fed a number that could not be anything but zero.
//
// So this asserts the SEAM rather than the three instances: anything that reduces a parental file
// to a ParentIndex must be able to accept that file's banner, and the figures a reconstruction
// states about itself must arrive at the two places that consume them. A fourth loader added
// later without the banner argument fails here rather than shipping silently.
{
  const PARENT_SOURCES = ['../../cli/om.ts', '../../audit/lab-zygote-check.ts', './Syngamy.tsx']
  for (const rel of PARENT_SOURCES) {
    const src = readFileSync(new URL(rel, import.meta.url), 'utf8')
    if (!src.includes('finishParent(')) continue
    assert.ok(/parseInferredBanner\(/.test(src),
      `${rel} builds a ParentIndex but never reads a reconstruction's banner, so a rebuilt `
      + 'reference loaded there reports 0% parental heterozygosity and silently loses the '
      + 'second-parent axis. Pass parseInferredBanner(head) as the third argument.')
  }

  // The scorer must pass BOTH figures on to the classifier. Named explicitly because they are
  // separate omissions with opposite effects: the floor makes absence over-called, the
  // heterozygosity makes everything keyed to it under-called.
  const scorer = readFileSync(new URL('./scoreSample.ts', import.meta.url), 'utf8')
  assert.ok(/spuriousAbsence: pat\.inferred\.spuriousAbsence/.test(scorer),
    "scoreSample must hand the reference's own absence floor to classify()")
  assert.ok(/inferredParent: true/.test(scorer),
    'scoreSample must tell classify() the parent is reconstructed, or the zero-heterozygosity '
    + 'limit misreports it as a monomorphic panel')

  // THE SUBSTITUTION ITSELF, end to end: a banner-bearing parent must arrive at classify() with
  // the banner's heterozygosity, and a measured one must be untouched.
  const banner = {
    hRetained: 0.17332, spuriousAbsence: 0.01374, contamination: 0.02749,
    products: 6, side: 'maternal',
  }
  const acc = emptyParent()
  // One homozygous row, which is all a reconstruction ever contains.
  collectParentRow({
    probesetId: 'AX-1', chrom: '1', pos: 1000, genotype: 'AA', baf: null, log2R: null,
  } as never, acc as never)
  const rebuilt = finishParent(acc as never, null, banner as never)
  assert.equal(rebuilt.measuredHeterozygosity, 0,
    'a reconstruction is homozygous by construction, so the file itself must read exactly zero')
  assert.equal(rebuilt.heterozygosity, 0.17332,
    'and the figure the reconstruction established must replace it')
  assert.ok(rebuilt.heterozygosity >= PANEL_HET_FLOOR,
    'which must clear PANEL_HET_FLOOR, or the second-parent axis stays withheld and the repair '
    + 'has achieved nothing. The guard is not relaxed; it is given the true figure')

  const measured = finishParent(acc as never, null)
  assert.equal(measured.heterozygosity, 0,
    'without a banner nothing is substituted: an unreadable reconstruction withholds rather than '
    + 'guessing, and a measured array is unaffected')
  assert.equal(measured.inferred, undefined, 'a measured array carries no banner')

  // An unreadable banner must NOT quietly become a measured array.
  const unreadable = finishParent(acc as never, null,
    { ...banner, hRetained: NaN } as never)
  assert.equal(unreadable.heterozygosity, 0,
    'a banner without a usable heterozygosity leaves the measured zero in place')
  assert.ok(unreadable.inferred, 'but the file is still known to be a reconstruction')
}

// WHICH OF THE TWO CAUSES, decided on the pipeline rather than on the helper. Chromosome 3 is
// handed back with almost every marker no-called while the parents' genotypes are left untouched,
// so the markers are on offer and only the sample failed to read them. The refusal has to say that,
// because the action it calls for is rerunning the array, and the other refusal's action is
// accepting that nothing here is readable.
{
  const everyNth = 50
  const blanked = rowsOf(sample).map((row, i) => (
    row.chrom === '3' && i % everyNth !== 0 ? { ...row, genotype: 'NC' as const } : row))
  const out = await runRows(blanked as never, true) as unknown as {
    chroms: { chrom: string; contribution?: unknown; contributionWithheld?: string }[]
    unmeasured: { chrom: string; why: string }[]
  }
  // It is not in `chroms`, and it must not be: 48 markers is not a measurement and its rate would
  // widen the null every other chromosome is read against. It has to be somewhere, though.
  assert.ok(!out.chroms.some((c) => c.chrom === '3'),
    'a chromosome at a fiftieth of its markers stays out of the measurement list')
  const three = out.unmeasured.find((u) => u.chrom === '3')
  assert.ok(three, 'but it is reported as unmeasured rather than vanishing from the result')
  const why = three!.why
  assert.ok(why, 'and it says why')
  assert.ok(/did not read them/.test(why!),
    'the parents offered these markers and the sample lost them, so the array is what is named. '
    + `Got: ${why}`)
  assert.ok(!/No array can do better/.test(why!),
    'and it must not tell the reader the parents are the limit when they are not')

  // Both counts appear, because the gap between them is the evidence for the claim.
  const nums = [...why!.matchAll(/([\d,]+)/g)].map((m) => Number(m[1].replace(/,/g, '')))
  assert.ok(nums.length >= 2, `both the offer and the return are stated. Got: ${why}`)
  assert.ok(Math.max(...nums) > 10 * Math.min(...nums.filter((n) => n > 0)),
    `the offer must be far above what came back, or this is not the array's failure. Got: ${why}`)

  // The untouched chromosomes are unaffected: this is a statement about chromosome 3 only.
  const one = out.chroms.find((c) => c.chrom === '1')
  assert.ok(one && !one.contributionWithheld,
    'a chromosome the sample read normally still carries its call')
  assert.equal(out.unmeasured.length, 1,
    'and no other chromosome is reported unmeasured')

  // A clean run reports nothing unmeasured, or the list would be noise on every array.
  const clean = await runRows(rowsOf(sample) as never, true) as unknown as {
    unmeasured: { chrom: string; why: string }[]
  }
  assert.equal(clean.unmeasured.length, 0,
    'an array that genotyped every chromosome reports none unmeasured')
}
console.log('scoreSample.check.ts: a refusal names the array when the array is what failed')

// ---------------------------------------------------------------- the report says what the card says
// The sex-chromosome line the page prints is in the report, and a sample not reportable prints no
// change in its bar and marks every row it lists as not reported.
{
  const { cardSays, sexLine } = await import('./pageSays.ts')
  const r = await run()
  r.chroms.push({ ...r.chroms[0], chrom: 'X' })
  r.sexChromosome = { ...(r.sexChromosome ?? {}), constitution: '47,XXY or a mix of cells', aneuploid: false,
    why: 'chrY is present; chrX reads two different X by its heterozygous calls' } as never
  const pdfOf = async () => Buffer.from(await (await buildReportPdf({
    files: [{ name: 'sample', size: 1, sha256: 'x', role: 'sample', markers: 1000, result: r }],
    donorHeterozygosity: 0.17, startedAt: null, generatedAt: '2026-01-01 00:00:00Z',
    tool: 'OriginMarker', reportId: 'CHECK', fromExamples: true,
  } as never)).arrayBuffer()).toString('latin1')
  const txt = await pdfOf()
  assert.ok(txt.includes(`(${sexLine(r.sexChromosome!)})`), 'the report prints the sex-chromosome line')
  assert.ok(txt.includes('47,XXY or a mix of cells: chrY is present'))

  // This fixture is one the run itself calls uninterpretable, so it is the not-reportable case.
  assert.equal(r.integrity?.level, 'uninterpretable')
  const says = cardSays(r, undefined, r.originClass, false)
  const nr = says.list.notReported!
  assert.ok(nr.rows.length > 0, 'the fixture must list rows for this to test anything')
  const out = await pdfOf()
  assert.equal(out.match(/\(not reported\)/g)?.length, nr.rows.length, 'every listed row marked not reported')
  assert.ok(out.includes('(not reportable)'), 'the summary column says not reportable')
  assert.ok(out.includes(`(${says.head.word})`), 'the bar is the headline word with no count after it')
  // With no parent loaded the sample line names no parent.
  const parentLine = async (x: typeof r) => Buffer.from(await (await buildReportPdf({
    files: [{ name: 'sample', size: 1, sha256: 'x', role: 'sample', markers: 1000, result: x }],
    donorHeterozygosity: 0.17, startedAt: null, generatedAt: '2026-01-01 00:00:00Z',
    tool: 'OriginMarker', reportId: 'CHECK', fromExamples: true,
  } as never)).arrayBuffer()).toString('latin1')
  const none = await parentLine({ ...r, noParent: true })
  assert.ok(none.includes('(sample no parent loaded)') && !none.includes('vs sperm donor'), 'no parent is named with none loaded')
  assert.ok((await parentLine(r)).includes('(sample vs sperm donor)'), 'a loaded father is still named')
  // With no parent at all nothing is inferred about the oocyte donor; with the father alone it is.
  assert.ok(!none.includes('inferred from paternal absence'), 'no inference with no parent loaded')
  const withFather = Buffer.from(await (await buildReportPdf({
    files: [{ name: 'sperm', size: 1, sha256: 'y'.repeat(64), role: 'donor', markers: 1000 },
      { name: 'sample', size: 1, sha256: 'x', role: 'sample', markers: 1000, result: r }],
    donorHeterozygosity: 0.17, startedAt: null, generatedAt: '2026-01-01 00:00:00Z',
    tool: 'OriginMarker', reportId: 'CHECK', fromExamples: true,
  } as never)).arrayBuffer()).toString('latin1')
  assert.ok(withFather.includes('inferred from paternal absence'), 'the father alone keeps the inference sentence')
  // The methods text too: with no parent it names no donor and infers nothing.
  const flatText = (t: string) => [...t.matchAll(/\(((?:\\.|[^\\)])*)\) Tj/g)].map((x) => x[1]).join(' ').replace(/\s+/g, ' ')
  const np = flatText(none)
  assert.ok(np.includes('No parental array was supplied.') && np.includes('No parental origin is named'), 'no-parent methods')
  assert.ok(!/inferred from|sperm donor's|No oocyte donor array|donor's obligate|paternal absence/.test(np),
    np.match(/.{60}(inferred from|sperm donor's|No oocyte donor array|donor's obligate|paternal absence).{40}/)?.[0])
  assert.ok(flatText(withFather).includes('No oocyte donor array was supplied. A maternal origin below is therefore inferred'),
    'the father alone keeps the methods sentence')
  console.log('  the report carries the sex-chromosome line, and a sample not reportable states no change')
}

// ---------------------------------------------------------------- the report's outputs: withheld, bars, methods
// The report lists every refusal and row not reported with its cause and their count, the bars of
// the registry that this setup applied and no other, and the methods paragraph with its citations.
{
  const { cardSays } = await import('./pageSays.ts')
  const { withheldRows, withheldHeading } = await import('./withheld.ts')
  const { barsApplied, barValue } = await import('./bars.ts')
  const { CHECKED } = await import('./methods.ts')
  const { reportMethods } = await import('./syngamyPdf.ts')
  const strings = (txt: string) => [...txt.matchAll(/\(((?:\\.|[^\\)])*)\) Tj/g)].map((m) => m[1]
    .replace(/\\([0-7]{3})|\\(.)/g, (_, o, c) => (o ? String.fromCharCode(parseInt(o, 8)) : c)))
  const report = async (r: Awaited<ReturnType<typeof run>>, father: boolean, tool = 'OriginMarker 9.9.9 (Check)') => {
    const files: never[] = [{ name: 'sample', size: 1, sha256: 'x', role: 'sample', markers: 1000, result: r } as never]
    if (father) files.unshift({ name: 'donor', size: 1, sha256: '0'.repeat(64), role: 'donor', markers: 1 } as never)
    const input = { files, donorHeterozygosity: 0.17, startedAt: null, generatedAt: '2026-01-01 00:00:00Z',
      tool, reportId: 'CHECK', fromExamples: false }
    const txt = Buffer.from(await (await buildReportPdf(input as never)).arrayBuffer()).toString('latin1')
    return { txt, s: strings(txt), all: strings(txt).join(' ').replace(/\s+/g, ' '), methods: reportMethods(input as never) }
  }
  const withF = await run()
  const none = { ...withF, noParent: true } as typeof withF
  for (const [r, father] of [[withF, true], [none, false]] as const) {
    const { txt, s, all, methods } = await report(r, father)
    const held = withheldRows(r, cardSays(r, undefined, r.originClass, false).list)
    assert.equal(held.length, (r.refusals ?? []).length
      + (cardSays(r, undefined, r.originClass, false).list.notReported?.rows.length ?? 0),
    'one withheld row per refusal and per row not reported')
    assert.ok(s.includes(withheldHeading(held.length)), 'the report carries the withheld count')
    for (const w of held) assert.ok(all.includes(w.what.split(' ').slice(0, 3).join(' ')), `withheld row printed: ${w.what}`)
    const applied = barsApplied({ parents: father ? 1 : 0, run: false, comparison: false })
    const listed = s.filter((x) => applied.some((b) => b.name === x))
    assert.deepEqual([...new Set(listed)].sort(), [...new Set(applied.map((b) => b.name))].sort(), 'the report lists exactly the bars applied')
    for (const b of applied) assert.ok(s.includes(barValue(b)), `${b.name} prints its value`)
    if (!father) {
      assert.ok(!s.includes('ABSENCE_MARGIN') && !s.includes('SURVIVOR_MIN'), 'no parental bar with no parent')
      assert.ok(methods.text.includes('No parental array was loaded.') && !/father|mother|parent's/.test(methods.text),
        'no parental sentence with no parent')
    } else assert.ok(methods.text.includes("The father's array was loaded"))
    assert.ok(all.includes(methods.text.split('. ')[1]), 'the methods paragraph is printed')
    assert.ok(methods.cites.every((c) => Object.values(CHECKED).includes(c)), 'citations only from the checked list')
    assert.ok(!methods.cites.includes(CHECKED.verdyck2022), 'the two-parent citation only with both parents')
    assert.ok(!/units of one embryo/.test(methods.text), 'no origin sentence for one cell')
    assert.ok(!/\/Image|\/XObject|\sBI\s/.test(txt), 'the report holds no image')
    assert.ok(!s.includes('Constants and where they come from'))
    const pages = (txt.match(/\/Type \/Page /g) ?? []).length
    assert.equal(s.filter((x) => x.includes('OriginMarker 9.9.9 (Check)') && x.includes(' | page ')).length, pages,
      'every page carries the version')
  }
  const bare = await report(withF, true, 'OriginMarker')
  assert.ok(bare.s.some((x) => x.includes('OriginMarker, version not reported | report')), 'no version given is said')
  console.log('  the report lists what was withheld, the bars applied and the methods paragraph')
}
