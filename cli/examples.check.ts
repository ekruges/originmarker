// The bundled examples still print what their own lines say. Run: node cli/examples.check.ts
//
// This guards an ANSWER rather than a build. A new person's first contact with the tool is an
// Examples button, and each card shows its example's line beside what the run printed, so a line
// that has stopped being true is the first thing they read wrong.
//
// It runs each set through cli/exampleRun.ts, which assembles the page's run from the page's own
// modules and words the card with pageSays, and pins every card. No fixtures: the files a user
// downloads are the files under test.
import assert from 'node:assert/strict'
import { existsSync, readFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { EXAMPLE_FILES, EXAMPLE_MARKERS, EXAMPLE_SETS } from '../web/src/examples.ts'
import { runExamples, summary } from './exampleRun.ts'

const DIR = new URL('../web/public/examples/', import.meta.url).pathname

// --- 1. every listed file is there, at the advertised density -----------------------------------
for (const e of EXAMPLE_FILES) {
  const path = `${DIR}${e.file}`
  assert.ok(existsSync(path), `${e.file} is listed in examples.ts but not present in public/examples`)
  assert.ok(e.file.startsWith(`${e.gsm}_`), `${e.file}: the name must start with its accession`)
  const rows = gunzipSync(readFileSync(path)).toString('utf8').split('\n').filter(Boolean).length - 1
  assert.equal(rows, EXAMPLE_MARKERS, `${e.gsm}: ${rows} markers, examples.ts advertises ${EXAMPLE_MARKERS}`)
}
const rank = { donor: 0, oocyte: 1, sample: 2 }
for (const set of EXAMPLE_SETS) {
  const order = set.examples.map((e) => rank[e.role])
  assert.deepEqual(order, [...order].sort(), `${set.id}: parents must load before samples`)
}

/**
 * Per sample, what its card must print: the headline word, the parent checks, the sex chromosomes,
 * and every row of the change list with the word on it. The quiet rows are pinned too, because an
 * unconfirmed row a reader sees is part of what the card says.
 */
const PIN: Record<string, Record<string, {
  head: string; father?: string; mother?: string; sex: string
  counted: string[]; once: string[]; seen: string[]
}>> = {
  two: {
    'GSM4472410_embryoD_TE.half.csv': {
      head: 'BOTH PARENTS', father: 'carries', mother: 'carries', sex: 'XX',
      counted: [], once: [], seen: [],
    },
    'GSM4472412_embryoE_TE.half.csv': {
      head: 'BOTH PARENTS', father: 'carries', mother: 'carries', sex: 'XY',
      counted: ['chr6 64.6-170.9Mb loss paternal | found two ways'],
      once: ['chr6 66.4-170.9Mb loss maternal | one reading'], seen: [],
    },
    'GSM4552364_embryoA_biopsy1.half.csv': {
      head: 'BOTH PARENTS', father: 'carries', mother: 'carries', sex: 'XX',
      counted: [
        'chr2 whole loss paternal | consistent',
        'chr5 149.4-180.7Mb loss maternal | found two ways',
        'chr6 112.0-170.7Mb loss maternal | found two ways',
        'chr7 whole loss maternal | consistent',
        'chr9 whole loss maternal | consistent',
        'chr11 whole loss paternal | consistent',
        'chr15 whole loss paternal | consistent',
        'chr16 47.8-87.7Mb loss paternal | found two ways',
        'chr17 whole loss maternal | consistent',
        'chr19 whole loss paternal | consistent',
      ],
      once: ['chr4 184.6-190.8Mb loss paternal | one reading'],
      seen: ['chr16 whole loss unclear | also seen'],
    },
  },
  father: {
    'GSM4472410_embryoD_TE.half.csv': {
      head: 'FATHER + SECOND PARENT', father: 'carries', sex: 'XX', counted: [], once: [], seen: [],
    },
    'GSM4472412_embryoE_TE.half.csv': {
      head: 'FATHER + SECOND PARENT', father: 'carries', sex: 'XY',
      counted: ['chr6 64.6-170.9Mb loss paternal | consistent'],
      once: ['chr6 66.4-170.9Mb loss maternal | one reading'],
      seen: ['chr6 whole loss unclear | also seen'],
    },
    'GSM4552364_embryoA_biopsy1.half.csv': {
      head: 'FATHER + SECOND PARENT', father: 'carries', sex: 'XX',
      counted: [
        'chr2 whole loss paternal | consistent',
        'chr7 whole loss maternal | consistent',
        'chr9 whole loss maternal | consistent',
        'chr11 whole loss paternal | consistent',
        'chr15 whole loss paternal | consistent',
        'chr16 47.8-87.7Mb loss paternal | consistent',
        'chr17 whole loss maternal | consistent',
        'chr19 whole loss paternal | consistent',
      ],
      once: [
        'chr4 184.6-190.8Mb loss paternal | one reading',
        'chr5 151.0-175.3Mb loss maternal | one reading',
        'chr6 112.0-170.9Mb loss maternal | one reading',
      ],
      seen: [],
    },
  },
  none: {
    'GSM5651119_z14_biopsy.half.csv': {
      head: 'TROPHECTODERM', sex: 'XY',
      counted: ['chr16 46.8-88.5Mb loss unclear | position only'], once: [], seen: [],
    },
    'GSM4472412_embryoE_TE.half.csv': {
      head: 'BULK', sex: 'XY',
      counted: ['chr6 66.4-170.9Mb loss unclear | position only'], once: [], seen: [],
    },
    'GSM4472410_embryoD_TE.half.csv': {
      head: 'TROPHECTODERM', sex: 'XX', counted: [], once: [], seen: [],
    },
  },
}

// --- 2. every set runs, and every card prints what its line says --------------------------------
for (const set of EXAMPLE_SETS) {
  const run = await runExamples(set.examples.map((e) => ({ path: `${DIR}${e.file}`, role: e.role })))
  const samples = set.examples.filter((e) => e.role === 'sample')
  const cards = run.cards.map(summary)
  assert.equal(cards.length, samples.length, `${set.id}: ${cards.length} cards for ${samples.length} samples`)

  // Two files labelled sperm: the run has to name the father, and name the right one.
  const candidates = set.examples.filter((e) => e.role === 'donor')
  if (candidates.length > 1) {
    assert.equal(run.father?.named, candidates[0].file.replace(/\.gz$/, ''),
      `${set.id}: the search says ${run.father?.line}`)
    const other = run.father!.rows.find((r) => r.name.startsWith('GSE186441_'))!
    assert.equal(other.carries, 0, `${set.id}: the second man carries in ${other.carries}`)
  }

  samples.forEach((e, i) => {
    const c = cards[i]
    const pin = PIN[set.id]?.[e.file.replace(/\.gz$/, '')]
    assert.ok(pin, `${set.id}/${e.gsm}: no pin for this example`)
    const where = `${set.id}/${e.gsm}`
    assert.equal(c.head, pin.head, `${where}: the headline reads ${c.head}`)
    // 'unread' is what the check reports for a parent that was not loaded at all.
    assert.equal(c.father, pin.father ?? 'unread', `${where}: the father's check reads ${c.father}`)
    assert.equal(c.mother, pin.mother, `${where}: the mother's check reads ${c.mother}`)
    assert.equal(c.sex, pin.sex, `${where}: the sex chromosomes read ${c.sex}`)
    assert.ok(c.decided, `${where}: the headline is not a result`)
    // A card that warns is not an example: the sets exist to show the tool reading, not refusing.
    assert.deepEqual(c.alerts, [], `${where}: ${c.alerts.join(' | ')}`)
    assert.deepEqual(c.counted, pin.counted, `${where}: counted ${c.counted.join('; ')}`)
    assert.deepEqual(c.once, pin.once, `${where}: seen once ${c.once.join('; ')}`)
    assert.deepEqual(c.seen, pin.seen, `${where}: read, not counted: ${c.seen.join('; ')}`)

    // The line starts with the headline the card prints, carries its count of changes, and names
    // no chromosome the card does not.
    assert.ok(e.expect.startsWith(pin.head), `${where}: line "${e.expect}" but the card reads ${pin.head}`)
    const n = pin.counted.length
    const said = /, (\d+|no) changes?[:,]/.exec(e.expect)
    assert.ok(said, `${where}: the line states no count of changes`)
    assert.equal(said[1] === 'no' ? 0 : Number(said[1]), n, `${where}: the card counts ${n}`)
    for (const m of e.expect.matchAll(/chromosome (\d+)/g)) {
      assert.ok(pin.counted.some((r) => r.startsWith(`chr${m[1]} `)),
        `${where}: the line names chromosome ${m[1]}, which the card does not`)
    }
  })
}

console.log(`examples.check.ts: ${EXAMPLE_SETS.length} example runs, `
  + `${EXAMPLE_FILES.length} files, every card prints what its line says`)
