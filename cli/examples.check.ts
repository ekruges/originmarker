// The bundled examples still produce what they claim. Run: node cli/examples.check.ts
//
// This is the only check in the repo that guards an ANSWER rather than a build. Every example
// carries an `expect` string making specific numeric claims: an opposite-homozygote rate, a second
// parental contribution, an inferred dropout. Those are prose, nothing verified them, and every
// threshold moved in this project could falsify one silently. A user's first contact with the tool
// is clicking one of these, so an example whose text no longer matches its output is the worst
// place for the numbers to be wrong.
//
// It runs the real modules over the real bundled files. No fixtures: the point is that the file a
// user downloads still behaves as advertised.
import assert from 'node:assert/strict'
import { readFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { EXAMPLES, EXAMPLE_MARKERS } from '../web/src/examples.ts'
import {
  headerMap, parseRow, emptyBafSums, accumulateBaf, accumulate, finishProfile,
} from '../web/src/ingest.ts'
import { isAutosome } from '../web/src/parentage.ts'
import { inferStage } from '../web/src/stage.ts'
import { emptyHet, addOneParent, hetCall } from '../web/src/obligateHet.ts'
import type { AB } from '../web/src/informativity.ts'

const DIR = new URL('../web/public/examples/', import.meta.url).pathname

interface Loaded { gt: Map<string, AB>; markers: number; callRate: number; hetRate: number }

function load(file: string): Loaded {
  const path = `${DIR}${file}`
  assert.ok(existsSync(path), `${file} is listed in examples.ts but not present in public/examples`)
  const raw = readFileSync(path)
  const text = (file.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8')
  const lines = text.split('\n')
  let h = -1
  for (let i = 0; i < 60; i += 1) if (lines[i] && !lines[i].startsWith('#')) { h = i; break }
  const map = headerMap(lines[h])
  assert.ok(map, `${file}: header not recognised by the shipped ingest`)
  const gt = new Map<string, AB>()
  const byChrom = new Map()
  const baf = emptyBafSums()
  let first = ''
  let markers = 0
  for (let i = h + 1; i < lines.length; i += 1) {
    const r = parseRow(lines[i], map as never)
    if (!r) continue
    markers += 1
    if (!first) first = r.probesetId
    accumulate(r as never, byChrom as never)
    accumulateBaf(r as never, baf as never)
    if (isAutosome(r.chrom) && r.genotype !== 'NC') gt.set(r.probesetId, r.genotype as AB)
  }
  const p = finishProfile(file, byChrom as never, baf as never, first)
  return { gt, markers, callRate: p.callRate, hetRate: p.hetRate }
}

const oppositeHom = (a: Map<string, AB>, b: Map<string, AB>) => {
  let n = 0
  let o = 0
  for (const [probe, ga] of a) {
    if (ga !== 'AA' && ga !== 'BB') continue
    const gb = b.get(probe)
    if (gb !== 'AA' && gb !== 'BB') continue
    n += 1
    if (ga !== gb) o += 1
  }
  return n ? o / n : NaN
}

const secondParent = (ref: Map<string, AB>, s: Map<string, AB>) => {
  const t = emptyHet()
  for (const [probe, pg] of ref) {
    const cg = s.get(probe)
    if (cg) addOneParent(pg as never, cg as never, t as never)
  }
  return hetCall(t as never, 1) as { ploidy: string; fraction: number }
}

// --- 1. every listed file exists, parses, and is the advertised size ----------------------------
const loaded = new Map<string, Loaded>()
for (const e of EXAMPLES) {
  const l = load(e.file)
  loaded.set(e.gsm, l)
  assert.equal(l.markers, EXAMPLE_MARKERS,
    `${e.gsm}: ${l.markers} markers, but examples.ts advertises ${EXAMPLE_MARKERS}`)
  assert.ok(l.gt.size > 50_000, `${e.gsm}: only ${l.gt.size} autosomal calls, which is not an array`)
}

// --- 2. the numeric claims in each `expect` string still hold ------------------------------------
//
// Only claims that can be checked from the files themselves are asserted. A claim about a rate is
// pulled out of the prose by pattern, so adding a number to an `expect` string automatically puts
// it under test rather than needing a second edit here.
const pct = (s: string, label: RegExp): number | null => {
  const m = label.exec(s)
  if (!m) return null
  // The first DEFINED group, not group 1: an alternation captures into different slots, and
  // reading group 1 blindly yields undefined and then NaN, which compares false against every
  // bound and so passes silently. That is the failure mode this whole file exists to prevent.
  const hit = m.slice(1).find((g) => g !== undefined)
  const v = hit === undefined ? NaN : Number(hit)
  return Number.isFinite(v) ? v : null
}

/** Pairs the run is expected to form: reference GSM -> sample GSMs it is the parent of. */
const FAMILIES: [string, string[]][] = [
  ['GSM4472397', ['GSM4472424']],
  ['GSM8826446', ['GSM8826445', 'GSM8826436']],
]

for (const [ref, kids] of FAMILIES) {
  const r = loaded.get(ref)
  if (!r) continue
  for (const kid of kids) {
    const k = loaded.get(kid)
    if (!k) continue
    const e = EXAMPLES.find((x) => x.gsm === kid)!
    const opp = oppositeHom(r.gt, k.gt)
    const link = secondParent(r.gt, k.gt)

    assert.ok(opp < 0.03,
      `${kid}: opposite-homozygote ${opp.toFixed(4)} against ${ref} is not a parent-child rate`)
    assert.equal(link.ploidy, 'biparental',
      `${kid}: no longer reads as a child of ${ref}, which its text claims`)

    // "a 1.11% opposite-homozygote rate"
    const claimedOpp = pct(e.expect, /([\d.]+)%\s+opposite-homozygote/)
    if (claimedOpp !== null) {
      assert.ok(Math.abs(100 * opp - claimedOpp) < 0.25,
        `${kid}: text claims ${claimedOpp}% opposite-homozygote, measured ${(100 * opp).toFixed(2)}%`)
    }
    // "second parental contribution at 9.0%" or "a 8.7% second contribution"
    const claimedSecond = pct(e.expect, /([\d.]+)%\s+of\s+[\d,]+\s+informative/)
      ?? pct(e.expect, /contribution at ([\d.]+)%/)
      ?? pct(e.expect, /([\d.]+)%\s+second contribution/)
    if (claimedSecond !== null) {
      assert.ok(Math.abs(100 * link.fraction - claimedSecond) < 0.4,
        `${kid}: text claims ${claimedSecond}% second contribution, measured `
        + `${(100 * link.fraction).toFixed(1)}%`)
    }
    // "the inferred dropout is 0.193"
    const m = /dropout (?:is |of )?(0\.\d+)/.exec(e.expect)
    if (m) {
      const got = inferStage({ hetRate: k.hetRate, callRate: k.callRate } as never)
      assert.ok(Math.abs(got.dropout - Number(m[1])) < 0.01,
        `${kid}: text claims dropout ${m[1]}, the stage module now infers ${got.dropout.toFixed(3)}`)
    }
  }
}

// --- 3. the unrelated examples are still unrelated ------------------------------------------------
//
// Half the set exists to show the tool saying no. If one of these quietly started reading as a
// child, the example would teach the opposite of its text.
for (const gsm of ['GSM4472407', 'GSM4472415', 'GSM8826446']) {
  const s = loaded.get(gsm)
  const r = loaded.get('GSM4472397')
  if (!s || !r) continue
  assert.ok(oppositeHom(r.gt, s.gt) > 0.03,
    `${gsm} now reads related to the sperm donor, but its text says unrelated`)
}

console.log(`examples.check.ts: ${EXAMPLES.length} examples load, parse and still produce the `
  + 'rates their own text claims')
