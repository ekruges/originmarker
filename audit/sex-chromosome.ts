/**
 * What chrX and chrY actually read on arrays of stated sex-chromosome constitution.
 *
 * The tool calls sex from chrX heterozygosity alone, so one X reads male whether the sample is
 * 46,XY or 45,X. Separating those needs a second measurement, and this is where the numbers for
 * it come from rather than from a textbook expectation: the CytoScan miscarriage series states a
 * karyotype per sample, including 45,X, so the two constitutions can be measured side by side on
 * one platform.
 *
 * Measures per array, against the array's own autosomes:
 *   xShift    median log2 ratio on chrX minus the autosomal median. One copy against two.
 *   xHet      chrX heterozygosity over autosomal heterozygosity. What callSex reads today.
 *   yShift    median log2 ratio on chrY minus the autosomal median.
 *   yCall     chrY call rate over the autosomal call rate.
 *
 * Run: OM_PROBES=<dir of .probes> OM_LABELS=<series labels json> \
 *        node --experimental-strip-types audit/sex-chromosome.ts
 */
import { createReadStream, existsSync, readFileSync, readdirSync } from 'node:fs'
import { createInterface } from 'node:readline'
import { join } from 'node:path'

const DIR = process.env.OM_PROBES ?? ''
const LABELS = process.env.OM_LABELS ?? ''
if (!DIR || !existsSync(DIR) || !LABELS || !existsSync(LABELS)) {
  console.log('sex-chromosome: needs OM_PROBES and OM_LABELS. NOT RUN.')
  process.exit(0)
}
const labels = new Map<string, string>(
  (JSON.parse(readFileSync(LABELS, 'utf8')) as { gsm: string; label: string }[])
    .map((r) => [r.gsm, r.label]),
)

const median = (xs: number[]): number =>
  (xs.length ? [...xs].sort((a, b) => a - b)[xs.length >> 1] : NaN)

/** Every twentieth autosomal value fixes the autosomal median far inside the gap this decides. */
const STRIDE = 20

interface Row {
  gsm: string; label: string
  xShift: number; xHet: number; yShift: number; yCall: number; markers: number
}

async function measure(path: string, gsm: string): Promise<Row> {
  const x: number[] = []; const y: number[] = []; const auto: number[] = []
  let xCalled = 0; let xHet = 0; let autoCalled = 0; let autoHet = 0
  let yCalled = 0; let yTotal = 0; let autoTotal = 0; let seen = 0
  const rl = createInterface({ input: createReadStream(path), crlfDelay: Infinity })
  let head: string[] | null = null
  let iChr = 1; let iLrr = 3; let iGt = 6
  for await (const line of rl) {
    const f = line.split('\t')
    if (!head) {
      head = f
      iChr = f.indexOf('chr'); iLrr = f.indexOf('log2R'); iGt = f.indexOf('genotype')
      continue
    }
    const c = f[iChr]
    const lrr = Number(f[iLrr])
    const gt = f[iGt]
    const called = gt === 'AA' || gt === 'AB' || gt === 'BB'
    if (c === 'X' || c === '23') {
      if (Number.isFinite(lrr)) x.push(lrr)
      if (called) { xCalled += 1; if (gt === 'AB') xHet += 1 }
    } else if (c === 'Y' || c === '24') {
      yTotal += 1
      if (called) yCalled += 1
      if (Number.isFinite(lrr)) y.push(lrr)
    } else if (/^\d+$/.test(c)) {
      autoTotal += 1
      if (called) { autoCalled += 1; if (gt === 'AB') autoHet += 1 }
      seen += 1
      if (seen % STRIDE === 0 && Number.isFinite(lrr)) auto.push(lrr)
    }
  }
  const autoRate = autoTotal ? autoCalled / autoTotal : NaN
  return {
    gsm,
    label: labels.get(gsm) ?? '?',
    xShift: median(x) - median(auto),
    xHet: (xCalled ? xHet / xCalled : NaN) / (autoCalled ? autoHet / autoCalled : NaN),
    yShift: median(y) - median(auto),
    yCall: (yTotal ? yCalled / yTotal : NaN) / autoRate,
    markers: autoTotal + xCalled + yTotal,
  }
}

const files = readdirSync(DIR).filter((f) => f.endsWith('.probes')).sort()
const cap = Number(process.env.OM_CAP ?? 9999)
console.log(`arrays: ${Math.min(files.length, cap)} of ${files.length}`)
const rows: Row[] = []
for (const f of files.slice(0, cap)) {
  rows.push(await measure(join(DIR, f), f.replace(/\..*$/, '')))
  if (rows.length % 25 === 0) console.log(`  measured ${rows.length}`)
}

const q = (xs: number[], p: number) => {
  const s = [...xs].filter(Number.isFinite).sort((a, b) => a - b)
  return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN
}
const f2 = (v: number) => (Number.isFinite(v) ? v.toFixed(3) : '-').padStart(7)

/** The stated constitutions worth separating, plus the normal arrays split by what chrX reads. */
const classOf = (r: Row): string => {
  if (/^45,X$/.test(r.label)) return '45,X (stated)'
  if (/Triploidy/i.test(r.label)) return 'triploidy (stated)'
  if (/No clinical Variations/i.test(r.label)) {
    return r.xHet < 0.15 ? 'normal, one X by het' : 'normal, two X by het'
  }
  return 'other stated abnormality'
}
const groups = new Map<string, Row[]>()
for (const r of rows) groups.set(classOf(r), [...(groups.get(classOf(r)) ?? []), r])

console.log('')
console.log('class                       n   xShift median  (10th-90th)      xHet   yShift    yCall')
for (const [k, g] of [...groups].sort()) {
  const xs = g.map((r) => r.xShift)
  console.log(`${k.padEnd(26)}${String(g.length).padStart(3)}  ${f2(q(xs, 0.5))}   `
    + `${f2(q(xs, 0.1))} ${f2(q(xs, 0.9))}  ${f2(q(g.map((r) => r.xHet), 0.5))} `
    + `${f2(q(g.map((r) => r.yShift), 0.5))} ${f2(q(g.map((r) => r.yCall), 0.5))}`)
}

// The question the rule turns on: does a stated 45,X separate from a normal male on intensity,
// on chrY, or on neither?
const one = groups.get('45,X (stated)') ?? []
const male = (groups.get('normal, one X by het') ?? [])
console.log('')
console.log(`45,X arrays: ${one.length}, normal arrays reading one X on heterozygosity: `
  + `${male.length}`)
if (one.length && male.length) {
  const sep = (k: keyof Row) => {
    const a = one.map((r) => r[k] as number).filter(Number.isFinite)
    const b = male.map((r) => r[k] as number).filter(Number.isFinite)
    return `${k}: 45,X ${f2(q(a, 0.1))} to ${f2(q(a, 0.9))}   `
      + `one-X normal ${f2(q(b, 0.1))} to ${f2(q(b, 0.9))}`
  }
  console.log(`  ${sep('xShift')}`)
  console.log(`  ${sep('yCall')}`)
  console.log(`  ${sep('yShift')}`)
}
console.log('')
for (const r of rows.filter((x) => /45,X/.test(x.label)).slice(0, 12)) {
  console.log(`  ${r.gsm}  xShift ${f2(r.xShift)}  xHet ${f2(r.xHet)}  `
    + `yShift ${f2(r.yShift)}  yCall ${f2(r.yCall)}`)
}
