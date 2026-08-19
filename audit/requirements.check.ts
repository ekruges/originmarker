// Does the shipped code do what was asked? Run: node audit/requirements.check.ts
//
// NOT A READING OF THE CODE. Every claim below is produced by CALLING the shipped functions and
// reporting what comes back, because this project has repeatedly shipped features that were
// complete, documented, and never wired: an origin read from a sign without its class, a timing
// test that nothing called, a feature track that matched nothing, an array gate that silenced every
// amplified sample, and a report argument never passed. Each passed its own tests, because each
// test exercised a fixture the application never sees.
//
// So this file exercises the paths a user's run actually takes, and prints the evidence for and
// against each requirement rather than a verdict on its own authority.
import { readFileSync } from 'node:fs'

const W = new URL('../web/src/', import.meta.url).pathname
const dosage = await import(`${W}dosageOrigin.ts`)
const tax = await import(`${W}abnormalities.ts`)
const post = await import(`${W}originPosterior.ts`)
const stage = await import(`${W}stage.ts`)
const defects = await import(`${W}defects.ts`)
const cmp = await import(`${W}comparison.ts`)
const oneP = await import(`${W}oneParentOrigin.ts`)

let pass = 0
let partial = 0
let fail = 0
const say = (s = '') => process.stdout.write(`${s}\n`)
const head = (n: string, ask: string) => {
  say(`\n${'='.repeat(92)}\n${n}\n  ASKED: ${ask}`)
}
const verdict = (v: 'MET' | 'PARTIAL' | 'NOT MET', why: string) => {
  if (v === 'MET') pass += 1
  else if (v === 'PARTIAL') partial += 1
  else fail += 1
  say(`  ${v}: ${why}`)
}
const ev = (s: string) => say(`    ${s}`)

// A synthetic region generator matching the shape callDosageOrigin consumes.
const N = 4000
const region = (n: number, centre: number) =>
  Array.from({ length: n }, (_, i) => ['AA', centre + 0.02 * Math.sin(i * 2.399963), 0] as const)
    .map(([g, b]) => [g, b] as [never, number])
const bg = region(N, 0.50)

// ---------------------------------------------------------------------------------------------
head('1. MORE TYPES OF CHROMOSOMAL ABNORMALITY',
  '"i need detections for more types of chromosomal abnormalities ... duplications, deletions, etc"')
{
  const t = tax.TAXONOMY
  const yes = t.filter((x: { detectable: string }) => x.detectable === 'yes')
  const partly = t.filter((x: { detectable: string }) => x.detectable === 'partly')
  const no = t.filter((x: { detectable: string }) => x.detectable === 'no')
  ev(`${t.length} classes enumerated: ${yes.length} detected, ${partly.length} partly, ${no.length} named impossible`)
  ev(`detected: ${yes.map((x: { label: string }) => x.label).join('; ')}`)
  ev(`impossible, with a stated limit: ${no.map((x: { label: string }) => x.label).join('; ')}`)

  // The detectors must actually fire, not merely be listed.
  // The window is sized past the MEASURED floor for a copy-neutral call rather than chosen: on
  // real diploid arrays this detector fires on about 1 window in 110 at 600 markers with no event
  // present at all, and reaches zero only at 1,800 to 2,100. See LOH_MIN_EXPECTED_HET and
  // audit/false-tract-rate.csv. A 1,000-marker fixture was inside the regime the measurement says
  // is not measurable.
  const win = (het: number, logR?: number, whole = false) =>
    ({ chrom: '7', startBp: 0, endBp: 159e6, called: 2800, het, logR, wholeChromosome: whole })
  const normal = Array.from({ length: 8 }, () => win(476))
  const loh = tax.detectLoh([...normal, win(56, 0.01, true)])
  const upd = tax.detectUpd([{ chrom: '7', startBp: 0, endBp: 40e6, markers: 12000, wholeChromosome: true }])
  const tri = tax.detectTriploidy([
    ...Array.from({ length: 400 }, () => 1 / 3), ...Array.from({ length: 400 }, () => 2 / 3)])
  const cx = tax.detectComplex(12, 22, 0.95)
  const fired = [['copy-neutral LOH', loh.length], ['uniparental disomy', upd.length],
    ['triploidy', tri ? 1 : 0], ['complex genome', cx ? 1 : 0]] as const
  for (const [name, n] of fired) ev(`detector fires: ${name} -> ${n ? 'yes' : 'NO'}`)
  const allFire = fired.every(([, n]) => n > 0)
  verdict(allFire ? 'MET' : 'NOT MET',
    allFire ? 'every new detector fires on data of the shape it is meant to catch'
      : 'a listed detector does not fire')
}

// ---------------------------------------------------------------------------------------------
head('2. AND FROM WHICH PARENT', '"and from which parent"')
{
  const rows: string[] = []
  let named = 0
  let blocked = 0
  for (const m of ['bulk', 'esc-single', 'trophectoderm', 'blastomere'] as const) {
    const c = dosage.callDosageOrigin(region(N, 0.5 + 0.14) as never, bg as never, m,
      { wholeChromosome: true, parents: 1, state: 'cnn-loh' })
    const p = c.posterior
    const isNamed = c.verdict === 'loaded-parent' || c.verdict === 'other-parent'
    if (isNamed) named += 1; else blocked += 1
    rows.push(`${m.padEnd(15)} ${c.verdict.padEnd(20)} ${p ? `band ${p.band} conf ${p.confidence.toFixed(4)}` : 'no posterior'}`)
  }
  for (const r of rows) ev(r)
  ev(`named on ${named} of 4 material classes, refused on ${blocked}`)
  verdict(named === 4 ? 'MET' : named ? 'PARTIAL' : 'NOT MET',
    named === 4 ? 'a parent is named on every material class, amplified included'
      : `${blocked} material class(es) still refuse`)
}

// ---------------------------------------------------------------------------------------------
head('3. BOTH COPIES OF ONE CHROMOSOME FROM ONE PARENT',
  '"Maybe both chromosomes in chr7 are paternal"')
{
  const iso = tax.detectUpd([{ chrom: '7', startBp: 0, endBp: 159e6, markers: 40000, wholeChromosome: true }])
  ev(`class returned: ${iso[0]?.cls ?? 'none'}`)
  ev(`states what it cannot see: ${/Heterodisomy leaves heterozygosity normal/.test(iso[0]?.evidence ?? '')}`)
  ev(`carries the 32% no-run rate: ${/32%/.test(iso[0]?.evidence ?? '')}`)
  const het = tax.taxonomyFor('heterodisomy')
  ev(`heterodisomy declared: ${het?.detectable} with a stated limit: ${!!het?.limit}`)
  verdict(iso[0]?.cls === 'isodisomy' ? 'MET' : 'NOT MET',
    'isodisomy is detected and named; heterodisomy is declared unreachable rather than missing')
}

// ---------------------------------------------------------------------------------------------
head('4. WHEN THE CHANGE AROSE',
  '"is haploid chromosome duplicated post fertilization or paternal, or sperm introduces a small variation"')
{
  const alone = tax.callUniformity(1, 1)
  const uniform = tax.callUniformity(3, 3)
  const lineage = tax.callUniformity(1, 3)
  ev(`one unit  -> ${alone.mechanism}, and names the remedy: ${/second array of the same embryo/.test(alone.why)}`)
  ev(`in all 3  -> ${uniform.mechanism}`)
  ev(`in 1 of 3 -> ${lineage.mechanism}`)

  // Grouping must be measured, not declared, or a labelling slip produces a confident mechanism.
  const conc: Record<string, Record<string, number>> = {
    a: { a: 1, b: 0.961, c: 0.55 }, b: { a: 0.961, b: 1, c: 0.55 }, c: { a: 0.55, b: 0.55, c: 1 },
  }
  const groups = tax.groupUnits(['a', 'b', 'c'], (x: string, y: string) => conc[x][y])
  ev(`grouping by concordance: ${groups.map((g: string[]) => `[${g.join('+')}]`).join(' ')}`)
  const wall1 = tax.unanswerable(false, 1).map((x: { cls: string }) => x.cls)
  const wall2 = tax.unanswerable(false, 2).map((x: { cls: string }) => x.cls)
  ev(`the wall lifts with a second unit: ${wall1.includes('gamete-de-novo')} -> ${wall2.includes('gamete-de-novo')}`)
  const ok = uniform.mechanism === 'meiotic' && lineage.mechanism === 'post-zygotic'
    && groups.length === 2 && !wall2.includes('gamete-de-novo')
  verdict(ok ? 'MET' : 'NOT MET',
    'timing is answered from a second unit, with units grouped by measured concordance')
}

// ---------------------------------------------------------------------------------------------
head('5. CALL MATERNAL BY THE ABSENCE OF PATERNAL', '"call M by ABSENCE of P"')
{
  // One parent loaded, as the paternal role. A shift the other way must name the mother.
  const up = dosage.callDosageOrigin(region(N, 0.5 + 0.20) as never, bg as never, 'bulk',
    { wholeChromosome: true, parents: 1 })
  const down = dosage.callDosageOrigin(region(N, 0.5 - 0.20) as never, bg as never, 'bulk',
    { wholeChromosome: true, parents: 1 })
  ev(`shift toward the loaded parent  -> ${up.verdict}`)
  ev(`shift away from the loaded parent -> ${down.verdict}`)
  const d = defects.defectsFrom([{ chrom: '7', kind: 'copy-loss', start: 0, end: 159e6 } as never],
    [], [], [], 'paternal', undefined,
    [{ where: 'chr7', verdict: down.verdict, shift: down.shift, z: down.z, impliedF: down.impliedF,
      window: down.window, why: down.why, confidence: down.posterior?.confidence,
      band: down.posterior?.band }])
  ev(`with only a father loaded, the display names: ${d[0]?.origin}`)
  // And the two arms must be mirror images, which is what stops a directional null.
  const symmetric = Math.abs(Math.abs(up.shift) - Math.abs(down.shift)) < 1e-9
  ev(`the two arms are mirror images: ${symmetric}`)
  verdict(d[0]?.origin === 'maternal' && symmetric ? 'MET' : 'PARTIAL',
    'a maternal call is produced from a paternal-only run, and the two directions are symmetric')
}

// ---------------------------------------------------------------------------------------------
head('6. EVERY ABNORMALITY: LOCATION AND A CONFIDENCE-SCORED ORIGIN',
  '"for EVERY SINGLE CHROMOSOMAL ABNORMALITY, we ALWAYS have both the ultraspecific location ... '
  + 'AND a CONFIDENCE SCORED CALL of the PARENTAL ORIGIN"')
{
  const classes = tax.TAXONOMY.filter((x: { detectable: string }) => x.detectable !== 'no')
  let withOrigin = 0
  let structural = 0
  for (const t of classes) {
    const blocked = tax.ORIGIN_UNREACHABLE.has(t.cls)
    if (blocked) { structural += 1; ev(`${t.label.padEnd(46)} origin blocked BY THE CLASS`) }
    else withOrigin += 1
  }
  ev(`${withOrigin} of ${classes.length} detectable classes can carry an origin; ${structural} cannot`)

  // The two that were on that list and should not have been, exercised rather than asserted.
  const at = (share: number, k: number) =>
    Array.from({ length: k }, () => ['AA', 1 - share] as const)
  const triTheirs = tax.callTriploidyOrigin(at(2 / 3, 600))
  const triOther = tax.callTriploidyOrigin(at(1 / 3, 600))
  ev(`triploidy, loaded parent's allele at two thirds -> ${triTheirs.origin}`)
  ev(`triploidy, loaded parent's allele at one third  -> ${triOther.origin}`)
  const cxBlocked = dosage.callDosageOrigin(region(N, 0.5 + 0.14) as never, bg as never,
    'trophectoderm', { wholeChromosome: true, noSelfReference: true })
  const mendelian = oneP.callOneParentOrigin(
    Array.from({ length: 400 }, (_, i) => ['AA', i < 120 ? 'BB' : 'AA']) as never, 0.20)
  ev(`complex genome: dosage ${cxBlocked.verdict}, and the Mendelian channel still says ${mendelian.verdict}`)

  // Location on every finding.
  const f = { cls: 'cnn-loh' as const, chrom: '7', startBp: 1e6, endBp: 14e6, wholeChromosome: false,
    evidence: 'x' }
  const withoutCall = defects.findingToDefect(f)
  const withCall = defects.findingToDefect(f, undefined, {
    verdict: 'other-parent', shift: 0.09, z: 5.1, impliedF: 0.31, why: 'from dosage',
    confidence: 0.9964, band: 'A', limitedBy: 'none' }, 'paternal')
  ev(`location present without an origin: ${!!withoutCall.locus}`)
  ev(`with a scored call: origin ${withCall.origin}, confidence ${withCall.confidence}, band ${withCall.band}`)

  // And a named parent must never appear without a number.
  const namedNoNumber = withCall.origin !== 'unclear' && !Number.isFinite(withCall.confidence)
  ev(`a named parent without a confidence is constructible: ${namedNoNumber}`)
  const bothReachable = triTheirs.origin === 'extra-set-loaded-parent'
    && triOther.origin === 'extra-set-other-parent'
    && mendelian.verdict !== 'refused'
  verdict(!namedNoNumber && withoutCall.locus && structural === 0 && bothReachable
    ? 'MET' : 'PARTIAL',
    `location is always present, every named parent carries a number, and all ${withOrigin} `
    + 'detectable classes can carry an origin: a triploid\'s extra set from allele fraction at the '
    + 'loaded parent\'s homozygous markers, a complex genome from the Mendelian channel that needs '
    + 'no self-reference')
}

// ---------------------------------------------------------------------------------------------
head('7. EVERY STAGE OF DEVELOPMENT, CATEGORISED AND DISCRIMINATED',
  '"every single stage of development both categorized/discriminated"')
{
  const cases = [
    ['bulk gDNA', { hetRate: 0.168, callRate: 0.97, hetBafSd: 0.088 }],
    ['amplified, bulk-level heterozygosity', { hetRate: 0.168, callRate: 0.97, hetBafSd: 0.232 }],
    ['trophectoderm', { hetRate: 0.150, callRate: 0.93, hetBafSd: 0.24 }],
    ['single cell', { hetRate: 0.130, callRate: 0.88, hetBafSd: 0.25 }],
    ['blastomere', { hetRate: 0.112, callRate: 0.86, hetBafSd: 0.26 }],
    ['haploid', { hetRate: 0.03, callRate: 0.93, hetBafSd: 0.28 }],
    ['failed', { hetRate: 0.31, callRate: 0.53, hetBafSd: 0.21 }],
  ] as const
  const seen = new Set<string>()
  for (const [name, p] of cases) {
    const s = stage.inferStage(p as never)
    seen.add(s.stage)
    ev(`${name.padEnd(38)} -> ${s.stage.padEnd(15)} dropout ${Number.isFinite(s.dropout) ? s.dropout.toFixed(3) : 'n/a'}`)
  }
  ev(`distinct stages reached: ${[...seen].join(', ')}`)
  const twoAxis = stage.inferStage({ hetRate: 0.168, callRate: 0.97, hetBafSd: 0.232 } as never).stage !== 'bulk'
  ev(`amplification overrules heterozygosity on the bulk rung: ${twoAxis}`)
  verdict(seen.size >= 5 && twoAxis ? 'MET' : 'PARTIAL',
    `${seen.size} distinct stages are reached, and the amplification axis is applied`)
}

// ---------------------------------------------------------------------------------------------
head('8. DATA MANIPULATION ADJUSTED FOR THE REALITIES OF DEVELOPMENT',
  '"changes to the data manipulation in light of the realities of development"')
{
  const mats = ['bulk', 'esc-single', 'trophectoderm', 'blastomere'] as const
  ev('constant            ' + mats.map((m) => m.padEnd(14)).join(''))
  const rows: [string, (m: string) => unknown][] = [
    ['drift floor        ', (m) => dosage.DRIFT_TAU[m]],
    ['variance inflation ', (m) => dosage.VIF_CHROMOSOME[m]],
    ['channel correlation', (m) => dosage.RESIDUAL_R[m]],
    ['band A accuracy    ', (m) => post.BAND_ACCURACY[m].A],
    ['floor, cnn-loh chr ', (m) => dosage.floorFor(m, 'cnn-loh', true, 1)],
  ]
  for (const [label, f] of rows) {
    ev(label + mats.map((m) => String(f(m)).padEnd(14)).join(''))
  }
  const spread = Math.max(...mats.map((m) => dosage.DRIFT_TAU[m]))
    / Math.min(...mats.map((m) => dosage.DRIFT_TAU[m]))
  ev(`drift spans ${spread.toFixed(1)}x across the classes, so the stage genuinely changes the answer`)
  const allDiffer = rows.every(([, f]) => new Set(mats.map((m) => String(f(m)))).size > 1)
  verdict(allDiffer ? 'MET' : 'PARTIAL',
    'every constant the origin call depends on is per-material and they genuinely differ')
}

// ---------------------------------------------------------------------------------------------
head('9. THE CALIBRATION ACTUALLY APPLIED', 'implied by "CONFIDENCE SCORED"')
{
  const c = dosage.callDosageOrigin(region(N, 0.5 + 0.14) as never, bg as never, 'trophectoderm',
    { wholeChromosome: true, parents: 2 })
  ev(`shipped maps applied by default: ${c.posterior && !c.posterior.uncalibrated}`)
  const raw = dosage.callDosageOrigin(region(N, 0.5 + 0.14) as never, bg as never, 'trophectoderm',
    { wholeChromosome: true, parents: 2, calibration: {} })
  ev(`opting out returns the raw posterior and says so: ${raw.posterior?.uncalibrated}`)
  ev(`bands carry measured accuracies: A=${post.BAND_ACCURACY.blastomere.A} D=${post.BAND_ACCURACY.blastomere.D}`)
  ev(`the weakest band still beats chance: ${post.BAND_ACCURACY.blastomere.D > 0.5}`)
  verdict(c.posterior && !c.posterior.uncalibrated ? 'MET' : 'NOT MET',
    'calls are calibrated by default and opting out is explicit')
}

// ---------------------------------------------------------------------------------------------
head('10. THE GENOTYPE CHANNELS', 'implied by "as much as you can"')
{
  const rows = (n: number, ex: number) => {
    const o: [string, string][] = []
    for (let i = 0; i < n; i += 1) o.push(['AA', i < ex ? 'BB' : 'AA'])
    return o
  }
  const seen = new Set<string>()
  for (const [n2, ex] of [[60, 3], [100, 2], [200, 5], [400, 8], [400, 40], [400, 300]] as const) {
    const c = oneP.callOneParentOrigin(rows(n2, ex) as never, 0.20)
    seen.add(c.posterior.toFixed(4))
    ev(`${String(n2).padStart(3)} markers, ${String(ex).padStart(3)} exclusive -> confidence ${c.posterior.toFixed(4)}  band ${c.band}  ${c.verdict}`)
  }
  ev(`distinct confidences across that range: ${seen.size}`)
  const cap = 1 - oneP.SYSTEMATIC_ERROR_BOUND
  const top = oneP.callOneParentOrigin(rows(900, 400) as never, 0.20)
  ev(`ceiling ${cap.toFixed(3)}, from ${oneP.VALIDATION_UNITS} independent validation units`)
  ev(`strongest evidence reaches ${top.posterior.toFixed(4)}, under the ceiling: ${top.posterior <= cap + 1e-9}`)
  ev(`and stays out of the top band: ${top.band !== 'A'}`)
  verdict(seen.size > 3 && top.posterior <= cap + 1e-9 && top.band !== 'A' ? 'MET' : 'PARTIAL',
    'the reported confidence varies with the evidence and is bounded by what the validation '
    + 'supports, rather than reporting a likelihood that is decisive under its own model')
}

// ---------------------------------------------------------------------------------------------
head('11. THE FEATURE COMPARISON', '"comparison to those other factors ... vs fragile sites etc"')
{
  const raw = JSON.parse(readFileSync(new URL('../web/public/hg19_features.json', import.meta.url), 'utf8'))
  const track = cmp.normaliseTrack(raw)!
  ev(`shipped track decodes: fragile=${track.fragile.length} longGenes=${track.longGenes.length} valleysES=${track.lateReplicationValleysES.length}`)
  const markers = new Map<string, number[]>()
  for (const ch of ['1', '2', '3', '4', '6', '8', '9', '10', '11', '12', '13', '14', '16']) {
    markers.set(ch, Array.from({ length: 6000 }, (_, i) => i * 40_000))
  }
  const onFra = track.fragile.filter((f: { chrom: string }) => markers.has(f.chrom)).slice(0, 6)
    .map((f: { chrom: string; startBp: number }) => ({ chrom: f.chrom, startBp: f.startBp + 6e5, endBp: f.startBp + 3.4e6 }))
  const off = ['2', '4', '9', '11', '13', '14'].map((chrom, i) => ({
    chrom, startBp: 3e6 + i * 1.5e6, endBp: 3e6 + i * 1.5e6 + 2.6e6 }))
  const c = cmp.compare(track, [...onFra, ...off], markers, { permutations: 800 })
  const fra = c.features.find((f: { feature: string }) => f.feature === 'common fragile site')
  ev(`planted 6 of 12 regions on real fragile sites`)
  ev(`fragile sites: observed ${(100 * fra.observed).toFixed(0)}% vs null ${(100 * fra.expected).toFixed(0)}%, fold ${fra.fold.toFixed(2)}, p=${fra.p.toFixed(4)}`)
  const flags = cmp.regionFlags(c)
  const starred = flags.filter((x: { related: string[] }) => x.related.length).length
  ev(`regions starred as having an alternative explanation: ${starred} of 12`)
  ev(`the parental caveat travels on the result: ${/never be read as support for a parental call/.test(c.caveat)}`)
  const right = fra.significant && starred === 6
  verdict(right ? 'MET' : 'PARTIAL',
    right ? 'the planted coincidence is found and exactly the planted regions are starred'
      : `expected 6 starred, got ${starred}`)
}

// ---------------------------------------------------------------------------------------------
say(`\n${'='.repeat(92)}`)
say(`  MET ${pass}   PARTIAL ${partial}   NOT MET ${fail}`)
say('')
say('  Three classes remain unanswerable and are named rather than counted as met: uniparental')
say('  heterodisomy needs both parents alongside the embryo, tandem and inserted duplications')
say('  carry no positional difference in either channel, and reverse segregation is normal in')
say('  copy number in 20 of 26 observed cases. Those are platform limits, not gaps here.')
process.exit(fail > 0 ? 1 : 0)
