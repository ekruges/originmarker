/**
 * What the shipped scoreSample reports for each array of a list, written out whole, so every stated
 * fact in audit/local/truth can be scored against it offline.
 *
 * The list is tab-separated: array id, array path, parent path, tag, and optionally the other
 * parent's path, `maternal` where the first parent is the mother, and the paths of the other arrays
 * of the same run, separated by `;`. A parent path of `-` is no parent. Each record also carries
 * `stageFree`, the stage with no parent consulted.
 *
 * Every claim is written whether the tool printed it or withheld it. A withheld claim is in
 * `refusals` with what it would have said.
 *
 * A file whose markers are not the proven chip's is not scored: its record carries `chip` with the
 * refusal and the causes, as the page prints them, and nothing else. See web/src/platform.ts.
 *
 * `card` holds the result card's change list as the overview counts it: the counted rows with their
 * pieces, the count of rows seen once, and whether the genome is too disturbed to reference; the
 * card's own count of changes and findings; and, on a sample not reportable, the rows listed apart.
 *
 * `surfaces` holds what the other surfaces print of the same sample: the overview's count of
 * changes and findings for it alone, the changes download's words and loaded-parent values, and,
 * on a sample not reportable or one whose sex line reads "47,XXY or a mix of cells", what its PDF
 * prints (pageSays, overview, syngamyTable, syngamyPdf); on a sample not reportable the run log's
 * gain, loss and stretch lines (`log`); and, where a whole call was withdrawn after the chromosome
 * was read, the note the chromosome table prints for it (`wholeWithdrawn`).
 *
 * `genomeVerdicts` holds the loaded parent's genome-wide verdict and, with both parents, the
 * mother's, read as the page reads it; `changes` holds each printed change as the result card lists
 * it, with its band and the second test's verdict over it (pageSays.ts).
 *
 * With OM_SURFACES=1 each row also checks that the card, the PDF and the figure say one count, one set
 * of coordinates and one pair of numbers (pageSays.surfaceProblems), one line per row on stderr.
 *
 * Run: node --experimental-strip-types audit/truth/run.ts <list.tsv> <shard> <shards> <out.jsonl>
 */
import { readFileSync, appendFileSync, existsSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'

const W = new URL('../../web/src/', import.meta.url).pathname
const ingest = await import(`${W}ingest.ts`)
const score = await import(`${W}scoreSample.ts`)
const keys = await import(`${W}answerKeys.ts`)
const inferredArray = await import(`${W}inferredArray.ts`)
const info = await import(`${W}informationRatio.ts`)
const stageMod = await import(`${W}stage.ts`)
const platform = await import(`${W}platform.ts`)
const parentage = await import(`${W}parentage.ts`)
const pageSays = await import(`${W}pageSays.ts`)
const sexing = await import(`${W}sexing.ts`)
const sexMaterial = await import(`${W}sexMaterial.ts`)
const overviewMod = await import(`${W}overview.ts`)
const table = await import(`${W}syngamyTable.ts`)
const reportPdf = await import(`${W}syngamyPdf.ts`)
const defectsMod = await import(`${W}defects.ts`)
const barsMod = await import(`${W}bars.ts`)
const withheldMod = await import(`${W}withheld.ts`)
const svgPdf = await import(`${W}svgPdf.ts`)
const plotSeries = await import(`${W}plotSeries.ts`)
const { gzipSync } = await import('node:zlib')
const buildInfo = readFileSync(new URL('../../originmarker/build_info.py', import.meta.url), 'utf8')
const TOOL = `OriginMarker ${/VERSION = "([^"]+)"/.exec(buildInfo)?.[1]} (${/CODENAME = "([^"]+)"/.exec(buildInfo)?.[1]})`
/** With OM_DUMP set, the result, the mother's result, the pair and the figure series of each record,
 *  for the figure gate, which renders the page's React figures where Vite is installed. */
const DUMP = process.env.OM_DUMP
const enc = (_k: string, v: unknown) => (v instanceof Map ? { __map: [...v] } : typeof v === 'number' && !Number.isFinite(v)
  ? { __num: String(v) } : v === undefined ? { __undef: 1 } : v)

const [listFile, shard, shards, out] = [process.argv[2], Number(process.argv[3]), Number(process.argv[4]), process.argv[5]]
const done = new Set<string>()
if (existsSync(out)) {
  for (const l of readFileSync(out, 'utf8').split('\n')) {
    if (l.startsWith('{')) { const r = JSON.parse(l); done.add(`${r.id}|${r.tag}`) }
  }
}
const lines = readFileSync(listFile, 'utf8').split('\n').filter(Boolean)
const lineOf = new WeakMap<string[], number>()
const list = lines.map((l, i) => { const f = l.split('\t'); lineOf.set(f, i + 1); return f }).filter((_, i) => i % shards === shard)

const text = (p: string) => { const raw = readFileSync(p); return (p.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8') }
function rows(path: string): any[] {
  const lines = text(path).split('\n')
  let h = 0; while (!lines[h] || lines[h].startsWith('#')) h++
  const map = ingest.headerMap(lines[h])
  if (!map) throw new ChipRefused(path, { text: platform.unreadableLines(lines, panel.freqB),
    causes: platform.missingRequired(lines[h]) })
  const o: any[] = []
  for (let i = h + 1; i < lines.length; i++) { const r = ingest.parseRow(lines[i], map); if (r) o.push(r) }
  checkChip(path, o, map)
  return o
}
const panel = keys.parsePanel(text(new URL(`../../web/public/${keys.FOUNDER_PANEL_PATH}`, import.meta.url).pathname))
class ChipRefused extends Error {
  file: string; call: any
  constructor(file: string, call: any) { super(call.text); this.file = file; this.call = call }
}
/** The same recognition the page runs, from the rows already read. Throws on another chip. */
function checkChip(path: string, rs: any[], map: any) {
  const b = ingest.emptyBuildSums()
  for (const r of rs) ingest.accumulateBuild(r, b)
  const call = platform.recogniseRows(rs, ingest.buildVerdict(b), map, panel.freqB)
  if (!call.proven) throw new ChipRefused(path, call)
}
const parents = new Map<string, any>()
function parentIndex(path: string) {
  if (path === '-') return score.noParent()
  if (parents.has(path)) return parents.get(path)
  const banner = inferredArray.parseInferredBanner(text(path).slice(0, inferredArray.MARK_BYTES)) ?? undefined
  const acc = score.emptyParent()
  const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); let first = ''
  for (const r of rows(path)) {
    if (!first) first = r.probesetId
    ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectParentRow(r, acc)
  }
  const idx = score.finishParent(acc, ingest.finishProfile(path, byChrom, bafSums, first).build.build, banner)
  parents.set(path, idx)
  return idx
}

const seg = (s: any) => ({ chrom: s.chrom, kind: s.kind, start: s.startBp, end: s.endBp, absentParent: s.absentParent,
  shortBy: s.shortBy,
  score: Number.isFinite(s.score) ? +s.score.toFixed(1) : undefined, lrrShift: Number.isFinite(s.lrrShift) ? +s.lrrShift.toFixed(3) : undefined,
  lrrZ: Number.isFinite(s.lrrZ) ? +s.lrrZ.toFixed(2) : undefined })
const held = (w: any) => (w && typeof w === 'object' && 'kind' in w && 'startBp' in w ? seg(w)
  : w && typeof w === 'object' ? (w.verdict ?? w.cls ?? w.state ?? 'object') : w)

// The called autosomal genotypes of another array of the same run, as the page keeps them: read by
// probe id, held as one byte per probe so a whole run's arrays fit in memory together.
const probeIndex = new Map<string, number>()
const CODES = ['', 'AA', 'AB', 'BB']
const callsOf = new Map<string, { get(id: string): string | undefined }>()
function calls(path: string): { get(id: string): string | undefined } {
  let m = callsOf.get(path)
  if (!m) {
    const rs = rows(path).filter((r) => /^\d+$/.test(r.chrom))
    for (const r of rs) if (!probeIndex.has(r.probesetId)) probeIndex.set(r.probesetId, probeIndex.size)
    const code = new Uint8Array(probeIndex.size)
    for (const r of rs) if (r.genotype !== 'NC') code[probeIndex.get(r.probesetId)!] = CODES.indexOf(r.genotype)
    m = { get: (id: string) => { const i = probeIndex.get(id); return i === undefined || i >= code.length ? undefined : CODES[code[i]] || undefined } }
    callsOf.set(path, m)
  }
  return m
}

for (const row of list) {
  const [id, path, parentPath, tag, otherPath, role, others] = row
  if (done.has(`${id}|${tag}`)) continue
  try {
    const pat = parentIndex(parentPath)
    const mat = otherPath ? parentIndex(otherPath) : null
    const acc = score.emptyCollected(pat, mat, panel as never)
    const byChrom = new Map(); const bafSums = ingest.emptyBafSums(); const builds = ingest.emptyBuildSums(); let first = ''
    const sameHom = keys.emptySameHom(); const homTally = info.emptyHomTally()
    for (const r of rows(path)) {
      if (!first) first = r.probesetId
      ingest.accumulate(r, byChrom); ingest.accumulateBaf(r, bafSums); score.collectRow(r, pat, mat, acc as never)
      ingest.accumulateBuild(r as never, builds)
      if (!/^\d+$/.test(r.chrom)) continue
      const p = panel.get(r.probesetId)
      if (p) { keys.addSameHom(sameHom, r.genotype, p, p); info.addHomozygous(homTally, r.genotype, p) }
    }
    const profile = ingest.finishProfile(id, byChrom, bafSums, first, builds)
    // The stage with no parent consulted: the founder panel's drop-in and information ratio only.
    const free = stageMod.inferStage({ ...profile, dropIn: keys.onParentScale(keys.dropIn(sameHom), panel.calibration),
      information: info.informationRatio(homTally) }) as any
    const res = await score.scoreSample({ acc, profile, pat, mat, soloRole: role === 'maternal' ? 'maternal' : 'paternal',
      sibs: (others ?? '').split(';').filter(Boolean).map(calls),
      sampleName: id, log: () => {} } as never) as any
    // The mother's genome-wide verdict, read as the page reads it, so the headline can be replayed.
    const y = sexing.sexCall((acc as any).sex)
    const maternal = (acc as any).tm && mat ? parentage.classify((acc as any).tm, mat.heterozygosity, { role: 'maternal', yBearing: y.yBearing }) : null
    if (process.env.OM_SURFACES && !res.noParent) {
      const cls = maternal ? parentage.pair(res, maternal, NaN).originClass : res.originClass
      const bad = pageSays.surfaceProblems(res, maternal ?? undefined, cls, tag.startsWith('recon'))
      process.stderr.write(`SURFACES ${tag} ${bad.length ? bad.join(' | ') : 'ok'}\n`)
    } else if (process.env.OM_SURFACES) {
      const bad = pageSays.surfaceProblems(res, undefined, res.originClass, false)
      process.stderr.write(`SURFACES ${tag} ${bad.length ? bad.join(' | ') : 'ok'}\n`)
    }
    // The result card's change list, as the overview counts it (web/src/overview.ts). Read last, after
    // every other field is taken.
    const paired = maternal && !res.noParent ? parentage.pair(res, maternal, NaN) : undefined
    const says = () => pageSays.cardSays(res, maternal ?? undefined, paired ? paired.originClass : res.originClass, tag.startsWith('recon'))
    const card = () => {
      const l = says().list
      return { counted: l.counted.map(cardRow), once: l.once.length, disturbed: !!l.disturbed,
        ...pageSays.cardCounts(l), countText: pageSays.countText(l), line: says().head.line,
        figure: { drawn: pageSays.drawnOf(res, says().defects), note: pageSays.figureNote(l, pageSays.drawnOf(res, says().defects)) ?? null },
        ...(l.notReported ? { notReported: l.notReported.rows.map(cardRow) } : {}) }
    }
    const surfaces = async () => {
      const s = says()
      const ov = new Map(overviewMod.overviewTable([overviewMod.overviewSampleOf(s.list, res)], new Map()))
      const tsv = table.syngamyTable([table.tableSampleOf(id, res, s)]).trim().split('\n').slice(1).map((l: string) => l.split('\t'))
      const at = (c: string) => tsv.map((x: string[]) => x[table.COLUMNS.indexOf(c)])
      const out: any = {
        overview: { changes: ov.get('changes with a gain or a loss'), findings: ov.get('findings (copy-neutral, uneven shares, complex genome)') },
        download: { rows: tsv.length, notReported: at('word').filter((w: string) => w === 'not reported').length,
          loadedParent: [...new Set(at('loaded_parent'))],
          ...(s.list.notReported ? { classes: [...new Set(at('class'))], parents: [...new Set(at('parent'))],
            evidence: [...new Set(at('evidence'))] } : {}) },
        // Every sentence the card prints of a row with neither a gain nor a loss.
        findingTexts: s.list.counted.filter((g: any) => pageSays.rowDirection(g) === 'none')
          .map((g: any) => [pageSays.whoseText(g), pageSays.rowSays(g, s.setup), pageSays.callTitle(g)].join(' | ')),
      }
      // The chromosome table's note where a whole call was withdrawn after it was written: it must
      // state the final call, as the change list does.
      if (!s.list.notReported) {
        const gone = (res.chroms ?? []).filter((c: any) => /^(LOST|GAINED): /.test(c.note ?? '') && !c.aneuploidy)
        if (gone.length) out.wholeWithdrawn = gone.map((c: any) => `chr${c.chrom}: ${pageSays.chromCells(false, c, s.defects).note}`)
      }
      if (s.list.notReported) {
        const nr = s.list.notReported
        // The run log's gain, loss and stretch lines of this sample, as the page prints them.
        out.log = pageSays.changeLogLines(res, id, res.role).map((x: any) => x[1])
        out.seenWhose = nr.rows.map((g: any) => pageSays.whoseText(g))
        // Every state cell of the page's and the report's tables of this sample, as each prints it.
        const P = pageSays
        const sexRows = (res.chroms ?? []).filter((c: any) => c.chrom === 'X' || c.chrom === 'Y' || c.chrom === '23')
        const chromRows = [...sexRows, ...(res.chroms ?? []), ...(maternal?.chroms ?? [])]
        const chromCells = chromRows.flatMap((c: any) => Object.values(P.chromCells(true, c)) as string[])
        const segCells = (res.segments ?? []).map((sg: any) => P.segmentKindText(true, sg, res.role))
        const gainCells = (res.gains ?? []).flatMap((g: any) => Object.values(P.gainCells(true, g)) as string[])
        const aneu = (res.chroms ?? []).filter((c: any) => c.aneuploidy || c.lostParent || c.mixed
          || c.alleleGain?.verdict === 'gain' || c.alleleGain?.verdict === 'imbalance')
        out.tables = {
          page: { chrom: chromCells, segment: segCells, gain: gainCells },
          pdf: { chrom: chromCells, aneuploidy: aneu.flatMap((c: any) => P.aneuploidyCells(true, c)), segment: segCells, gain: gainCells,
            finding: (res.findings ?? []).map((f: any) => P.stateWord(true, '', P.kindSeen(f.cls))),
            origin: [...(res.oneParent ?? []), ...(res.dosageCalls ?? [])].map((o: any) => P.verdictText(true, o.verdict)) },
          rows: { chrom: chromRows.length, maternalChrom: maternal?.chroms?.length ?? 0, aneuploidy: aneu.length,
            segment: (res.segments ?? []).length, gain: (res.gains ?? []).length, finding: (res.findings ?? []).length,
            origin: (res.oneParent ?? []).length + (res.dosageCalls ?? []).length },
        }
        // Every sentence a card not reportable prints about the sample and its rows.
        out.notReportedTexts = [s.head.line, s.title, pageSays.figureNote(s.list, pageSays.drawnOf(res, s.defects)) ?? '',
          ...defectsMod.runAlerts(res).flatMap((a: any) => [a.headline, a.body]),
          pageSays.disturbedText(nr, false),
          ...nr.rows.flatMap((g: any) => [pageSays.whoseText(g), pageSays.rowSays(g, s.setup), pageSays.callTitle(g),
            ...g.pieces.map((p: any) => pageSays.pieceText(p, g)), ...pageSays.pdfCells(g, s.setup),
            pageSays.seenMarkText({ ...g, dir: pageSays.rowDirection(g) }).sentence])].filter(Boolean)
      }
      {
        const files: any[] = [{ name: 'sample', size: 1, sha256: 'x', role: 'sample', markers: 1, result: res,
          maternal: maternal ?? undefined, paired, inferred: false }]
        if (!res.noParent && res.role === 'paternal') files.unshift({ name: 'donor', size: 1, sha256: '0'.repeat(64), role: 'donor', markers: 1 })
        if (!res.noParent && (maternal || res.role === 'maternal')) files.unshift({ name: 'oocyte', size: 1, sha256: '0'.repeat(64), role: 'oocyte', markers: 1 })
        const txt = Buffer.from(await (await reportPdf.buildReportPdf({
          files,
          donorHeterozygosity: NaN, startedAt: null, generatedAt: '', tool: TOOL, reportId: 'GATE', fromExamples: false,
          runArrays: (others ?? '').split(';').filter(Boolean).length,
        } as never)).arrayBuffer()).toString('latin1')
        // Every string the PDF draws, unescaped; a wrapped sentence is its lines joined by spaces.
        const strings = [...txt.matchAll(/\(((?:\\.|[^\\)])*)\) Tj/g)].map((m) => m[1]
          .replace(/\\([0-7]{3})|\\(.)/g, (_, o, c) => (o ? String.fromCharCode(parseInt(o, 8)) : c)))
        const flat = (x: string) => x.replace(/[^\x20-\x7e]/g, '').replace(/\s+/g, ' ').trim()
        const all = flat(strings.join(' '))
        out.pdf = { notReported: strings.filter((x) => x === 'not reported').length,
          bar: strings.includes(`${s.head.word}${pageSays.countText(s.list)}`),
          sexLine: res.sexChromosome ? all.includes(flat(pageSays.sexLine(res.sexChromosome))) : null,
          inferred: all.includes('inferred from paternal absence'), noParentLine: all.includes('sample no parent loaded'),
          vsDonor: all.includes('sample vs donor'),
          noParentMethods: all.includes('No parental array was supplied.'),
          inferMethods: all.includes('No oocyte donor array was supplied. A maternal origin below is therefore inferred'),
          pairedMethods: all.includes('An oocyte donor array was supplied'),
          donorText: (all.match(/.{0,40}(sperm donor's|donor's obligate|inferred from).{0,40}/g) ?? []),
          findingsBanner: all.includes('Findings, neither a gain nor a loss, in 1 of 1 sample(s)'),
          // Strings drawn whole that are a state word of a reportable sample's tables.
          stateStrings: s.list.notReported ? strings.filter((x) => /^(LOSS|GAIN|GAIN \(ALLELES\)|PARENT ABSENT|IMBALANCE|MIXED|DNA absent|extra copies|(paternal|maternal) alleles absent|alleles absent, both parents)$/.test(x)
            || /^(LOST|GAINED):/.test(x)) : [],
          seenHeadings: s.list.notReported ? strings.filter((x) => x.endsWith(': readings seen on this array, not reported')).length : 0,
          changeBanner: (all.match(/Chromosomal change in 1 of 1 sample\(s\): (.*?)\. Coordinates/) ?? [])[1] ?? null }
        // What the outputs print: the withheld count, the bars listed, the methods paragraph, the file itself.
        const input = { files, tool: TOOL, runArrays: (others ?? '').split(';').filter(Boolean).length }
        const facts = svgPdf.pdfFacts(txt)
        const names = new Set(barsMod.BARS.map((b: any) => b.name))
        const methods = reportPdf.reportMethods(input as never)
        out.outputs = {
          withheldHeading: strings.find((x) => x.startsWith('Withheld, and why: ')) ?? null,
          withheldRows: withheldMod.withheldRows(res, s.list).length,
          barsListed: [...new Set(strings.filter((x) => names.has(x)))],
          barsApplied: [...new Set(barsMod.barsApplied(reportPdf.reportSetup(input as never)).map((b: any) => b.name))],
          barValues: barsMod.barsApplied(reportPdf.reportSetup(input as never)).map((b: any) => barsMod.barValue(b)),
          typedTable: all.includes('Constants and where they come from'),
          methods: methods.text, cites: methods.cites.map((c: any) => c.doi),
          methodsPrinted: all.includes(flat(methods.text)),
          xrefOk: facts.xrefOk, images: facts.images, pages: facts.pages,
          footers: strings.filter((x) => x.includes(` | ${TOOL} | report `) && x.includes(' | page ')).length,
        }
        if (DUMP) {
          const series = { ...plotSeries.plotData(acc as never, { loaded: pat.absent ? null : { role: role === 'maternal' ? 'maternal' : 'paternal', gt: pat.gt },
            second: mat?.gt ?? null }), seen: pageSays.notReportable(res) }
          const line = lineOf.get(row) ?? 0
          appendFileSync(`${DUMP}/dump-${tag}-${line}.json.gz`, gzipSync(JSON.stringify({ line, tag, others: (others ?? '').split(';').filter(Boolean).length,
            res, maternal: maternal ?? null, paired: paired ?? null, series, tool: TOOL }, enc)))
        }
      }
      return out
    }
    const cardRow = (g: any) => ({ chrom: g.chrom, gain: g.gain, whole: g.whole, origin: g.origin, start: g.startBp, end: g.endBp,
      word: g.word, pieces: g.pieces.map((p: any) => ({ kind: p.kind, start: p.startBp, end: p.endBp })) })
    appendFileSync(out, JSON.stringify({
      id, tag, parent: parentPath.split('/').pop(),
      callRate: profile.callRate, hetRate: profile.hetRate,
      stage: { stage: res.stage?.stage, refusedBy: res.stage?.refusedBy ?? '', ploidy: res.stage?.ploidy ?? '',
        dropIn: res.stage?.dropInCorrected ?? null },
      stageFree: { stage: free.stage, refusedBy: free.refusedBy ?? '', ploidy: free.ploidy ?? '',
        dropIn: free.dropInCorrected ?? null },
      information: info.informationRatio(homTally),
      parentGenome: res.parentGenome ?? null, parentGenomeOther: res.parentGenomeOther ?? null,
      zygosity: res.zygosity, originClass: res.originClass,
      sex: res.sexChromosome?.constitution ?? '',
      // The chrY reading and what the page prints of sex: the line's reason, the sample's sex word,
      // the sperm type and the parent checks.
      y: { yBearing: y.yBearing, callRatio: y.callRatio, lrrShift: y.lrrShift, yBafShare: y.yBafShare ?? null,
        bafFloor: y.bafFloor ?? null, why: y.why ?? null },
      sexWhy: res.sexChromosome?.why ?? '', sexWord: sexMaterial.shownSex(profile.sex, res.sexMaterial, res.sexChromosome),
      spermType: res.spermType ?? null, parentSanity: res.parentSanity ?? null, matSanity: res.parentSanityOther ?? null,
      parentSex: { pat: pat.sex ?? null, mat: mat?.sex ?? null },
      chroms: (res.chroms ?? []).map((c: any) => ({ chrom: c.chrom, aneuploidy: c.aneuploidy ?? '',
        lostParent: c.lostParent ?? '', aneuploidyParent: c.aneuploidyParent ?? '', gain: c.alleleGain?.verdict ?? '', gainParent: c.alleleGain?.parent ?? '',
        lrrShift: c.lrrShift ?? null, lrrZ: c.lrrZ ?? null, state: c.copyState?.state ?? '' })),
      segments: (res.segments ?? []).map(seg),
      findings: (res.findings ?? []).map((f: any) => ({ cls: f.cls, chrom: f.chrom, start: f.startBp, end: f.endBp,
        whole: !!f.wholeChromosome, originBlocked: f.originBlocked ? true : undefined, flag: f.flag })),
      origins: [...(res.oneParent ?? []).map((o: any) => ({ by: 'one-parent', where: o.where, verdict: o.verdict, parent: o.parent ?? null })),
        ...(res.dosageCalls ?? []).map((d: any) => ({ by: 'dosage', where: d.where, verdict: d.verdict }))],
      refusals: (res.refusals ?? []).map((x: any) => ({ claim: x.claim, chrom: x.chrom ?? '', gate: x.gate,
        cause: x.cause, withheld: held(x.withheld) })),
      genomeVerdicts: { loaded: res.verdict ?? null, mother: maternal?.verdict ?? null },
      // Each printed change with the band it carries and the second test's verdict over it.
      changes: pageSays.defectsForResult(res).map((d: any) => ({ chrom: d.chrom, start: d.startBp, end: d.endBp,
        kind: d.kind, origin: d.origin, step: d.step ?? null, band: d.band ?? null, second: d.second?.verdict ?? null })),
      ...(res.motherFromRun ? { motherFromRun: { ...res.motherFromRun,
        drawArrays: res.motherFromRun.drawArrays.map((i: number) => (others ?? '').split(';').filter(Boolean)[i]) } } : {}),
      card: card(),
      surfaces: await surfaces(),
    }) + '\n')
  } catch (e) {
    if (e instanceof ChipRefused) {
      appendFileSync(out, JSON.stringify({ id, tag, chip: { proven: false, file: e.file.split('/').pop(),
        knownShare: e.call.knownShare ?? null, causes: e.call.causes, text: e.call.text } }) + '\n')
    } else appendFileSync(out, JSON.stringify({ id, tag, error: String(e).slice(0, 300) }) + '\n')
  }
}
appendFileSync(out, `DONE ${shard}\n`)
