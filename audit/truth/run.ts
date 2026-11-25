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
 * `copyReference` is the result's own. Where the row names other arrays of its run, `sisterGains` is
 * the page's sister reading (sisterGain.ts) over this array and those others, each other scored as
 * its own row would be: the same parents, and the rest of the run as its sibs. `from` is the index in
 * the `others` column of the array that shows the loss; `sisterUnread` lists others that could not be
 * scored, which the page leaves out of the run the same way. A run's arrays are held as their window
 * medians, losses and one byte per probe, never as whole results.
 *
 * With OM_SURFACES=1 each row also checks that the card, the PDF and the figure say one count, one set
 * of coordinates and one pair of numbers (pageSays.surfaceProblems), one line per row on stderr.
 *
 * Run: node --experimental-strip-types audit/truth/run.ts <list.tsv> <shard> <shards> <out.jsonl>
 */
import { readFileSync, appendFileSync, existsSync } from 'node:fs'
import { ChipRefused, panel, calls, chromOf, collectArray, scoreArray as scoreFresh } from './score-array.ts'
import { openStore, openList, treeKey, only } from '../gate/store.ts'

const W = new URL('../../web/src/', import.meta.url).pathname
const keys = await import(`${W}answerKeys.ts`)
const info = await import(`${W}informationRatio.ts`)
const stageMod = await import(`${W}stage.ts`)
const parentage = await import(`${W}parentage.ts`)
const pageSays = await import(`${W}pageSays.ts`)
const sexing = await import(`${W}sexing.ts`)
const sexMaterial = await import(`${W}sexMaterial.ts`)
const overviewMod = await import(`${W}overview.ts`)
const table = await import(`${W}syngamyTable.ts`)
const reportPdf = await import(`${W}syngamyPdf.ts`)
const defectsMod = await import(`${W}defects.ts`)
const S = await import(`${W}sisterGain.ts`)

const [listFile, shard, shards, out] = [process.argv[2], Number(process.argv[3]), Number(process.argv[4]), process.argv[5]]
const done = new Set<string>()
if (existsSync(out)) {
  for (const l of readFileSync(out, 'utf8').split('\n')) {
    if (l.startsWith('{')) { const r = JSON.parse(l); done.add(`${r.id}|${r.tag}`) }
  }
}
const list = readFileSync(listFile, 'utf8').split('\n').filter(Boolean).map((l) => l.split('\t'))
  .filter((_, i) => i % shards === shard)

const seg = (s: any) => ({ chrom: s.chrom, kind: s.kind, start: s.startBp, end: s.endBp, absentParent: s.absentParent,
  shortBy: s.shortBy,
  score: Number.isFinite(s.score) ? +s.score.toFixed(1) : undefined, lrrShift: Number.isFinite(s.lrrShift) ? +s.lrrShift.toFixed(3) : undefined,
  lrrZ: Number.isFinite(s.lrrZ) ? +s.lrrZ.toFixed(2) : undefined })
const held = (w: any) => (w && typeof w === 'object' && 'kind' in w && 'startBp' in w ? seg(w)
  : w && typeof w === 'object' ? (w.verdict ?? w.cls ?? w.state ?? 'object') : w)

// With OM_STORE, a scored array is read from the store (audit/gate/store.ts) where one is kept for
// the same scoring closure and inputs, and kept there where not. An array whose inputs are not all
// on the open list (OM_OPEN) is scored and never stored. What the record reads of the
// parents and the collected tallies is held by name. The per-marker LRR tracks the sister reading
// reads are kept for an array scored with a run, and collected again from the array, not scored,
// for one scored without.
const store = process.env.OM_STORE
  ? openStore(process.env.OM_STORE, treeKey(new URL('../..', import.meta.url).pathname, 'rows'), openList(process.env.OM_OPEN))
  : null
async function scoreArray(id: string, path: string, parentPath: string, otherPath: string | undefined,
  role: string | undefined, sibPaths: string[]) {
  const inputs = [path, parentPath, otherPath ?? '', ...sibPaths]
  if (!store || !store.allows(inputs)) return scoreFresh(id, path, parentPath, otherPath, role, sibPaths)
  const key = store.key(['run.ts', id, role ?? null, otherPath === undefined], inputs)
  let v = store.get(key)
  if (!v) {
    const s = await scoreFresh(id, path, parentPath, otherPath, role, sibPaths)
    v = { res: s.res, profile: s.profile, sameHom: s.sameHom, homTally: s.homTally, maternal: s.maternal,
      acc: { sex: (s.acc as any).sex, ...(sibPaths.length ? { cnByChrom: (s.acc as any).cnByChrom } : {}) },
      pat: { sex: s.pat.sex }, mat: s.mat ? { sex: s.mat.sex, heterozygosity: s.mat.heterozygosity } : null }
    store.put(key, v)
    return s
  }
  return { res: v.res, profile: v.profile, sameHom: v.sameHom, homTally: v.homTally, maternal: v.maternal, pat: only(v.pat), mat: v.mat && only(v.mat),
    acc: only(v.acc, { cnByChrom: () => collectArray(id, path, parentPath, otherPath).acc.cnByChrom }) }
}

// Each array's part in the sister reading, keyed by what it was scored with, so an array read once
// in a shard is not scored again for every row of its run. The key holds the sibs as a set, which
// assumes the order of a run's arrays changes no loss call.
//
// An other array that has its own row later in this shard, with the same parents and the same run,
// is scored as that row and its record kept, so it is scored once and its own record is the one
// the reading used.
const unitOf = new Map<string, any>()
const unitKey = (path: string, parentPath: string, otherPath: string | undefined, role: string | undefined,
  sibPaths: string[]) => [path, parentPath, otherPath ?? '', role ?? '', [...sibPaths].sort().join(';')].join('|')
const sisterUnit = (res: any, acc: any, path: string) => S.sisterUnitOf(res, S.lrrWindows(acc.cnByChrom), calls(path), chromOf)
const kept = new Map<string, { rec: any; unit: any }>()
async function otherUnit(path: string, parentPath: string, otherPath: string | undefined, role: string | undefined,
  sibPaths: string[]) {
  const key = unitKey(path, parentPath, otherPath, role, sibPaths)
  if (!unitOf.has(key)) {
    let u = null
    try {
      const own = rowOf.get(key)
      if (own) {
        const [id, , , tag, , , others] = own
        const got = await recordOf(id, path, parentPath, tag, otherPath, role, (others ?? '').split(';').filter(Boolean))
        kept.set(`${id}|${tag}`, got)
        u = got.unit
      } else {
        const { acc, res } = await scoreArray(path.split('/').pop() as string, path, parentPath, otherPath, role, sibPaths)
        u = sisterUnit(res, acc, path)
      }
    } catch { u = null }
    unitOf.set(key, u)
  }
  return unitOf.get(key)
}

/** The record of one row, and this array's own part in the sister reading where its run has others. */
async function recordOf(id: string, path: string, parentPath: string, tag: string, otherPath: string | undefined,
  role: string | undefined, run: string[]) {
  const { pat, mat, acc, profile, sameHom, homTally, res, maternal } = await scoreArray(id, path, parentPath, otherPath, role, run)
  // The stage with no parent consulted: the founder panel's drop-in and information ratio only.
  const free = stageMod.inferStage({ ...profile, dropIn: keys.onParentScale(keys.dropIn(sameHom), panel.calibration),
    information: info.informationRatio(homTally) }) as any
  // The chrY reading. The mother's genome-wide verdict comes with the scored array (score-array.ts).
  const y = sexing.sexCall((acc as any).sex)
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
        donorHeterozygosity: NaN, startedAt: null, generatedAt: '', tool: 'OriginMarker', reportId: 'GATE', fromExamples: false,
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
    }
    return out
  }
  const cardRow = (g: any) => ({ chrom: g.chrom, gain: g.gain, whole: g.whole, origin: g.origin, start: g.startBp, end: g.endBp,
    word: g.word, pieces: g.pieces.map((p: any) => ({ kind: p.kind, start: p.startBp, end: p.endBp })) })
  const rec: any = {
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
      drawArrays: res.motherFromRun.drawArrays.map((i: number) => run[i]) } } : {}),
    copyReference: res.copyReference ?? null,
    card: card(),
    surfaces: await surfaces(),
  }
  if (!run.length) return { rec, unit: null }
  const unit = sisterUnit(res, acc, path)
  unitOf.set(unitKey(path, parentPath, otherPath, role, run), unit)
  return { rec, unit }
}

// This shard's rows still to write that have a run, by the key their array would be scored under.
const rowOf = new Map<string, string[]>()
for (const row of list) {
  const [id, path, parentPath, tag, otherPath, role, others] = row
  const run = (others ?? '').split(';').filter(Boolean)
  if (run.length && !done.has(`${id}|${tag}`)) rowOf.set(unitKey(path, parentPath, otherPath, role, run), row)
}

for (const [id, path, parentPath, tag, otherPath, role, others] of list) {
  if (done.has(`${id}|${tag}`)) continue
  try {
    const run = (others ?? '').split(';').filter(Boolean)
    const got = kept.get(`${id}|${tag}`)
    kept.delete(`${id}|${tag}`)
    const { rec, unit } = got ?? await recordOf(id, path, parentPath, tag, otherPath, role, run)
    if (unit) {
      // Each other array of the run, scored as its own row would be. The units are read in this
      // array's order, so a reading's `from` maps back to its place in `others`.
      const units = [unit]; const from: number[] = []; const unread: number[] = []
      for (const [k, p] of run.entries()) {
        const u = await otherUnit(p, parentPath, otherPath, role, [path, ...run].filter((q) => q !== p))
        if (u) { units.push(u); from.push(k) } else unread.push(k)
      }
      rec.sisterGains = S.sisterGainsFor(units, 0).map((g: any) => ({ chrom: g.chrom, start: g.startBp,
        end: g.whole ? null : g.endBp, whole: g.whole, shift: +g.shift.toFixed(4), z: +g.z.toFixed(3),
        windows: g.windows, from: from[g.from - 1] }))
      if (unread.length) rec.sisterUnread = unread
    }
    appendFileSync(out, JSON.stringify(rec) + '\n')
  } catch (e) {
    if (e instanceof ChipRefused) {
      appendFileSync(out, JSON.stringify({ id, tag, chip: { proven: false, file: e.file.split('/').pop(),
        knownShare: e.call.knownShare ?? null, causes: e.call.causes, text: e.call.text } }) + '\n')
    } else appendFileSync(out, JSON.stringify({ id, tag, error: String(e).slice(0, 300) }) + '\n')
  }
}
if (store) process.stderr.write(`STORE ${JSON.stringify(store.stats)}\n`)
appendFileSync(out, `DONE ${shard}\n`)
