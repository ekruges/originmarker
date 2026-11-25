/**
 * The fragile site a printed break falls in, by Table 3 of Georgakilas et al. 2014 (hg19).
 *
 * The same rule as the break test in audit/truth/fragile_sites.py, so the page and the test agree:
 * a break is an edge of a piece of a change the card counts (a gain or a loss of part of a
 * chromosome; not a row seen once, a finding, or anything on a sample not reportable), an edge that is not within 0.5 Mb of the chip's first or last marker on
 * its chromosome, and not the unplaced edge of an arm-level call; it is in a site when
 * start <= coordinate <= end. The sites must equal audit/truth/fragile-sites-georgakilas2014.tsv
 * (fragileSites.check.ts compares them). Nothing here is read by the scoring.
 */
import { CITATIONS } from './citations.ts'

export const FRAGILE_CITE = CITATIONS.georgakilas_2014

/** Georgakilas et al. 2014, Table 3: name, chromosome, start, end (GRCh37/hg19). */
export const FRAGILE_SITES: readonly (readonly [string, string, number, number])[] = [
  ['FRA1E', '1', 97749961, 98119925], ['FRA1H', '1', 216110179, 226742379],
  ['FRA2Ccen', '2', 18524760, 19272430], ['FRA2Ctel', '2', 14934889, 15681539],
  ['FRA2G', '2', 169498510, 170313244], ['FRA2H', '2', 186716111, 187251132],
  ['FRA3B', '3', 59623632, 63846635], ['FRA4F', '4', 90208191, 97312924],
  ['FRA6E', '6', 160275245, 166079958], ['FRA6F', '6', 111577970, 112568919],
  ['FRA6H', '6', 27804038, 37156663], ['FRA7B', '7', 1, 12392296],
  ['FRA7E', '7', 80508751, 84935939], ['FRA7G', '7', 115894865, 116072877],
  ['FRA7H', '7', 130413791, 130857950], ['FRA7I', '7', 144671106, 146121417],
  ['FRA7K', '7', 110657856, 111031412], ['FRA8C', '8', 124285269, 128421245],
  ['FRA9G', '9', 17135038, 17503917], ['FRA10F', '10', 125391503, 128256865],
  ['FRA11E', '11', 32086458, 34028916], ['FRA11G', '11', 113688799, 118157704],
  ['FRA13A', '13', 35546088, 36184897], ['FRA13E', '13', 73285774, 76386498],
  ['FRA16D', '16', 78420359, 78731553], ['FRAXB', 'X', 6595111, 7548235],
]

/** The chip's first and last marker on each chromosome (hg19), from its probe annotation. */
const CHIP_ENDS: Record<string, readonly [number, number]> = {
  1: [86028, 249222527], 2: [11944, 243041411], 3: [72365, 197848556], 4: [59350, 190906015],
  5: [14782, 180698413], 6: [203397, 170907734], 7: [41892, 159124173], 8: [164984, 146296404],
  9: [16967, 141101939], 10: [93190, 135440226], 11: [193146, 134945120], 12: [190980, 133831319],
  13: [19020095, 115103150], 14: [19071913, 107287663], 15: [20044342, 102397317],
  16: [85629, 90170095], 17: [6157, 81103682], 18: [11358, 78015180], 19: [110783, 59097933],
  20: [63231, 62960292], 21: [10873592, 48099610], 22: [16057417, 51224208], X: [60425, 155233115],
}

/** An edge this close to the chip's first or last marker is the chromosome's end, not a break. */
export const END_BP = 500_000

/** The site a printed edge falls in, or undefined when the edge is not a break or is in no site. */
export function fragileSiteAt(chrom: string, x: number): string | undefined {
  const ends = CHIP_ENDS[chrom]
  if (!ends || x - ends[0] <= END_BP || ends[1] - x <= END_BP) return undefined
  return FRAGILE_SITES.find(([, c, s, e]) => c === chrom && s <= x && x <= e)?.[0]
}

/**
 * The plain line for a piece whose edges fall in a site, e.g. "Break at 61.2 Mb in fragile site
 * FRA3B (Georgakilas et al. 2014, Table 3)", or undefined. An arm-level call has no placed edge.
 */
export function fragileSiteText(p: { chrom: string; startBp: number; endBp: number; arm?: string }): string | undefined {
  if (p.arm) return undefined
  const hits = [p.startBp, p.endBp]
    .map((x) => [x, fragileSiteAt(p.chrom, x)] as const)
    .filter(([, s]) => s)
    .map(([x, s]) => `${(x / 1e6).toFixed(1)} Mb in fragile site ${s}`)
  if (!hits.length) return undefined
  return `Break at ${hits.join(', and at ')} (${FRAGILE_CITE.authors} ${FRAGILE_CITE.year}, Table 3)`
}

const BREAK_KINDS = new Set(['copy-loss', 'copy-gain', 'parental-absence'])

type Piece = { chrom: string; startBp: number; endBp: number; kind: string; arm?: string }

/** The pieces whose edges are breaks: pieces of the card's counted rows that are not whole, gains and
 *  losses only, none on a sample not reportable. */
export function countedBreakPieces<P extends Piece>(list: {
  counted: readonly { whole: boolean; pieces: readonly P[] }[]; notReported?: unknown
}): P[] {
  if (list.notReported) return []
  return list.counted.filter((g) => !g.whole).flatMap((g) => g.pieces).filter((p) => BREAK_KINDS.has(p.kind))
}

/** One line per counted piece with a break in a site, e.g. "chr3: Break at 61.2 Mb in fragile site FRA3B (...)". */
export function fragileLines(list: Parameters<typeof countedBreakPieces>[0]): string[] {
  return countedBreakPieces(list).flatMap((p) => {
    const t = fragileSiteText(p)
    return t ? [`chr${p.chrom}: ${t}`] : []
  })
}
