/**
 * The carrier's own genotypes, read in the browser and never sent anywhere.
 *
 * This is Layer B steps 1 and 2 done for you: the panel ranks candidates on expected
 * heterozygosity, which is a POPULATION prior and says nothing about the person being
 * tested, and the protocol then tells you to genotype that carrier and keep only the
 * markers where they are actually heterozygous. Given their VCF or array export, that
 * filter is a lookup.
 *
 * IT NEVER LEAVES THE BROWSER. The terms promise no patient data is submitted or retained,
 * and a carrier's WGS is the most identifying file in this workflow, so the parsing is here
 * rather than behind an upload endpoint. There is no endpoint to send it to.
 *
 * What it cannot do is phase. Knowing the carrier is heterozygous at a marker does not say
 * which of their two chromosomes each allele sits on, and short-read WGS cannot bridge the
 * tens of kilobases to the pathogenic variant either. That still needs a relative, or reads
 * long enough to span the distance. See the docs section on choosing a platform.
 *
 * Self-check: node src/genotypes.check.ts
 */

/** What the file says about this carrier at one site.
 *
 *  `absent` is its own answer and must never be read as hom-ref: a variants-only VCF omits
 *  every site where the sample matches the reference, so absence there usually means
 *  hom-ref, and in a no-call region it means nothing was measured. Those are different
 *  facts and only an all-sites/gVCF file distinguishes them. */
export type Call = 'het' | 'hom' | 'nocall' | 'absent'

export interface SiteCall {
  call: Call
  /** The genotype as the file wrote it: "0/1", "AG". Displayed, never re-derived. */
  text: string
  /** Where in the file this came from, for the build cross-check. */
  pos: number
  chrom: string
}

export type FileFormat = 'vcf' | 'array'

export interface GenotypeSet {
  format: FileFormat
  /** Keyed by rsID (lowercased) and by "chrom:pos". rsID is preferred at lookup because it
   *  survives a coordinate build difference; position is the fallback for sites the file
   *  left unnamed. */
  byRsid: Map<string, SiteCall>
  byPos: Map<string, SiteCall>
  /** VCF sample columns found, and which one was read. */
  samples: string[]
  sample: string | null
  linesRead: number
  sitesKept: number
}

/** '11', 'chr11', 'CHR11' are one chromosome. */
export const normChrom = (c: string) => String(c).trim().toLowerCase().replace(/^chr/, '')

const posKey = (chrom: string, pos: number) => `${normChrom(chrom)}:${pos}`

/**
 * A GT field to a call.
 *
 * Phase separators are read but the phase is discarded: '0|1' and '0/1' are both simply
 * heterozygous here. A VCF's phase set is local to the caller's own blocks and says nothing
 * about which parental chromosome carries what, which is the only phase this app cares
 * about and the one thing it never claims to know.
 */
export function callFromGT(gt: string): Call {
  const alleles = gt.trim().split(/[/|]/)
  if (alleles.length < 2) return 'nocall'
  if (alleles.some((a) => a === '.' || a === '')) return 'nocall'
  return alleles.every((a) => a === alleles[0]) ? 'hom' : 'het'
}

/**
 * An array export's genotype string to a call: "AG" heterozygous, "AA" homozygous.
 *
 * Heterozygosity is strand-invariant, which is what makes array data usable here at all:
 * a file reported on the opposite strand turns AG into TC, and both are still one allele
 * of each. That is the whole question this module asks, so the A/T and C/G ambiguity that
 * plagues array merging cannot produce a wrong het/hom answer. It would matter for allele
 * IDENTITY, which is a phasing concern and is not decided here.
 */
export function callFromArray(g: string): Call {
  const s = g.trim().toUpperCase().replace(/[^ACGTDI-]/g, '')
  if (!s || s.length < 2 || s.includes('-')) return 'nocall'
  const a = s.slice(0, 1), b = s.slice(1, 2)
  return a === b ? 'hom' : 'het'
}

export interface ParseOpts {
  /** Only sites on this chromosome inside [lo, hi] are kept. A whole-genome file is mostly
   *  irrelevant to one panel, and holding it would cost hundreds of megabytes of memory. */
  chrom: string
  lo: number
  hi: number
  /** Which VCF sample column to read, by name. Defaults to the first. */
  sample?: string
}

/** One line of a VCF or an array export, folded into the set. Exported for the self-check. */
export function ingestLine(line: string, set: GenotypeSet, opts: ParseOpts): void {
  if (!line || line.startsWith('##')) return

  // VCF header: remember the sample columns so a multi-sample file can be pointed at one.
  if (line.startsWith('#CHROM')) {
    const cols = line.split(/\t/)
    set.format = 'vcf'
    set.samples = cols.slice(9).map((s) => s.trim()).filter(Boolean)
    set.sample = opts.sample && set.samples.includes(opts.sample)
      ? opts.sample
      : (set.samples[0] ?? null)
    return
  }
  if (line.startsWith('#')) return           // array exports comment their header with '#'

  const f = line.split(/[\t ]+/)
  if (f.length < 4) return

  let chrom: string, pos: number, rsid: string, call: Call, text: string

  if (set.format === 'vcf') {
    // #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <samples...>
    if (f.length < 10 || !set.sample) return
    chrom = f[0]; pos = Number(f[1]); rsid = f[2]
    const gtIndex = f[8].split(':').indexOf('GT')
    if (gtIndex < 0) return
    const col = 9 + Math.max(0, set.samples.indexOf(set.sample))
    const gt = (f[col] ?? '').split(':')[gtIndex] ?? ''
    text = gt
    call = callFromGT(gt)
  } else {
    // rsid  chromosome  position  genotype   (23andMe / AncestryDNA style array export)
    rsid = f[0]; chrom = f[1]; pos = Number(f[2])
    text = f.slice(3).join('')
    call = callFromArray(text)
  }

  if (!Number.isFinite(pos)) return
  set.linesRead++
  // The window is the only filter. A site outside it cannot be a marker in this panel.
  if (normChrom(chrom) !== normChrom(opts.chrom) || pos < opts.lo || pos > opts.hi) return

  const site: SiteCall = { call, text, pos, chrom }
  set.sitesKept++
  if (rsid && rsid !== '.' && /^rs\d+$/i.test(rsid)) set.byRsid.set(rsid.toLowerCase(), site)
  set.byPos.set(posKey(chrom, pos), site)
}

export const emptySet = (): GenotypeSet => ({
  format: 'array',
  byRsid: new Map(), byPos: new Map(),
  samples: [], sample: null, linesRead: 0, sitesKept: 0,
})

/** Parse whole text. The streaming path in the component shares ingestLine with this. */
export function parseGenotypes(text: string, opts: ParseOpts): GenotypeSet {
  const set = emptySet()
  // A VCF announces itself; anything else is read as an array export.
  if (/^##fileformat=VCF/im.test(text)) set.format = 'vcf'
  for (const line of text.split(/\r?\n/)) ingestLine(line, set, opts)
  return set
}

/** What the file says about one marker, and how it was matched. */
export interface MarkerCall extends SiteCall {
  matchedBy: 'rsid' | 'position'
}

/**
 * Look a marker up. rsID first: it is stable across coordinate builds, so it still finds
 * the site in a GRCh37 file, where matching on position would land on the wrong one or on
 * nothing. Position is the fallback for a file that leaves the ID column empty.
 */
export function lookup(
  m: { rsid?: string | null; chrom: string; pos: number }, set: GenotypeSet,
): MarkerCall | null {
  const byId = m.rsid ? set.byRsid.get(m.rsid.toLowerCase()) : undefined
  if (byId) return { ...byId, matchedBy: 'rsid' }
  const byPos = set.byPos.get(posKey(m.chrom, m.pos))
  return byPos ? { ...byPos, matchedBy: 'position' } : null
}

export interface BuildCheck {
  checked: number
  agreed: number
  /** True when rsID-matched sites systematically sit at different coordinates, which means
   *  the file is on another assembly. */
  mismatch: boolean
  offsets: number[]
}

/**
 * Is this file on the same assembly as the panel?
 *
 * The panel is GRCh38 (R6). A GRCh37 export intersects by rsID perfectly well and by
 * POSITION lands somewhere else entirely, so a file matched mostly by position would return
 * confident genotypes for the wrong sites. rsID-matched sites give the answer for free:
 * their coordinates should agree, and if they consistently do not, the file is another
 * build. Reported rather than corrected, since lifting over here would be a second source
 * of coordinates and R1 keeps that to the live APIs.
 */
export function checkBuild(
  markers: { rsid?: string | null; chrom: string; pos: number }[], set: GenotypeSet,
): BuildCheck {
  let checked = 0, agreed = 0
  const offsets: number[] = []
  for (const m of markers) {
    if (!m.rsid) continue
    const hit = set.byRsid.get(m.rsid.toLowerCase())
    if (!hit) continue
    checked++
    if (hit.pos === m.pos) agreed++
    else offsets.push(hit.pos - m.pos)
  }
  // A handful of disagreements is ordinary (an indel's anchor base, a rebuilt rsID). A
  // majority disagreeing is an assembly difference.
  return { checked, agreed, offsets: offsets.slice(0, 5),
           mismatch: checked >= 4 && agreed / checked < 0.5 }
}

/** Per-side counts of markers this carrier is ACTUALLY heterozygous at: the panel's own
 *  coverage question asked again against a real person rather than a population prior. */
export function informativeCoverage(
  markers: { dist: number; rsid?: string | null; chrom: string; pos: number }[],
  set: GenotypeSet,
) {
  let lower = 0, higher = 0, hom = 0, absent = 0, nocall = 0
  for (const m of markers) {
    const c = lookup(m, set)
    if (!c) { absent++; continue }
    if (c.call === 'nocall') { nocall++; continue }
    if (c.call === 'hom') { hom++; continue }
    if (m.dist > 0) higher++
    else lower++
  }
  return { lower, higher, het: lower + higher, hom, absent, nocall }
}
