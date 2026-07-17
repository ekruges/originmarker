// Types mirror the fixed API contract + panelbuilder dataclasses exactly.
// Relative base: the app is served from /originmarker/, never hardcode an origin.
const BASE = 'api'

export interface Health {
  ok: boolean
  version: string
  // Release identity, e.g. 'Build 2.1 "Synapsis"'. Not `build` below, which is the GENOME
  // build and labels coordinates.
  release: string
  release_codename: string
  release_gloss: string
  gnomad_dataset: string
  build: string
  // null until the background warm-up lands. Render as '-', never as the string 'null'.
  ensembl_release: number | null
  map_source: string
  ldlink_enabled: boolean
  nl_enabled: boolean
  /** Two flags, because they fail apart: no primer3 means no pair at all, while no UCSC key
   *  means a pair that exists and is NOT VERIFIED. Never collapse them, or the UI implies
   *  the second was checked. Optional: a server older than the primer module sends neither,
   *  and absent must read as off rather than as unknown-so-probably-on. */
  primers_enabled?: boolean
  insilico_pcr_enabled?: boolean
  /** The settings this server designs under, i.e. asdict(primers.DEFAULTS). Optional, and
   *  the primer form only draws where it is present: the numbers must come from the engine
   *  that will use them, and a copy kept here would be a second set of defaults to drift. */
  primer_defaults?: PrimerSettings | null
  // Canonical text, served so the UI never paraphrases it.
  disclaimer: string
  layer_b_steps: string[]
}

export interface VariantRecord {
  query: string
  rsid: string | null
  gene: string | null
  strand: number | null
  chrom: string
  pos_grch38: number
  vcf_ref: string
  vcf_alt: string
  clinical_significance: string | null
  review_status: string | null
  clinvar_accession: string | null
  build: string
  pos_grch37: number | null
  build_note: string | null
}

export interface Rarity {
  gnomad_af_genome: number | null
  gnomad_ac_genome: number | null
  gnomad_an_genome: number | null
  thousand_genomes_ac: number | null
  population_LD_usable: boolean
  reason: string
}

// Engine emits these verbatim (panelbuilder._tier / annotate); do not shorten.
export const TIERS = ['A_core(<2kb)', 'B_near(2-30kb)', 'C_flank(30kb+)'] as const
export type Tier = (typeof TIERS)[number]
export const TIER_LABEL: Record<Tier, string> = {
  'A_core(<2kb)': 'A · core <2 kb',
  'B_near(2-30kb)': 'B · near 2–30 kb',
  'C_flank(30kb+)': 'C · flank 30 kb+',
}

export interface Marker {
  rsid: string
  variant_id: string
  chrom: string
  pos: number
  ref: string
  alt: string
  af: number
  maf: number
  /** Global expected heterozygosity 2pq: a population prior, never a carrier's genotype. */
  het: number
  het_max_pop: number
  /** Signed bp from the variant: negative = lower GRCh38 coordinate, positive = higher. */
  dist: number
  /** "lower coord" | "higher coord". Not tel/cen: that mapping depends on the centromere,
   *  and nothing here knows where it is. */
  side: string
  tier: Tier
  per_pop_maf: Record<string, number>
  ensembl_pos_check: string | null
  cm: number | null
  recomb_fraction: number | null
  hotspot_between: boolean | null
  map_approx: boolean | null
  /** The engine's verdict that this marker meets ESHRE's structural flanking criteria
   *  (panelbuilder.FLANKING_CRITERIA). Never recomputed here: the engine applies the rule
   *  once, and a mirror is a place to drift. Absent on a panel built before the rule
   *  shipped, which is silence, not false. */
  meets_eshre_flanking_criteria?: boolean
  /** The candidate genotyping pair for this marker, or the design's own account of why
   *  there is none. Absent means NO DESIGN WAS ATTEMPTED: only markers in scope are
   *  designed for. A design that ran and failed is a PrimerResult carrying `error`.
   *  Absence and failure are different facts and must never render alike. */
  primer?: PrimerResult | null
}

export const isUpper = (m: Marker) => m.dist > 0

// --- primers ------------------------------------------------------------------
// Mirrors primers.py: PrimerSettings, Primer and PrimerResult, name for name. The engine
// owns every knob, every range and every word of every warning; nothing here restates one.

/** Every knob on primers.PrimerSettings, in its names. The whole set, not the subset the
 *  form draws: `params.primer_settings` is an asdict() of all of them, and a rebuild hands
 *  the whole object back so a knob this UI never draws is carried rather than reset. */
export const PRIMER_SETTING_KEYS = [
  'opt_tm', 'min_tm', 'max_tm', 'max_pair_diff_tm',
  'min_size', 'opt_size', 'max_size',
  'min_gc', 'max_gc', 'gc_clamp',
  'salt_monovalent', 'salt_divalent', 'dntp_conc', 'dna_conc',
  'min_product', 'max_product',
  'max_poly_x', 'max_self_any_th', 'max_self_end_th', 'max_hairpin_th',
  'max_ns_accepted', 'tm_formula', 'salt_corrections',
  'mask_maf', 'target_pad',
] as const
export type PrimerSettingKey = (typeof PRIMER_SETTING_KEYS)[number]

/** primers.PrimerSettings. Every field is a number; the engine validates the ranges and
 *  words its own refusals, so nothing here re-judges a value it is about to send. */
export type PrimerSettings = Record<PrimerSettingKey, number>

/** panelbuilder.PRIMER_SCOPES. "starred" met the flanking criteria, "recommended" is the
 *  whole shortlist, "none" designs nothing. */
export const PRIMER_SCOPES = ['starred', 'recommended', 'none'] as const
export type PrimerScope = (typeof PRIMER_SCOPES)[number]

/** primers.SIZE_CAP. primer3's own cap sits above this, so a longer oligo is accepted and
 *  returns a Tm from outside the range its model is defined over: a plausible number, which
 *  is worse here than an error. */
export const PRIMER_SIZE_CAP = 35

/** primers.Primer. `seq` is what gets ordered, 5'->3' on the strand it binds from (R7);
 *  `pos` is 1-based GRCh38 of its leftmost template base, plus strand (R6). */
export interface Primer {
  seq: string
  pos: number
  idx: number
  length: number
  tm: number
  gc: number
}

/**
 * primers.PrimerResult: a pair with its warnings welded on, or a failure that says why.
 *
 * `fwd`/`rev` are null whenever `error` is set, and `warnings` is never empty for a pair.
 * That is the engine's invariant and this UI keeps its half of it: a primer never renders
 * without what is wrong with it.
 */
export interface PrimerResult {
  fwd?: Primer | null
  rev?: Primer | null
  product_size?: number | null
  product_start?: number | null
  /** What was kept clear, as primers.MaskSite. Rendered as a count beside `mask_note`,
   *  which is the engine's own statement of the mask's reach. */
  masked?: unknown[]
  /** The mask's reach, at both lengths. The counts differ per marker; the sentences after
   *  them are identical on every row, so only `short` belongs in the table. */
  mask_note?: PrimerWarning | null
  /** The engine's words about this pair, verbatim, at both lengths. Never empty for a pair,
   *  so an empty one is a contract this UI cannot vouch for and must not paper over. */
  warnings?: PrimerWarning[]
  /** A STATE CODE, never prose: `not_checked`, `one_product`, or one of app/ispcr.py's
   *  other verdicts. Absent and null both mean never checked, which is not a pass. The
   *  words that go with it are in `warnings`; rendering this field is rendering "danger". */
  insilico_pcr?: string | null
  error?: string | null
}

/**
 * primers.Note: one thing wrong with a pair, at both lengths, with somewhere to read more.
 *
 * `short` is what the table shows and `long` is what the exports print. Both are the
 * engine's; nothing here writes either, and nothing here shortens `long` on its own.
 */
export interface PrimerWarning {
  code: string
  short: string
  long: string
  docs: string
}

export const PCR_NOT_CHECKED = 'not_checked'
export const PCR_ONE_PRODUCT = 'one_product'
/** Asked, and the answer could not be read: a timeout, an unreadable page, a spent quota.
 *  app/ispcr.py keeps this apart from `danger` because the remedies differ, and so must
 *  this: one says redesign the pair, the other says ask again later. */
export const PCR_UNKNOWN = 'unknown'

/**
 * Is this pair's verdict one to shout about?
 *
 * `unknown` is NOT: UCSC never answered, so it contradicted nothing, and a red "do not
 * order this pair without redesigning it" over a timeout is a finding this page invented.
 * It is not clean either, which `pcrUnchecked` is what says so. The two predicates move
 * together or a quota stop renders green.
 *
 * Anything else, a state this UI has never heard of included, IS dangerous: an unreadable
 * verdict must never render as a clean one, because that is the one mistake that turns an
 * unverified pair into a verified one.
 */
export const pcrDangerous = (d: PrimerResult) => {
  const s = d.insilico_pcr
  return !(s == null || s === PCR_NOT_CHECKED || s === PCR_ONE_PRODUCT || s === PCR_UNKNOWN)
}

/**
 * The verdict's own words, or null where there is nothing to shout.
 *
 * jobs.py welds the verdict onto the pair as the warning whose code IS the state, so this
 * finds it by that code rather than by position. The fallback is not decoration: a
 * dangerous state whose words did not arrive still has to shout, and a silent red row
 * would be read as a clean one.
 */
export const pcrDanger = (d: PrimerResult): PrimerWarning | null => {
  if (!pcrDangerous(d)) return null
  const w = (d.warnings ?? []).find((x) => x.code === d.insilico_pcr)
  if (w) return w
  const said = `The verdict is "${d.insilico_pcr}" and its explanation did not reach this `
    + `page. Treat the pair as unverified and check it at UCSC before ordering.`
  return { code: d.insilico_pcr ?? 'unknown', short: said, long: said, docs: PRIMER_DOCS }
}

/** primers.PRIMER_DOCS. Where a note points when it carries no route of its own. */
export const PRIMER_DOCS = '#/docs/primers'

/**
 * Not checked against the genome, which is not a pass: no product anywhere else in GRCh38
 * has been ruled out.
 *
 * `unknown` belongs here, and this is the half that is easy to forget. It is not dangerous,
 * so if only `pcrDangerous` learns the token, a pair UCSC timed out on falls through both
 * predicates and renders as the green "one product" it never earned. app/ispcr.py's rule is
 * that only a parsed single product is a pass; everything else is some flavour of
 * unverified, and this is where that lands.
 */
export const pcrUnchecked = (d: PrimerResult) =>
  d.insilico_pcr == null
  || d.insilico_pcr === PCR_NOT_CHECKED
  || d.insilico_pcr === PCR_UNKNOWN

/** A design ran and produced a pair. `error` and a pair are mutually exclusive upstream;
 *  read the pair, not the absence of an error, so a result carrying both still warns. */
export const hasPair = (d: PrimerResult) => !!(d.fwd || d.rev)

/** The Marker field the star rule writes to, as panelbuilder.FLANKING_CRITERIA["field"]
 *  names it. Read through this so `flankingRule`'s name check governs the actual access. */
const STAR_FIELD = 'meets_eshre_flanking_criteria'

/** The engine judged this marker and said yes. Anything else, absent included, is not a
 *  star: only `true` may print one. */
export const starred = (m: Marker) => m[STAR_FIELD] === true

/** The engine's own statement of the star rule: its field, its numbers and its words. The
 *  UI renders `legend` verbatim and never restates the rule. */
export interface FlankingCriteria {
  field: string
  max_dist_bp: number
  min_per_side: number
  /** The key beside the star. */
  legend: string
  /** One line of criteria, for the hover. */
  summary?: string
  /** Docs route for a reader who clicks the star. */
  docs_href?: string
  /** The full wording, for print. */
  note: string[]
}

/**
 * The star rule this panel was built under, or null if the star cannot be honestly drawn.
 *
 * Null on an older panel that has no rule, and on a rule whose verdict no longer lands in
 * the field read here: a renamed field would otherwise show as a legend above rows that
 * quietly lost their stars, which reads as "none qualified". Both are silence, and silence
 * must render as nothing at all rather than as a rule the page cannot show the results of.
 */
export const flankingRule = (p: Provenance): FlankingCriteria | null => {
  const fc = p.flanking_criteria
  return fc && fc.field === STAR_FIELD && fc.legend ? fc : null
}

/**
 * Two independent sources disagree about where this marker IS.
 *
 * Three distinct states, collapsing any two is an error: null is silence (never checked,
 * not a pass), 'ok' is agreement, anything else is "MISMATCH:<pos>" carrying the position
 * Ensembl claims.
 */
export const posMismatch = (m: Marker) =>
  m.ensembl_pos_check != null && m.ensembl_pos_check !== 'ok'

export interface Coverage {
  lower_count: number
  higher_count: number
  lower_core_near: number
  higher_core_near: number
  /** Shortlisted markers meeting the flanking criteria, per side. The engine also flags a
   *  side that falls under its own minimum into `flags`, so a reader is told without these
   *  being rendered. Absent on a panel built before the rule shipped. */
  lower_flanking_count?: number
  higher_flanking_count?: number
  flags: string[]
}

export interface Provenance {
  sources: { clinvar?: string; ensembl?: string; gnomad?: string; genetic_map?: string }
  build: string
  window_bp: number
  common_maf: number
  ancestry_rank: string | null
  candidate_n: number
  requested_build: string
  /** The star rule, in the engine's words. Read it through `flankingRule`, never directly:
   *  its presence alone does not mean the star is safe to draw. */
  flanking_criteria?: FlankingCriteria
  /** The Ensembl release THIS panel was built against. Optional: older panels carry none
   *  and must render as unknown. Never substitute the live server's release, which is a
   *  fact about now, not about this panel. */
  ensembl_release?: number | null
  /** Oldest source response used, i.e. how old the DATA is. The cache has no TTL, so this
   *  is not built_utc. */
  queried_utc: string
  built_utc?: string
  source_responses_from_cache?: number
  source_responses_from_network?: number
  elapsed_s: number
  disclaimer: string
  layer_b_steps: string[]
}

export interface PanelResult {
  variant: VariantRecord
  rarity: Rarity
  candidates: Marker[]
  recommended: Marker[]
  coverage: Coverage
  params: Record<string, unknown>
  provenance: Provenance
}

/** The primer scope and settings THIS panel was built under, as `params` records them. */
export interface PrimerBuild {
  scope: PrimerScope
  settings: PrimerSettings
}

/**
 * What this panel was built under, or null if it does not record it.
 *
 * Null on a panel from a server with no primer module, which renders as nothing at all
 * rather than as an empty form. A settings block missing a knob is null too: a form seeded
 * from a gap would show a number nobody chose and then rebuild on it.
 */
export const primerBuild = (r: PanelResult): PrimerBuild | null => {
  const scope = r.params?.primer_scope
  const s = r.params?.primer_settings as Partial<PrimerSettings> | undefined
  if (!s || !(PRIMER_SCOPES as readonly string[]).includes(scope as string)) return null
  return PRIMER_SETTING_KEYS.every((k) => typeof s[k] === 'number')
    ? { scope: scope as PrimerScope, settings: s as PrimerSettings }
    : null
}

export interface ResolveResponse {
  variant: VariantRecord
  rarity: Rarity
  transcript_sense: string
  clinvar_url: string
  // `ld_banner` is omitted on purpose: the server may still send it, but rarity has
  // exactly one verdict and it is `rarity.reason`. Do not add it back.
}

export interface StructuredQuery {
  variant: string
  gene?: string | null
  window_bp?: number
  build?: string
  ancestry?: string | null
  common_maf?: number
  cross_check?: boolean
  /** Which markers get a pair. Inputs, beside the other inputs. */
  primer_scope?: PrimerScope
  /** Knob overrides, as primers.PrimerSettings names them. Partial: omitted asks for the
   *  server's own defaults, and that is the only way to ask for them. A copy held on this
   *  side drifts, and the panel would then report the copy as what it was built under. */
  primer_settings?: Partial<PrimerSettings>
  /**
   * Run the UCSC check as part of this build rather than leaving it to the primer box.
   *
   * The wire shape here is app.main.PanelIn, not pb.StructuredQuery: this says what to do
   * after the build, so the server pops it before constructing the query. Costs about 15
   * seconds per designed pair and spends the same per-IP budget as the button.
   */
  verify_primers?: boolean
}

/**
 * The /api/panel body: a StructuredQuery plus how the variant got chosen.
 *
 * Separate from StructuredQuery, which mirrors the pb dataclass the engine runs on, and
 * that dataclass has no provenance field. This is the wire shape only. Note what it still
 * cannot carry: a coordinate. There is no field for one on either side of the hop.
 */
export interface PanelRequest extends StructuredQuery {
  nl_text?: string
  /** The model that named `variant` out of its own memory. null is a claim too: it says no
   *  model chose this variant. Never set it from anything but `NLResponse.used_llm`. */
  nl_model?: string | null
}

export interface NLResponse {
  query: StructuredQuery
  used_llm: boolean
  note: string
  /** The text the user typed. */
  text?: string
  /** The model that supplied `query.variant`. Server sends it only when `used_llm`. */
  model?: string | null
  /** Gene symbols the USER named. Absent or empty means they named none: that is silence,
   *  not agreement, and nothing may be inferred from it. */
  named_genes?: string[]
}

/**
 * Fold the parse provenance into the query /api/panel will receive.
 *
 * `used_llm` is the sole authority on whether a model chose the variant: a model name in
 * the response cannot promote a free regex parse into a model's choice, and its absence on
 * the model path cannot demote one. A query a human typed carries no provenance at all.
 */
export const withNlProvenance = (q: StructuredQuery, nl: NLResponse | null): PanelRequest =>
  nl ? { ...q, nl_text: nl.text, nl_model: nl.used_llm ? nl.model ?? null : null } : q

/**
 * The user named a gene, and the variant that resolved sits in a different one.
 *
 * Three states, and collapsing any two is an error: naming no gene is silence (there is
 * nothing to disagree with), a record with no gene of its own cannot disagree either, and
 * only a gene that is present on both sides and differs is a mismatch. Compared case-
 * insensitively; an alias the user typed reads as a mismatch, which is the safe direction.
 */
export const geneMismatch = (resolved: string | null, named?: string[]): boolean => {
  const want = (named ?? []).map((g) => (g ?? '').trim().toUpperCase()).filter(Boolean)
  const got = (resolved ?? '').trim().toUpperCase()
  return want.length > 0 && got !== '' && !want.includes(got)
}

// The tags panelbuilder emits, as panelbuilder.TAGS defines them: a line tagged with one
// the engine has and this list lacks renders as INFO, which reads as routine and hides
// what the tag was for. Anything genuinely unknown still renders, as INFO: a tag is a hint
// about a line, never a licence to drop one.
export const LOG_TAGS = ['FETCH', 'CACHE', 'INFO', 'WARN', 'SKIP', 'DONE'] as const
export type LogTag = (typeof LOG_TAGS)[number]

export interface LogLine {
  tag: LogTag
  text: string
}

/**
 * Normalise one build-log frame. Returns null for a frame that is not a log line.
 *
 * app/jobs.py owns the wire shape and the API passes it through unreshaped, so this is the
 * only place that assumes anything about it. A frame with no text is not a line (which is
 * what keeps a `progress` frame from being read as one); a missing or unknown tag is.
 */
export const asLogLine = (d: unknown): LogLine | null => {
  if (!d || typeof d !== 'object') return null
  const o = d as Record<string, unknown>
  const text = typeof o.text === 'string' ? o.text
    : typeof o.message === 'string' ? o.message : null
  if (text === null) return null
  const tag = typeof o.tag === 'string' ? o.tag.trim().toUpperCase() : ''
  return { tag: (LOG_TAGS as readonly string[]).includes(tag) ? (tag as LogTag) : 'INFO', text }
}

export interface JobStatus {
  status: 'running' | 'done' | 'error'
  result?: PanelResult
  error?: string
  /** The whole log so far, not a delta: the poller's copy of the SSE `log` stream, for a
   *  proxy that buffers SSE. Declared, not trusted: every line goes through asLogLine. */
  log?: LogLine[]
}

export interface LDResponse {
  r2: number
  dprime: number
  pop: string
  note: string
  caveat: string
}

export const ANCESTRIES = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'SAS', 'MID'] as const
export type Ancestry = (typeof ANCESTRIES)[number]

/** Carries the server's `detail` so the UI can fail loudly rather than guess. */
export class ApiError extends Error {
  status: number
  constructor(message: string, status: number) {
    super(message)
    this.name = 'ApiError'
    this.status = status
  }
}

async function req<T>(path: string, init?: RequestInit): Promise<T> {
  let res: Response
  try {
    res = await fetch(`${BASE}${path}`, {
      ...init,
      headers: { 'Content-Type': 'application/json', ...init?.headers },
    })
  } catch {
    throw new ApiError('Cannot reach the OriginMarker API. Is the backend running?', 0)
  }
  if (!res.ok) {
    let detail = `${res.status} ${res.statusText}`
    try {
      const body = await res.json()
      if (body?.detail) detail = typeof body.detail === 'string' ? body.detail : JSON.stringify(body.detail)
    } catch {
      /* non-JSON error body; keep the status line */
    }
    throw new ApiError(detail, res.status)
  }
  return res.json() as Promise<T>
}

export const api = {
  health: () => req<Health>('/health'),
  resolve: (variant: string, build?: string) =>
    req<ResolveResponse>('/resolve', { method: 'POST', body: JSON.stringify({ variant, build }) }),
  panel: (q: PanelRequest) => req<{ job_id: string }>('/panel', { method: 'POST', body: JSON.stringify(q) }),
  job: (id: string) => req<JobStatus>(`/panel/${encodeURIComponent(id)}`),
  nl: (text: string) => req<NLResponse>('/nl', { method: 'POST', body: JSON.stringify({ text }) }),
  genes: (q: string) => req<{ symbol: string; description: string }[]>(`/genes?q=${encodeURIComponent(q)}`),
  ld: (a: string, b: string, pop: string) =>
    req<LDResponse>(`/ld?a=${encodeURIComponent(a)}&b=${encodeURIComponent(b)}&pop=${encodeURIComponent(pop)}`),
  streamUrl: (id: string) => `${BASE}/panel/${encodeURIComponent(id)}/stream`,
  exportUrl: (id: string, ext: string) => `${BASE}/export/${encodeURIComponent(id)}.${ext}`,
  /** Check a finished panel's pairs against the whole genome. Never automatic: UCSC allows
   *  one request every 15 seconds, so this is a deliberate act with a real wait, and the
   *  panel is already usable without it. `rsids` omitted means every designed pair. */
  verify: (panelJobId: string, rsids?: string[]) =>
    req<{ job_id: string }>(`/panel/${encodeURIComponent(panelJobId)}/verify`,
      { method: 'POST', body: JSON.stringify({ rsids: rsids ?? null }) }),
  verifyJob: (id: string) => req<VerifyStatus>(`/verify/${encodeURIComponent(id)}`),
}

/** A verification run. `verdicts` are app/ispcr.py's own dicts, keyed by rsID: it words
 *  every verdict, and a second wording here would be the one that drifts into a pass. */
export interface VerifyStatus {
  status: 'running' | 'done' | 'error'
  stage: string
  fraction: number
  log?: unknown[]
  verdicts: Record<string, { state: string; note: string; products?: unknown[] }>
  error?: string | null
}

// --- outbound record links, built from rsid/coord (never from memory) ---------
export const links = {
  dbsnp: (rsid: string) => `https://www.ncbi.nlm.nih.gov/snp/${rsid}`,
  gnomad: (m: Pick<Marker, 'chrom' | 'pos' | 'ref' | 'alt'>) =>
    `https://gnomad.broadinstitute.org/variant/${m.chrom}-${m.pos}-${m.ref}-${m.alt}?dataset=gnomad_r4`,
  ensembl: (rsid: string) => `https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=${rsid}`,
  ucsc: (chrom: string, pos: number) =>
    `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr${chrom.replace(/^chr/, '')}:${pos - 100}-${pos + 100}`,
  clinvar: (acc: string) => `https://www.ncbi.nlm.nih.gov/clinvar/variation/${acc.replace(/^VCV0*/, '')}/`,
}
