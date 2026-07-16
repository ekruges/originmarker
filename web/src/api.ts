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
}

export const isUpper = (m: Marker) => m.dist > 0

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

export interface JobStatus {
  status: 'running' | 'done' | 'error'
  result?: PanelResult
  error?: string
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
