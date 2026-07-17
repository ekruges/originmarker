import type { LogLine } from './api'

/** What the exported file's footer identifies: which build, which instance, which moment.
 *  Only fields that live in the app; the browser-side ones (time, url, agent) are the ctx. */
export interface LogMeta {
  release?: string
  jobId?: string | null
  // Null-tolerant so the whole api.Provenance object is structurally assignable here; the
  // composer guards every field with `!= null` regardless.
  provenance?: {
    build?: string
    ensembl_release?: string | number | null
    queried_utc?: string
    window_bp?: number
    common_maf?: number
    candidate_n?: number
    elapsed_s?: number
    sources?: { gnomad?: string; genetic_map?: string }
  } | null
}

/** The three facts only the browser knows at the moment of export. Passed in rather than
 *  read inside, so the composer is pure and the self-check can pin the footer. */
export interface ExportCtx {
  now: string
  url: string
  agent: string
}

/**
 * The whole log as a .txt: every line, then one footer line carrying all the build/debug
 * info available, so a log pasted into a bug report says which deploy, which panel, and
 * which moment it came from. The footer is `#`-prefixed so it reads as a comment, and it is
 * one line by design: whoever receives the file greps it, they do not read it.
 *
 * Pure and JSX-free on purpose: BuildLog.check.ts runs under node, which cannot parse the
 * component's JSX, so the part worth pinning lives here where a fast test can reach it.
 */
export function buildLogText(lines: LogLine[], meta: LogMeta, ctx: ExportCtx): string {
  const body = lines.map((l) => `[${l.tag}] ${l.text}`).join('\n')
  const p = meta.provenance
  const footer = [
    'originmarker build/debug',
    `release=${meta.release ?? 'unknown'}`,
    `job=${meta.jobId ?? 'none'}`,
    `exported=${ctx.now}`,
    `url=${ctx.url}`,
    p?.build ? `genome=${p.build}` : null,
    p?.sources?.gnomad ? `gnomad=${p.sources.gnomad}` : null,
    p?.ensembl_release != null ? `ensembl=${p.ensembl_release}` : null,
    p?.sources?.genetic_map ? `map=${p.sources.genetic_map}` : null,
    p?.queried_utc ? `queried=${p.queried_utc}` : null,
    p?.window_bp != null ? `window_bp=${p.window_bp}` : null,
    p?.common_maf != null ? `maf=${p.common_maf}` : null,
    p?.candidate_n != null ? `candidates=${p.candidate_n}` : null,
    p?.elapsed_s != null ? `built_s=${p.elapsed_s}` : null,
    `agent=${ctx.agent}`,
  ].filter(Boolean).join(' | ')
  return `${body}\n\n# ${footer}\n`
}
