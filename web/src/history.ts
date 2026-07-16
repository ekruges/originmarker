// Local query history: what was typed, what it resolved to, how many candidates it found.
//
// This app has no account and keeps no server-side record of who asked what. This file is
// the only thing that remembers a query at all, it remembers it in one browser, and every
// entry has to be erasable from the UI that shows it. It must not become the first place
// the app accumulates a record of someone's work.

const KEY = 'originmarker.history'

/** Enough to recognise a session's work, few enough to read without scrolling. */
export const HISTORY_CAP = 12

export interface Entry {
  /** Verbatim, as typed: it is what the user will recognise in the list. */
  text: string
  /** What the parse resolved `text` to, when it resolved to something. */
  variant?: string
  /** Candidates in the panel built from this query. Absent until a panel is built, and
   *  absence is a state the row must render: a query that never built still belongs here. */
  count?: number
  /** ms since the epoch. */
  at: number
}

/** localStorage, or nothing. Reading `window.localStorage` throws outright where the origin
 *  refuses storage, so even the lookup is guarded: losing the history is acceptable, losing
 *  the search box is not. */
function store(): Storage | null {
  try {
    return window.localStorage
  } catch {
    return null
  }
}

/** One stored record, or nothing. Anything another version wrote is garbage to this one,
 *  and garbage is dropped rather than rendered. */
function clean(v: unknown): Entry[] {
  if (!v || typeof v !== 'object') return []
  const o = v as Record<string, unknown>
  const text = typeof o.text === 'string' ? o.text.trim() : ''
  if (!text) return []
  const e: Entry = { text, at: typeof o.at === 'number' && Number.isFinite(o.at) ? o.at : 0 }
  if (typeof o.variant === 'string' && o.variant.trim()) e.variant = o.variant.trim()
  // A count is a number of markers, so a negative or fractional one did not come from a
  // panel and must not be shown as one.
  if (typeof o.count === 'number' && Number.isInteger(o.count) && o.count >= 0) e.count = o.count
  return [e]
}

/** Every stored query, most recent first. Never throws: unreadable storage is an empty
 *  history, which is a state the UI already has to handle. */
export function readHistory(): Entry[] {
  const s = store()
  if (!s) return []
  let raw: string | null
  try {
    raw = s.getItem(KEY)
  } catch {
    return []
  }
  if (!raw) return []
  try {
    const parsed: unknown = JSON.parse(raw)
    return Array.isArray(parsed) ? parsed.flatMap(clean).slice(0, HISTORY_CAP) : []
  } catch {
    return []
  }
}

/** Never throws: storage can be full or refused. Returns what was handed in, so the caller
 *  renders what it meant to store whether or not the store took it. */
function write(es: Entry[]): Entry[] {
  try {
    store()?.setItem(KEY, JSON.stringify(es))
  } catch {
    /* quota or refused: the list is still correct for this page's lifetime */
  }
  return es
}

const same = (e: Entry, text: string) => e.text === text

/**
 * Note a submitted query. One entry per distinct text, most recent first.
 * Matched verbatim: HGVS is case-sensitive, so two texts differing in case are two queries.
 */
export function recordQuery(text: string, variant?: string): Entry[] {
  const t = text.trim()
  if (!t) return readHistory()
  const prev = readHistory()
  // The count carries over: it belongs to this text, and re-typing the text does not unknow
  // what the last panel built from it held. A fresh build overwrites it.
  const e: Entry = { ...prev.find((x) => same(x, t)), text: t, at: Date.now() }
  if (variant?.trim()) e.variant = variant.trim()
  return write([e, ...prev.filter((x) => !same(x, t))].slice(0, HISTORY_CAP))
}

/**
 * Record how many candidates a built panel held, against the query that asked for it.
 * Keyed on the resolved variant because a builder knows the variant, not the text someone
 * typed to reach it. A count for a variant no entry claims is dropped, never invented.
 */
export function noteCount(variant: string, count: number): Entry[] {
  const prev = readHistory()
  const v = variant.trim()
  if (!v || !Number.isInteger(count) || count < 0) return prev
  const i = prev.findIndex((e) => e.variant === v)
  if (i < 0) return prev
  const next = [...prev]
  next[i] = { ...next[i], count }
  return write(next)
}

/** Drop one query. */
export function forgetQuery(text: string): Entry[] {
  return write(readHistory().filter((e) => !same(e, text.trim())))
}

/** Drop the lot, and the key with it: an empty array left behind is still a record that
 *  someone was here. */
export function clearHistory(): Entry[] {
  try {
    store()?.removeItem(KEY)
  } catch {
    /* refused: nothing was readable through this module anyway */
  }
  return []
}
