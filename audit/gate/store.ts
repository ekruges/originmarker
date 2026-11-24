/**
 * A content-addressed store of scored results, shared by every lane and the baseline.
 *
 * A key is a sha256 over the scoring closure's hash (audit/gate/closure.ts, read from the tree's
 * GATE-KEYS.json), the node that runs it, and the bytes and names of every input file, with the
 * other arguments the scoring takes. It is a hash, never a file name. A value is the scored objects,
 * serialized whole by node:v8 (Maps, typed arrays, NaN and shared references survive), with long
 * arrays of integers held as deltas, gzipped, written to a temporary name and renamed.
 *
 * Inputs: an array is read only where OM_OPEN lists its path, one absolute path a line, written by
 * audit/gate/seal.py from units.tsv. A path not listed, a sealed or exposed array among them, is
 * refused before a byte of it is read, whether to score it, to key it or to store it.
 *
 * The closure's hash is trusted only for the tree it was computed on: every file GATE-KEYS.json
 * lists must have the bytes it names, or the store refuses to open.
 */
import { readFileSync, writeFileSync, renameSync, mkdirSync, existsSync } from 'node:fs'
import { createHash } from 'node:crypto'
import { gzipSync, gunzipSync } from 'node:zlib'
import { serialize, deserialize } from 'node:v8'
import { join } from 'node:path'

const sha = (b: Buffer | string) => createHash('sha256').update(b).digest('hex')

/** The scoring hash of a stage, checked against the bytes of the tree it runs from. */
export function treeKey(treeRoot: string, stage: string, part = 'score'): string {
  const f = join(treeRoot, 'GATE-KEYS.json')
  if (!existsSync(f)) throw new Error('GATE-KEYS.json missing: ship the tree with audit/gate/ship.sh')
  const k = JSON.parse(readFileSync(f, 'utf8'))[stage]?.[part]
  if (!k) throw new Error(`GATE-KEYS.json has no ${stage}.${part}`)
  for (const [p, h] of Object.entries(k.files as Record<string, string>)) {
    const q = join(treeRoot, p)
    if (!existsSync(q) || sha(readFileSync(q)) !== h) throw new Error(`GATE-KEYS.json does not match this tree at ${p}`)
  }
  return k.hash
}

/** The paths an open list allows. */
export function openList(file: string | undefined): Set<string> {
  if (!file || !existsSync(file)) throw new Error('OM_OPEN must name the open list (audit/gate/seal.py)')
  return new Set(readFileSync(file, 'utf8').split('\n').filter(Boolean))
}

const INT = '\u0000om-int32-deltas'
/** The differences of an array of integers, where each fits an int32; null otherwise. */
function deltas(x: unknown[]): Int32Array | null {
  const d = new Int32Array(x.length)
  let prev = 0
  for (let i = 0; i < x.length; i++) {
    const v = x[i]
    if (typeof v !== 'number' || !Number.isInteger(v) || Object.is(v, -0)) return null
    const e = v - prev
    if (e < -(2 ** 31) || e >= 2 ** 31) return null
    d[i] = e; prev = v
  }
  return d
}
/** A copy with every long array of int32 integers held as deltas; shared references stay shared. */
export function pack(x: unknown, seen = new Map<object, unknown>()): unknown {
  if (!x || typeof x !== 'object') return x
  if (seen.has(x)) return seen.get(x)
  if (Array.isArray(x)) {
    const d = x.length >= 64 && Object.keys(x).length === x.length ? deltas(x) : null
    if (d) { const o = { [INT]: d }; seen.set(x, o); return o }
    const o: unknown[] = new Array(x.length)
    seen.set(x, o)
    for (let i = 0; i < x.length; i++) if (i in x) o[i] = pack(x[i], seen)
    return o
  }
  if (x instanceof Map) {
    const o = new Map(); seen.set(x, o)
    for (const [k, v] of x) o.set(pack(k, seen), pack(v, seen))
    return o
  }
  if (x instanceof Set) {
    const o = new Set(); seen.set(x, o)
    for (const v of x) o.add(pack(v, seen))
    return o
  }
  if (ArrayBuffer.isView(x) || x instanceof Date || x instanceof RegExp || x instanceof Error) return x
  const o: Record<string, unknown> = {}
  seen.set(x, o)
  for (const k of Object.keys(x)) o[k] = pack((x as any)[k], seen)
  return o
}
export function unpack(x: unknown, seen = new Map<object, unknown>()): unknown {
  if (!x || typeof x !== 'object') return x
  if (seen.has(x)) return seen.get(x)
  if (Array.isArray(x)) {
    const o: unknown[] = new Array(x.length); seen.set(x, o)
    for (let i = 0; i < x.length; i++) if (i in x) o[i] = unpack(x[i], seen)
    return o
  }
  if (x instanceof Map) { const o = new Map(); seen.set(x, o); for (const [k, v] of x) o.set(unpack(k, seen), unpack(v, seen)); return o }
  if (x instanceof Set) { const o = new Set(); seen.set(x, o); for (const v of x) o.add(unpack(v, seen)); return o }
  if (ArrayBuffer.isView(x) || x instanceof Date || x instanceof RegExp || x instanceof Error) return x
  const d = (x as any)[INT]
  if (d instanceof Int32Array && Object.keys(x).length === 1) {
    const o = new Array(d.length); let s = 0
    for (let i = 0; i < d.length; i++) { s += d[i]; o[i] = s }
    seen.set(x, o)
    return o
  }
  const o: Record<string, unknown> = {}; seen.set(x, o)
  for (const k of Object.keys(x)) o[k] = unpack((x as any)[k], seen)
  return o
}

export interface Store {
  /** Whether every input is on the open list; key() refuses any that is not. */
  allows(inputs: string[]): boolean
  key(parts: unknown[], inputs: string[]): string
  get(key: string): any | undefined
  put(key: string, value: unknown): void
  readonly stats: { hits: number; misses: number; bytes: number }
}

export function openStore(dir: string, scoreHash: string, open: Set<string>): Store {
  const fileHash = new Map<string, string>()
  const allowed = (p: string) => { if (!open.has(p)) throw new Error('an input not on the open list is refused') }
  const where = (k: string) => join(dir, 'objects', k.slice(0, 2), `${k}.v8.gz`)
  const stats = { hits: 0, misses: 0, bytes: 0 }
  return {
    stats,
    allows: (inputs) => inputs.every((p) => p === '-' || p === '' || open.has(p)),
    key(parts, inputs) {
      const h = createHash('sha256').update(`om-store-1\n${scoreHash}\n${process.version}\n${JSON.stringify(parts)}\n`)
      for (const p of inputs) {
        if (p === '-' || p === '') { h.update(`${p}\n`); continue }
        allowed(p)
        let f = fileHash.get(p)
        if (!f) { f = sha(readFileSync(p)); fileHash.set(p, f) }
        h.update(`${p}\n${f}\n`)
      }
      return h.digest('hex')
    },
    get(key) {
      const f = where(key)
      if (!existsSync(f)) { stats.misses++; return undefined }
      stats.hits++
      return unpack(deserialize(gunzipSync(readFileSync(f))))
    },
    put(key, value) {
      const f = where(key)
      mkdirSync(join(dir, 'objects', key.slice(0, 2)), { recursive: true })
      const b = gzipSync(serialize(pack(value)), { level: 6 })
      const tmp = `${f}.${process.pid}.tmp`
      writeFileSync(tmp, b)
      renameSync(tmp, f)
      stats.bytes += b.length
    },
  }
}

/**
 * An object that answers only the fields it holds and throws on any other, so a record built from
 * a stored result cannot silently read a field the store did not keep. `lazy` fields are computed
 * on first read.
 */
export function only<T extends object>(held: T, lazy: Record<string, () => unknown> = {}): T {
  const memo = new Map<string, unknown>()
  return new Proxy(held, {
    get(t, k, r) {
      if (typeof k !== 'string' || k === 'then' || k in t) return Reflect.get(t, k, r)
      if (k in lazy) { if (!memo.has(k)) memo.set(k, lazy[k]()); return memo.get(k) }
      throw new Error(`a stored result does not hold ${k}`)
    },
  })
}
