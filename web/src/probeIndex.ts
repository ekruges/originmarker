/**
 * A number for each probe id, so the per-marker tally reads and writes arrays instead of hashing the
 * same id into several Maps.
 *
 * The number is the order an id was first seen in this process; it is never written into a result.
 * One numbering serves every chip loaded in a session: it is a dictionary from id to number, not a
 * layout of a chip, so two chips share it without either reading the other's values. Every array
 * kept by number belongs to one Map object (the panel, its frequencies, one parent's genotypes), so
 * a second chip's panel or parent has arrays of its own.
 */
const NUMBER = new Map<string, number>()
const IDS: string[] = []

/** The probe's number, given the first time an id is seen. */
export function probeNumber(id: string): number {
  let k = NUMBER.get(id)
  if (k === undefined) { k = IDS.length; NUMBER.set(id, k); IDS.push(id) }
  return k
}

const NONE: unique symbol = Symbol('absent')
const tables = new WeakMap<ReadonlyMap<string, unknown>, unknown[]>()

/**
 * A reader of m by probe number: what m.get(id) returns, taken from m once a probe and kept. The Map
 * must not change after it is first read here: a parent's genotypes, the founder panel and its
 * frequencies are finished before any sample is scored against them.
 */
export function readerOf<V>(m: ReadonlyMap<string, V>): (k: number, id: string) => V | undefined {
  let a = tables.get(m) as (V | typeof NONE)[] | undefined
  if (!a) { a = []; tables.set(m, a) }
  const t = a
  return (k, id) => {
    const v = t[k]
    if (v === undefined) {
      const w = m.get(id)
      t[k] = w === undefined ? NONE : w
      return w
    }
    return v === NONE ? undefined : v as V
  }
}

/**
 * A Map from probe id that holds its values by probe number. It behaves as a Map: entries in the
 * order their ids were first set, a set of a present id keeps its place, and setAt(k, v) is set(id, v)
 * for the id numbered k, without hashing the id again.
 */
export class ProbeMap<V> implements Map<string, V> {
  private vals: V[] = []
  private at: number[] = []
  private order: number[] = []
  readonly [Symbol.toStringTag] = 'Map'
  get size(): number { return this.order.length }
  setAt(k: number, v: V): this {
    if (this.at[k] === undefined) { this.at[k] = this.order.length; this.order.push(k) }
    this.vals[k] = v
    return this
  }
  set(id: string, v: V): this { return this.setAt(probeNumber(id), v) }
  get(id: string): V | undefined {
    const k = NUMBER.get(id)
    return k === undefined || this.at[k] === undefined ? undefined : this.vals[k]
  }
  has(id: string): boolean {
    const k = NUMBER.get(id)
    return k !== undefined && this.at[k] !== undefined
  }
  delete(id: string): boolean {
    const k = NUMBER.get(id)
    if (k === undefined || this.at[k] === undefined) return false
    const i = this.at[k]
    this.order.splice(i, 1)
    for (let j = i; j < this.order.length; j++) this.at[this.order[j]] = j
    this.at[k] = undefined as unknown as number
    this.vals[k] = undefined as unknown as V
    return true
  }
  clear(): void { this.vals = []; this.at = []; this.order = [] }
  /** The value for probe number k, as get gives it for that number's id. */
  getAt(k: number): V | undefined { return this.at[k] === undefined ? undefined : this.vals[k] }
  /** The entries with their probe numbers, in the order entries gives them. */
  *numbered(): Generator<[number, string, V]> { for (const k of this.order) yield [k, IDS[k], this.vals[k]] }
  *entries(): MapIterator<[string, V]> { for (const k of this.order) yield [IDS[k], this.vals[k]] }
  *keys(): MapIterator<string> { for (const k of this.order) yield IDS[k] }
  *values(): MapIterator<V> { for (const k of this.order) yield this.vals[k] }
  [Symbol.iterator](): MapIterator<[string, V]> { return this.entries() }
  forEach(f: (v: V, id: string, m: Map<string, V>) => void, self?: unknown): void {
    for (const k of this.order) f.call(self, this.vals[k], IDS[k], this)
  }
}
