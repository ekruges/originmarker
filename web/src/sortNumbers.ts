/**
 * Sorts an array of numbers ascending, in place, with the result `xs.sort((a, b) => a - b)` gives,
 * and returns it.
 *
 * Where every element is a number that is neither NaN nor -0, the values go through a Float64Array's
 * native numeric sort and back: for such values the comparator sort's order is fixed by value
 * alone, and a double round-trips a Float64Array exactly. Any other array (NaN, -0, or a
 * non-number, where the comparator sort's order also depends on where an element started) is
 * sorted by the comparator itself.
 */
export function sortNumbers<T extends number[]>(xs: T): T {
  const n = xs.length
  for (let i = 0; i < n; i++) {
    const v = xs[i]
    if (typeof v !== 'number' || v !== v || (v === 0 && 1 / v < 0) || !(i in xs)) return xs.sort((a, b) => a - b)
  }
  const t = Float64Array.from(xs).sort()
  for (let i = 0; i < n; i++) xs[i] = t[i]
  return xs
}

/** Whether every element is a number that is neither NaN nor -0, so order is fixed by value alone. */
function plain(xs: readonly unknown[]): boolean {
  for (let i = 0; i < xs.length; i++) {
    const v = xs[i]
    if (typeof v !== 'number' || v !== v || (v === 0 && 1 / v < 0) || !(i in xs)) return false
  }
  return true
}

/** Hoare selection in place: after it, t[k] is the value an ascending sort puts at k. */
function select(t: Float64Array, k: number): number {
  let lo = 0, hi = t.length - 1
  while (hi > lo) {
    const p = t[(lo + hi) >> 1]
    let i = lo, j = hi
    while (i <= j) {
      while (t[i] < p) i++
      while (t[j] > p) j--
      if (i <= j) { const s = t[i]; t[i] = t[j]; t[j] = s; i++; j-- }
    }
    if (k <= j) hi = j
    else if (k >= i) lo = i
    else return t[k]
  }
  return t[k]
}

/**
 * The value at index k of `[...xs].sort((a, b) => a - b)`, without sorting: by selection where every
 * element is a plain number (the value at a sorted position does not depend on how ties are
 * ordered), by that sort otherwise. xs is not changed.
 */
export function kthNumber(xs: readonly number[], k: number): number {
  if (!plain(xs)) return [...xs].sort((a, b) => a - b)[k]
  return select(Float64Array.from(xs), k)
}

/** The values at indexes k - 1 and k of the ascending sort, as kthNumber reads them. */
export function kthPair(xs: readonly number[], k: number): [number, number] {
  if (!plain(xs)) { const q = [...xs].sort((a, b) => a - b); return [q[k - 1], q[k]] }
  const t = Float64Array.from(xs)
  const hi = select(t, k)
  let lo = -Infinity
  for (let i = 0; i < k; i++) if (t[i] > lo) lo = t[i]
  return [lo, hi]
}
