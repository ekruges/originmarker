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
