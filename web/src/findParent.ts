/**
 * Which of several candidate fathers is the father of a run of arrays.
 *
 * Each candidate is read against each array by the parent check of parentGenome.ts, unchanged: the
 * same markers, the same bars, the same verdict. Nothing here has a threshold of its own.
 *
 * One array names a candidate when exactly one carries. A run names the parent when exactly one
 * candidate carries on any of its arrays. An array that carries none of them abstains: a maternal
 * pronucleus, a polar body, a cell that lost this parent's chromosomes and a noisy array all read
 * that way, and none of them speaks against any candidate. When two or more candidates carry, none
 * is named: the check does not tell a parent from a parent's close relative, and a cell carries its
 * mother's genome as surely as its father's.
 */
import { addCarry, emptyCarryTally, parentGenome, type CarryTally, type ParentGenome } from './parentGenome.ts'
import { isAutosome } from './parentage.ts'

/** Per candidate, per autosome, the parent check's tally against one array. */
export type CandidateTallies = Map<string, CarryTally>[]

export const emptyCandidateTallies = (candidates: number): CandidateTallies =>
  Array.from({ length: candidates }, () => new Map())

/**
 * One row of the array into every candidate's tally. The same markers scoreSample tallies for a
 * labelled parent: autosomal, with a founder frequency.
 */
export function addCandidateRow(
  t: CandidateTallies, candidates: readonly ReadonlyMap<string, string>[],
  r: { chrom: string; probesetId: string; genotype: string }, freqB: ReadonlyMap<string, number>,
): void {
  if (!isAutosome(r.chrom)) return
  const fq = freqB.get(r.probesetId)
  if (fq === undefined) return
  for (let i = 0; i < candidates.length; i += 1) {
    let m = t[i].get(r.chrom)
    if (!m) t[i].set(r.chrom, m = emptyCarryTally())
    addCarry(m, candidates[i].get(r.probesetId), r.genotype, fq)
  }
}

export const candidateGenomes = (t: CandidateTallies): ParentGenome[] =>
  t.map((m) => parentGenome([...m.values()]))

/** One array: the candidates whose genome it carries. */
export const carriers = (gs: readonly ParentGenome[]): number[] =>
  gs.flatMap((g, i) => (g.verdict === 'carries' ? [i] : []))

export interface ParentSearch {
  /** Per candidate, over the arrays read: how many carry it, read without it, or do not show it. */
  counts: { carries: number; none: number; unshown: number }[]
  arrays: number
  /** Arrays that carry no candidate, and arrays that carry two or more. */
  carryNone: number
  carrySeveral: number
  outcome: 'named' | 'none' | 'several'
  /** Index of the named candidate. */
  named: number | null
}

export function findParent(arrays: readonly (readonly ParentGenome[])[], candidates: number): ParentSearch {
  const counts = Array.from({ length: candidates }, () => ({ carries: 0, none: 0, unshown: 0 }))
  let carryNone = 0
  let carrySeveral = 0
  for (const gs of arrays) {
    gs.forEach((g, i) => {
      if (g.verdict === 'carries') counts[i].carries += 1
      else if (g.verdict === 'none') counts[i].none += 1
      else counts[i].unshown += 1
    })
    const c = carriers(gs).length
    if (c === 0) carryNone += 1
    if (c > 1) carrySeveral += 1
  }
  const carried = counts.flatMap((c, i) => (c.carries ? [i] : []))
  return {
    counts, arrays: arrays.length, carryNone, carrySeveral,
    outcome: carried.length === 1 ? 'named' : carried.length === 0 ? 'none' : 'several',
    named: carried.length === 1 ? carried[0] : null,
  }
}

/** The known limit, as the tool states it. */
export const FATHER_SEARCH_LIMIT = 'The check cannot tell a parent from a close relative of that '
  + 'parent. Two candidates who are close relatives of each other, such as brothers or a father and '
  + 'his son, can both carry, and then neither is named. A cell carries its mother\'s genome too, so '
  + 'label her file oocyte, not sperm.'

export interface ParentSearchText {
  rows: { name: string; carries: number; none: number; unshown: number }[]
  line: string
  abstain: string
  limit: string
}

/** What the page and the command line print. `names` are the candidate fathers' file names, in order. */
export function fatherSearchText(s: ParentSearch, names: readonly string[]): ParentSearchText {
  const rows = names.map((name, i) => ({ name, ...s.counts[i] }))
  const of = `${s.arrays} array${s.arrays === 1 ? '' : 's'}`
  const line = s.outcome === 'named'
    ? `${names[s.named!]} is the father: its genome is in ${rows[s.named!].carries} of ${of}, and no `
      + 'other candidate\'s is in any. The samples are scored against it.'
    : s.outcome === 'none'
      ? `None of the ${names.length} candidates is shown to be the father: no array carries any of `
        + 'their genomes. The samples are scored with no father loaded.'
      : `No father is named: ${rows.filter((r) => r.carries).map((r) => `${r.name} is in ${r.carries}`)
        .join(', ')} of ${of}. The samples are scored with no father loaded. To score against one `
        + 'of them, leave only that file labelled sperm.'
  const abstain = s.carryNone
    ? `${s.carryNone} of ${of} carry none of the candidates and count for or against none: a `
      + 'maternal pronucleus, a polar body, a cell that lost the father\'s chromosomes, or an array '
      + 'whose calls are noise reads that way.'
    : ''
  return { rows, line, abstain, limit: FATHER_SEARCH_LIMIT }
}
