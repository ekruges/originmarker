/**
 * Which stages of the gate a change must run, from the code alone.
 *
 * Each stage is a script; its records depend on what that script can reach (closure.ts). A stage
 * with a scoring part (the real rows: score-array.ts) has two closures, its scoring and its whole.
 * Between a base and a head commit:
 *   score     the scoring closure changed: every array is scored again (the store keeps the results);
 *   surfaces  only the whole changed: records are rebuilt from the stored results;
 *   run       a stage with no scoring part whose closure changed, or a stage the base does not have;
 *   nothing   neither changed: the base's records stand.
 * The checks (scripts/checks.sh, scripts/release-check.sh) run on every change.
 *
 * Run: node --experimental-strip-types audit/gate/plan.ts <base> <head>
 *      node --experimental-strip-types audit/gate/plan.ts --keys [<commit>]   GATE-KEYS.json for a tree
 */
import { closure, commitTree, worktree, type Tree } from './closure.ts'

export interface Stage { name: string; entry: string; score?: string; jobs: string }
export const STAGES: Stage[] = [
  { name: 'rows', entry: 'audit/truth/run.ts', score: 'audit/truth/score-array.ts', jobs: 'audit/truth/run.ts' },
  { name: 'built', entry: 'audit/validate/known-answer-run.ts', jobs: 'audit/validate/known-answer-run.ts' },
  { name: 'suite-census', entry: 'audit/refusal-census.ts', jobs: 'audit/refusal-census.ts' },
  { name: 'suite-invariance', entry: 'audit/validate/invariance.ts', jobs: 'audit/validate/invariance.ts' },
  { name: 'suite-gains', entry: 'audit/validate/constructed-gains.ts', jobs: 'audit/validate/constructed-gains.ts' },
]

export type Keys = Record<string, { all: { hash: string; files: Record<string, string> }; score?: { hash: string; files: Record<string, string> } } | null>

/** Each stage's closures in a tree; null where the tree has no such script. */
export function keysOf(tree: Tree, stages = STAGES): Keys {
  const out: Keys = {}
  for (const s of stages) {
    if (!tree.read(s.entry)) { out[s.name] = null; continue }
    const all = closure(s.entry, [s.entry], tree)
    const score = s.score && tree.read(s.score) ? closure(s.entry, [s.score], tree) : undefined
    out[s.name] = { all: { hash: all.hash, files: all.files }, ...(score ? { score: { hash: score.hash, files: score.files } } : {}) }
  }
  return out
}

export type Decision = 'score' | 'surfaces' | 'run' | 'nothing'
export function decide(base: Keys[string], head: Keys[string]): Decision | null {
  if (!head) return null
  if (!base) return head.score ? 'score' : 'run'
  if (head.score && head.score.hash !== base.score?.hash) return 'score'
  if (head.all.hash !== base.all.hash) return head.score ? 'surfaces' : 'run'
  return 'nothing'
}

export function plan(baseTree: Tree, headTree: Tree, stages = STAGES) {
  const b = keysOf(baseTree, stages), h = keysOf(headTree, stages)
  return stages.map((s) => ({ stage: s.name, decision: decide(b[s.name], h[s.name]) }))
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const [a, b] = process.argv.slice(2)
  if (a === '--keys') {
    process.stdout.write(JSON.stringify(keysOf(b ? commitTree(b) : worktree()), null, 1) + '\n')
  } else if (a && b) {
    for (const { stage, decision } of plan(commitTree(a), commitTree(b))) if (decision) console.log(`${stage}\t${decision}`)
    console.log('checks\trun')
  } else {
    console.error('usage: plan.ts <base> <head> | plan.ts --keys [<commit>]')
    process.exit(2)
  }
}
