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
 * With --ran <stage>=<called.tsv> (audit/gate/coverage.py, from a full run of that stage on the base),
 * the stage's closures leave out functions that run never called (closure.ts), so a change only to
 * code the stage does not run is planned as nothing.
 *
 * Run: node --experimental-strip-types audit/gate/plan.ts [--ran <stage>=<called.tsv>]... <base> <head>
 *      node --experimental-strip-types audit/gate/plan.ts --keys [<commit>]   GATE-KEYS.json for a tree
 */
import { readFileSync } from 'node:fs'
import { closure, commitTree, declarations, worktree, type Tree } from './closure.ts'

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
export function keysOf(tree: Tree, stages = STAGES, ran: Record<string, Set<string>> = {}): Keys {
  const out: Keys = {}
  for (const s of stages) {
    if (!tree.read(s.entry)) { out[s.name] = null; continue }
    const all = closure(s.entry, [s.entry], tree, ran[s.name])
    const score = s.score && tree.read(s.score) ? closure(s.entry, [s.score], tree, ran[s.name]) : undefined
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

/** The names of the functions a coverage list says were called, read in the tree it was made on. */
export function ranNames(called: string, tree: Tree): Set<string> {
  const out = new Set<string>()
  const decls = new Map<string, ReturnType<typeof declarations>>()
  for (const line of readFileSync(called, 'utf8').split('\n')) {
    const [rel, a, b] = line.split('\t')
    if (!rel || !tree.read(rel)) continue
    if (!decls.has(rel)) decls.set(rel, declarations(rel, tree))
    const d = decls.get(rel)!.find((x) => x.start <= Number(a) && Number(b) <= x.end)
    if (d?.fn) for (const n of d.names) out.add(`${rel}#${n}`)
  }
  return out
}

export function plan(baseTree: Tree, headTree: Tree, stages = STAGES, ran: Record<string, Set<string>> = {}) {
  const b = keysOf(baseTree, stages, ran), h = keysOf(headTree, stages, ran)
  return stages.map((s) => ({ stage: s.name, decision: decide(b[s.name], h[s.name]) }))
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const args = process.argv.slice(2)
  const ranArgs: string[] = []
  while (args[0] === '--ran') { ranArgs.push(args[1]); args.splice(0, 2) }
  const [a, b] = args
  if (a === '--keys') {
    process.stdout.write(JSON.stringify(keysOf(b ? commitTree(b) : worktree()), null, 1) + '\n')
  } else if (a && b) {
    const baseTree = commitTree(a)
    const ran = Object.fromEntries(ranArgs.map((x) => { const [st, f] = x.split('='); return [st, ranNames(f, baseTree)] }))
    for (const { stage, decision } of plan(baseTree, commitTree(b), STAGES, ran)) if (decision) console.log(`${stage}\t${decision}`)
    console.log('checks\trun')
  } else {
    console.error('usage: plan.ts <base> <head> | plan.ts --keys [<commit>]')
    process.exit(2)
  }
}
