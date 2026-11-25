/**
 * The scoring closure: what changes its hash and what does not, on a small program held in memory,
 * and on the repo's own scoring.
 *
 * Run: node --experimental-strip-types audit/gate/closure.check.ts
 */
import assert from 'node:assert/strict'
import { closure, worktree, type Tree } from './closure.ts'
import { decide } from './plan.ts'

const mem = (files: Record<string, string>): Tree => ({
  read: (rel) => (rel in files ? Buffer.from(files[rel]) : null),
  list: (dir) => Object.keys(files).filter((p) => p.startsWith(`${dir}/`)),
})
const BASE: Record<string, string> = {
  'audit/e.ts': `import { f } from '../web/src/a.ts'
const W = new URL('../web/src/', import.meta.url).pathname
const N = await import(\`\${W}b.ts\`)
console.log(f(1), N.x())`,
  'web/src/a.ts': `import { PATH } from './c.ts'
let count = 0
/** f's doc */
export function f(n: number): number { count++; return g(n) + PATH.length }
function g(n: number) { return n * 2 }
export function h() { return 'wording' }
export function w() { count = 5 }
export const T = { a: 1 }
export function readsT() { return T.a }`,
  'web/src/b.ts': `export const x = () => 1
export const y = () => 'other wording'
export interface I { p: number }`,
  'web/src/c.ts': `export const PATH = 'panel/d.txt'`,
  'web/public/panel/d.txt': 'panel v1',
  'web/public/other.txt': 'not named',
}
const hashOf = (over: Record<string, string> = {}) => closure('audit/e.ts', ['audit/e.ts'], mem({ ...BASE, ...over })).hash
const h0 = hashOf()
const edit = (f: string, a: string, b: string) => { assert.ok(BASE[f].includes(a), a); return { [f]: BASE[f].replace(a, b) } }

// A declaration nothing reached names does not change the hash; one that is reached does.
assert.equal(hashOf(edit('web/src/a.ts', `'wording'`, `'new wording'`)), h0, 'an unreached function')
assert.notEqual(hashOf(edit('web/src/a.ts', 'n * 2', 'n * 3')), h0, 'a reached function')
// Comments, layout and types do not; module state written by an unreached function does.
assert.equal(hashOf(edit('web/src/a.ts', `/** f's doc */`, `/** f's new doc */`)), h0, 'a comment')
assert.equal(hashOf(edit('web/src/a.ts', 'g(n: number)', 'g(n: any)')), h0, 'a type annotation')
assert.equal(hashOf(edit('web/src/b.ts', 'p: number', 'p: string')), h0, 'a type-only declaration')
assert.equal(hashOf(edit('web/src/a.ts', 'return g(n) + PATH.length', 'return  g( n )  +  PATH.length')), h0, 'layout')
assert.notEqual(hashOf(edit('web/src/a.ts', 'count = 5', 'count = 6')), h0, 'a writer of reached module state')
// A namespace reached by member follows only the members used.
assert.equal(hashOf(edit('web/src/b.ts', `'other wording'`, `'changed'`)), h0, 'an unused namespace member')
assert.notEqual(hashOf(edit('web/src/b.ts', '() => 1', '() => 2')), h0, 'a used namespace member')
// Code that runs on load in any module of the program is reached.
assert.notEqual(hashOf({ 'web/src/b.ts': BASE['web/src/b.ts'] + '\nglobalThis.k = 1' }), h0, 'a load-time statement')
// Shipped data named by a reached string enters by its bytes; data nothing names does not.
assert.notEqual(hashOf({ 'web/public/panel/d.txt': 'panel v2' }), h0, 'named data')
assert.equal(hashOf({ 'web/public/other.txt': 'changed' }), h0, 'unnamed data')
// A path taken from the environment brings in every shipped file.
const envd = { 'audit/e.ts': BASE['audit/e.ts'] + '\nconst p = process.env.OM_PANEL_FILE' }
assert.notEqual(hashOf({ ...envd, 'web/public/other.txt': 'changed' }), hashOf(envd), 'an environment path')
// An import() the closure cannot follow is an error, never a short closure.
assert.throws(() => hashOf({ 'audit/e.ts': 'const q = "b.ts"\nawait import(q)' }), /cannot resolve/)
assert.throws(() => hashOf({ 'audit/e.ts': `import { f } from '../web/src/gone.ts'` }), /not found/)

// A run's called functions: one it named but never called does not count; one it called does; a
// value it reads always does.
{
  const prog = { ...BASE, 'web/src/a.ts': BASE['web/src/a.ts'] + `
export function rare(n: number) { return n > 1e9 ? h() : 0 }
export const LIMIT = 3` }
  prog['audit/e.ts'] = prog['audit/e.ts'].replace('import { f }', 'import { f, rare, LIMIT }') + '\nconsole.log(rare(LIMIT))'
  const ran = new Set(['web/src/a.ts#f', 'web/src/a.ts#g', 'web/src/a.ts#rare', 'web/src/b.ts#x'])
  const at = (over: Record<string, string> = {}) => closure('audit/e.ts', ['audit/e.ts'], mem({ ...prog, ...over }), ran).hash
  const g0 = at()
  const ed = (a: string, b: string) => ({ 'web/src/a.ts': prog['web/src/a.ts'].replace(a, b) })
  assert.equal(at(ed(`'wording'`, `'other'`)), g0, 'h was never called')
  assert.notEqual(at(ed('n > 1e9', 'n > 1e8')), g0, 'rare was called')
  assert.notEqual(at(ed('LIMIT = 3', 'LIMIT = 4')), g0, 'a value read by called code')
  assert.notEqual(at(ed('n * 2', 'n * 4')), g0, 'g was called')
}

// The repo's scoring: the scorer is reached, a sentence the page prints is not.
{
  const c = closure('audit/truth/run.ts', ['audit/truth/score-array.ts'], worktree())
  const src = (await import('node:fs')).readFileSync(new URL('../../web/src/pageSays.ts', import.meta.url), 'utf8')
  const parserPath = (await import('node:module')).createRequire(new URL('../../web/package.json', import.meta.url).pathname).resolve('rolldown/parseAst')
  const { parseAst } = await import(parserPath)
  const body = parseAst(src, { lang: 'ts' }).body
  const named = (n: string) => body.findIndex((s: any) => src.slice(s.start, s.end).match(new RegExp(`^export (?:function|const) ${n}\\b`)))
  const reached = new Set(c.reached.get('web/src/pageSays.ts'))
  assert.ok(reached.has(named('changeLogLines')), 'changeLogLines is scoring')
  assert.ok(!reached.has(named('figureNote')), 'figureNote is a surface')
  assert.ok(c.reached.get('web/src/scoreSample.ts')?.length, 'scoreSample is reached')
  assert.ok(!c.reached.has('web/src/syngamyPdf.ts'), 'the PDF is a surface')
}

// The planner's decision.
const k = (all: string, score?: string) => ({ all: { hash: all, files: {} }, ...(score ? { score: { hash: score, files: {} } } : {}) })
assert.equal(decide(k('a', 's'), k('b', 's2')), 'score')
assert.equal(decide(k('a', 's'), k('b', 's')), 'surfaces')
assert.equal(decide(k('a', 's'), k('a', 's')), 'nothing')
assert.equal(decide(k('a'), k('b')), 'run')
assert.equal(decide(k('a'), k('a')), 'nothing')
assert.equal(decide(null, k('a', 's')), 'score')
assert.equal(decide(k('a'), k('b', 's')), 'score')
assert.equal(decide(k('a'), null), null)
console.log('closure and planner checks pass')
