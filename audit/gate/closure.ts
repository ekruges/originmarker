/**
 * What a script's output can depend on, read from the code: the top-level statements it can reach,
 * and the data files those statements name.
 *
 * A program is an entry script and every module it imports, statically or by `import()`, followed to
 * the end. Its roots are every statement of the root files, and every top-level statement of every
 * other module of the program that runs code when the module loads (anything but an import, an export
 * list, a type, a function or class declaration, or a variable whose value is a function, a literal,
 * or a literal object or array). From the roots, each identifier a statement names is followed to
 * the declarations it can mean: a declaration of that name in the same module, or, through an
 * import, the export of that name in the module it comes from. `N.x` on a namespace (`import * as N`
 * or `const N = await import(...)`) follows the export `x`; any other use of the namespace takes the
 * whole module. Names are matched without scope, so a local that shadows a top-level name adds that
 * declaration too: the closure can only be larger than the truth.
 *
 * Module state: where a reached declaration is a mutable binding (`let`, `var`, a class, or a `const`
 * holding an object, array, Map, Set or regular expression), every statement of the program that
 * writes it (scan below) is reached too.
 *
 * Not seen: a mutable object handed to a function, or handed out by one, and written there.
 *
 * Data: a string in a reached statement that is the path of a file under web/public, as the tool
 * names its shipped data, adds that file's bytes. A reached statement that takes a file, directory
 * or path from the environment (`process.env.*FILE`, `*DIR`, `*PATH`) adds every file under
 * web/public, since the job line may name any of them. A package import in a module with a reached
 * statement adds web/package-lock.json.
 *
 * The hash is over each reached statement's syntax tree, in file order, with its file's path, leaving
 * out positions and types: comments, layout, type annotations and type-only declarations, which type
 * stripping erases, do not enter it.
 */
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { createHash } from 'node:crypto'
import { createRequire } from 'node:module'
import { execFileSync } from 'node:child_process'
import { dirname, join, normalize, relative, resolve } from 'node:path'
import { pathToFileURL } from 'node:url'

export const ROOT = resolve(dirname(new URL(import.meta.url).pathname), '../..')

const parserPath = createRequire(join(ROOT, 'web/package.json')).resolve('rolldown/parseAst')
const { parseAst } = await import(pathToFileURL(parserPath).href)

/** Files of a tree, by path relative to the repo root. */
export interface Tree { read(rel: string): Buffer | null; list(dir: string): string[] }

export function worktree(root = ROOT): Tree {
  const walk = (d: string): string[] => readdirSync(join(root, d)).flatMap((n) => {
    const p = join(d, n)
    return statSync(join(root, p)).isDirectory() ? walk(p) : [p]
  })
  return {
    read: (rel) => (existsSync(join(root, rel)) && statSync(join(root, rel)).isFile() ? readFileSync(join(root, rel)) : null),
    list: (dir) => (existsSync(join(root, dir)) ? walk(dir) : []),
  }
}

/** The tree of a commit, read from the object store. */
export function commitTree(commit: string, repo = ROOT): Tree {
  const blobs = new Map<string, string>()
  for (const l of execFileSync('git', ['-C', repo, 'ls-tree', '-r', '-z', commit], { maxBuffer: 1 << 28 }).toString().split('\0')) {
    const m = l.match(/^\d+ blob ([0-9a-f]+)\t(.*)$/s)
    if (m) blobs.set(m[2], m[1])
  }
  return {
    read: (rel) => { const b = blobs.get(rel); return b ? execFileSync('git', ['-C', repo, 'cat-file', 'blob', b], { maxBuffer: 1 << 28 }) : null },
    list: (dir) => [...blobs.keys()].filter((p) => p.startsWith(`${dir}/`)),
  }
}

type Node = { type: string; start: number; end: number; [k: string]: any }
interface Stmt {
  i: number; canon: string; names: string[]; refs: Set<string>; writes: Set<string>; strings: string[]
  effect: boolean; mutable: boolean; typeOnly: boolean; exported: boolean; envPath: boolean; whole: string[]
}
interface Mod {
  file: string; stmts: Stmt[]; byName: Map<string, number[]>; pkg: boolean
  imports: Map<string, { from: string; name: string }>
  exports: Map<string, string>; reexports: Map<string, { from: string; name: string }>; stars: string[]
  deps: string[]
}

const PURE = new Set(['ArrowFunctionExpression', 'FunctionExpression', 'Literal', 'TemplateLiteral'])
function pure(n: Node | null | undefined): boolean {
  if (!n) return true
  if (n.type === 'TSAsExpression' || n.type === 'TSSatisfiesExpression' || n.type === 'TSNonNullExpression') return pure(n.expression)
  if (n.type === 'Literal') return !n.regex
  if (n.type === 'TemplateLiteral') return n.expressions.length === 0
  if (n.type === 'UnaryExpression') return pure(n.argument)
  if (n.type === 'ObjectExpression') return n.properties.every((p: Node) => p.type === 'Property' && !p.computed && pure(p.value))
  if (n.type === 'ArrayExpression') return n.elements.every((e: Node) => e && e.type !== 'SpreadElement' && pure(e))
  if (n.type === 'BinaryExpression' || n.type === 'LogicalExpression') return pure(n.left) && pure(n.right)
  if (n.type === 'ConditionalExpression') return pure(n.test) && pure(n.consequent) && pure(n.alternate)
  if (n.type === 'NewExpression') return n.callee.type === 'Identifier' && ['Map', 'Set', 'WeakMap', 'Array'].includes(n.callee.name)
    && n.arguments.length <= 1 && n.arguments.every((x: Node) => x.type === 'ArrayExpression' && pure(x))
  return PURE.has(n.type)
}
/** A value whose object can be written after it is made. */
function holdsObject(n: Node | null | undefined): boolean {
  if (!n) return false
  if (n.type === 'TSAsExpression' || n.type === 'TSSatisfiesExpression') return holdsObject(n.expression)
  if (n.type === 'ArrowFunctionExpression' || n.type === 'FunctionExpression' || n.type === 'TemplateLiteral' || n.type === 'UnaryExpression') return false
  if (n.type === 'Literal') return !!n.regex
  return true
}

function patternNames(p: Node | null, out: string[]) {
  if (!p) return
  if (p.type === 'Identifier') out.push(p.name)
  else if (p.type === 'ObjectPattern') for (const q of p.properties) patternNames(q.type === 'RestElement' ? q.argument : q.value, out)
  else if (p.type === 'ArrayPattern') for (const q of p.elements) patternNames(q, out)
  else if (p.type === 'RestElement') patternNames(p.argument, out)
  else if (p.type === 'AssignmentPattern') patternNames(p.left, out)
}

/** Keys that hold only types, which type stripping erases. */
const TYPE_KEYS = new Set(['typeAnnotation', 'returnType', 'typeParameters', 'typeArguments', 'superTypeArguments', 'implements'])
const TYPE_ONLY = new Set(['TSInterfaceDeclaration', 'TSTypeAliasDeclaration', 'TSDeclareFunction'])

/** A statement as the code it runs: its syntax tree with positions and types left out, so comments,
 * layout and type annotations do not change it. */
function canonical(n: Node): string {
  return JSON.stringify(n, (k, v) => (k === 'start' || k === 'end' || TYPE_KEYS.has(k) ? undefined
    : typeof v === 'bigint' ? `${v}n` : v))
}

/** Methods that change the object they are called on. */
const MUTATORS = new Set(['set', 'add', 'delete', 'clear', 'push', 'pop', 'shift', 'unshift', 'splice', 'sort', 'reverse',
  'fill', 'copyWithin', 'exec', 'test'])
const WRITES_FIRST = new Set(['assign', 'defineProperty', 'defineProperties', 'setPrototypeOf'])
/** The variable an assignment target, or the object of a method call, belongs to. */
function rootOf(e: Node | null | undefined): string | null {
  while (e && (e.type === 'MemberExpression' || e.type === 'TSNonNullExpression' || e.type === 'TSAsExpression'
    || e.type === 'ParenthesizedExpression' || e.type === 'ChainExpression')) e = e.type === 'MemberExpression' ? e.object : e.expression
  return e?.type === 'Identifier' ? e.name : null
}

/**
 * Identifiers a node names, `N.x` kept as one ref, strings it holds, the `import()`s in it, and the
 * variables it writes: assigned, updated, deleted from, changed by a method (MUTATORS) or by
 * Object.assign and its kin. A variable handed to a function that writes it is not seen.
 */
function scan(n: Node, refs: Set<string>, strings: string[], dyn: Node[], env: { path: boolean } = { path: false },
  writes: Set<string> = new Set()) {
  const write = (e: Node | null | undefined) => { const r = rootOf(e); if (r) writes.add(r) }
  const visit = (x: any, parent?: Node, key?: string) => {
    if (!x || typeof x !== 'object') return
    if (Array.isArray(x)) { for (const y of x) visit(y, parent, key); return }
    if (typeof x.type !== 'string') return
    const t = x.type
    if (t === 'AssignmentExpression') write(x.left)
    if (t === 'UpdateExpression' || (t === 'UnaryExpression' && x.operator === 'delete')) write(x.argument)
    if (t === 'CallExpression' && x.callee.type === 'MemberExpression' && !x.callee.computed) {
      const m = x.callee.property.name
      if (MUTATORS.has(m)) write(x.callee.object)
      if (WRITES_FIRST.has(m) && rootOf(x.callee.object) === 'Object') write(x.arguments[0])
    }
    if (t === 'Identifier') {
      const skip = parent && !parent.computed && (
        (parent.type === 'MemberExpression' && key === 'property')
        || ((parent.type === 'Property' || parent.type === 'MethodDefinition' || parent.type === 'PropertyDefinition'
          || parent.type === 'TSPropertySignature' || parent.type === 'TSMethodSignature' || parent.type === 'TSEnumMember') && key === 'key'
          && !(parent.type === 'Property' && parent.shorthand))
        || (parent.type === 'TSQualifiedName' && key === 'right')
        || ((parent.type === 'LabeledStatement' || parent.type === 'BreakStatement' || parent.type === 'ContinueStatement') && key === 'label'))
      if (!skip) refs.add(x.name)
      return
    }
    if (t === 'MemberExpression' && !x.computed && x.property.type === 'Identifier' && x.object.type === 'MemberExpression'
      && x.object.object.type === 'Identifier' && x.object.object.name === 'process' && x.object.property.name === 'env'
      && /FILE|DIR|PATH/.test(x.property.name)) env.path = true
    if (t === 'MemberExpression' && !x.computed && x.object.type === 'Identifier' && x.property.type === 'Identifier') {
      refs.add(`${x.object.name}.${x.property.name}`)
      return
    }
    if (t === 'Literal' && typeof x.value === 'string') strings.push(x.value)
    if (t === 'TemplateElement' && x.value) strings.push(x.value.cooked ?? x.value.raw)
    if (t === 'ImportExpression') dyn.push(x)
    for (const [k, v] of Object.entries(x)) if (k !== 'type' && k !== 'start' && k !== 'end' && !TYPE_KEYS.has(k)) visit(v, x, k)
  }
  visit(n)
}

/** The `import()` a declarator's value is, awaited or not. */
const importInit = (d: Node): Node | null => {
  const v = d.init?.type === 'AwaitExpression' ? d.init.argument : d.init
  return v?.type === 'ImportExpression' ? v : null
}

function parseModule(file: string, src: string): Mod {
  const ast = parseAst(src, { lang: 'ts' }, file)
  const urlDirs = new Map<string, string>()
  const resolveSpec = (spec: string): string | 'pkg' | null => {
    if (spec.startsWith('node:')) return null
    if (!spec.startsWith('.') && !spec.startsWith('/')) return 'pkg'
    return normalize(join(dirname(file), spec))
  }
  const dynTarget = (x: Node): string => {
    const s = x.source
    if (s.type === 'Literal' && typeof s.value === 'string') { const r = resolveSpec(s.value); if (r && r !== 'pkg') return r }
    if (s.type === 'TemplateLiteral' && s.expressions.length === 1 && s.quasis[0].value.cooked === ''
      && s.expressions[0].type === 'Identifier' && urlDirs.has(s.expressions[0].name)) {
      return normalize(join(urlDirs.get(s.expressions[0].name)!, s.quasis[1].value.cooked))
    }
    throw new Error(`${file}: an import() the closure cannot resolve: ${src.slice(s.start, s.end)}`)
  }
  const m: Mod = { file, stmts: [], byName: new Map(), pkg: false, imports: new Map(), exports: new Map(), reexports: new Map(), stars: [], deps: [] }
  const dep = (spec: string) => { const r = resolveSpec(spec); if (r === 'pkg') m.pkg = true; else if (r) { m.deps.push(r); return r } return null }
  ast.body.forEach((s0: Node, i: number) => {
    let s = s0
    const st: Stmt = { i, canon: canonical(s), names: [], refs: new Set(), writes: new Set(), strings: [], effect: false, mutable: false,
      typeOnly: false, exported: s.type.startsWith('Export'), envPath: false, whole: [] }
    m.stmts.push(st)
    if (s.type === 'ImportDeclaration') {
      const from = dep(s.source.value)
      for (const sp of s.specifiers) {
        if (!from) continue
        const name = sp.type === 'ImportNamespaceSpecifier' ? '*' : sp.type === 'ImportDefaultSpecifier' ? 'default'
          : (sp.imported.name ?? sp.imported.value)
        m.imports.set(sp.local.name, { from, name })
      }
      return
    }
    if (s.type === 'ExportAllDeclaration') {
      const from = dep(s.source.value)
      if (from) { if (s.exported) m.reexports.set(s.exported.name ?? s.exported.value, { from, name: '*' }); else m.stars.push(from) }
      return
    }
    if (s.type === 'ExportNamedDeclaration' && !s.declaration) {
      const from = s.source ? dep(s.source.value) : null
      for (const sp of s.specifiers) {
        const local = sp.local.name ?? sp.local.value, exported = sp.exported.name ?? sp.exported.value
        if (s.source) { if (from) m.reexports.set(exported, { from, name: local }) } else m.exports.set(exported, local)
      }
      return
    }
    if (s.type === 'ExportDefaultDeclaration') {
      const d = s.declaration
      if (d.id?.name) { st.names.push(d.id.name); m.exports.set('default', d.id.name) } else st.names.push('*default')
      if (!d.id) m.exports.set('default', '*default')
      s = d
    } else if (s.type === 'ExportNamedDeclaration') {
      s = s.declaration
    }
    if (TYPE_ONLY.has(s.type) || s.declare) {
      st.typeOnly = true
      if (s.id?.name) st.names.push(s.id.name)
      else if (s.type === 'VariableDeclaration') for (const d of s.declarations) patternNames(d.id, st.names)
    } else if (s.type === 'VariableDeclaration') {
      for (const d of s.declarations) {
        patternNames(d.id, st.names)
        const init = d.init?.type === 'AwaitExpression' ? d.init.argument : d.init
        if (init?.type === 'NewExpression' && init.callee.type === 'Identifier' && init.callee.name === 'URL'
          && init.arguments[0]?.type === 'Literal' && d.id.type === 'Identifier' && init.arguments[1]?.type === 'MemberExpression') {
          urlDirs.set(d.id.name, normalize(join(dirname(file), init.arguments[0].value)))
        }
        if (init?.type === 'MemberExpression' && init.object.type === 'NewExpression' && init.object.callee.name === 'URL'
          && init.object.arguments[0]?.type === 'Literal' && d.id.type === 'Identifier') {
          urlDirs.set(d.id.name, normalize(join(dirname(file), init.object.arguments[0].value)))
        }
        if (init?.type === 'ImportExpression' && d.id.type === 'Identifier' && s.kind === 'const') {
          const to = dynTarget(init); m.deps.push(to); m.imports.set(d.id.name, { from: to, name: '*' })
          scan(init.source, st.refs, st.strings, [])
          st.effect = true
          continue
        }
        if (!pure(d.init)) st.effect = true
        if (s.kind !== 'const' || holdsObject(d.init)) st.mutable = true
      }
    } else if (s.type === 'FunctionDeclaration' || s.type === 'TSEnumDeclaration' || s.type === 'TSModuleDeclaration') {
      if (s.id?.name) st.names.push(s.id.name)
    } else if (s.type === 'ClassDeclaration') {
      if (s.id?.name) st.names.push(s.id.name)
      st.mutable = true
    } else if (st.names.length === 0) {
      st.effect = true
    }
    const dyn: Node[] = []
    if (st.typeOnly) {
      // nothing it names runs
    } else if (!(s.type === 'VariableDeclaration' && s.kind === 'const' && s.declarations.some((d: Node) => importInit(d) && d.id.type === 'Identifier'))) {
      const env = { path: false }
      scan(s, st.refs, st.strings, dyn, env, st.writes)
      st.envPath = env.path
    } else {
      const env = { path: false }
      for (const d of s.declarations) if (!(importInit(d) && d.id.type === 'Identifier')) scan(d, st.refs, st.strings, dyn, env, st.writes)
      st.envPath = env.path
    }
    for (const x of dyn) { const to = dynTarget(x); m.deps.push(to); st.whole.push(to) }
  })
  m.stmts.forEach((st) => { for (const n of st.names) { const a = m.byName.get(n) ?? []; a.push(st.i); m.byName.set(n, a) } })
  return m
}

export interface Closure { hash: string; files: Record<string, string>; reached: Map<string, number[]> }

const sha = (b: Buffer | string) => createHash('sha256').update(b).digest('hex')

/**
 * The closure of `roots` (whole files) in the program of `entry`. `files` maps every file of the
 * program, and every data file reached, to the sha256 of its bytes.
 */
export function closure(entry: string, roots: string[], tree: Tree = worktree()): Closure {
  const mods = new Map<string, Mod>()
  const load = (f: string) => {
    if (mods.has(f)) return
    const b = tree.read(f)
    if (!b) throw new Error(`${f}: imported but not found`)
    const m = parseModule(f, b.toString('utf8'))
    mods.set(f, m)
    for (const d of m.deps) load(d)
  }
  load(entry)
  for (const r of roots) load(r)

  const exportOf = (file: string, name: string, seen = new Set<string>()): [string, number][] => {
    const k = `${file}|${name}`
    if (seen.has(k)) return []
    seen.add(k)
    const m = mods.get(file)!
    if (name === '*') return m.stmts.filter((s) => s.names.length).map((s) => [file, s.i])
    const local = m.exports.get(name)
    if (local) return (m.byName.get(local) ?? []).map((i) => [file, i])
    const own = m.byName.get(name)
    const re = m.reexports.get(name)
    if (re) return exportOf(re.from, re.name, seen)
    if (own && m.stmts[own[0]].exported) return own.map((i) => [file, i])
    return m.stars.flatMap((s) => exportOf(s, name, seen))
  }
  const targets = (m: Mod, ref: string): [string, number][] => {
    const dot = ref.indexOf('.')
    const head = dot < 0 ? ref : ref.slice(0, dot)
    const out: [string, number][] = (m.byName.get(head) ?? []).map((i) => [m.file, i])
    const imp = m.imports.get(head)
    if (imp) {
      if (imp.name === '*' && dot >= 0) out.push(...exportOf(imp.from, ref.slice(dot + 1)))
      else out.push(...exportOf(imp.from, imp.name))
    }
    return out
  }

  const reached = new Map<string, Set<number>>()
  const queue: [string, number][] = []
  const reach = (f: string, i: number) => {
    let s = reached.get(f); if (!s) reached.set(f, s = new Set())
    if (!s.has(i)) { s.add(i); queue.push([f, i]) }
  }
  for (const r of roots) for (const s of mods.get(r)!.stmts) reach(r, s.i)
  for (const m of mods.values()) if (m.file !== entry) for (const s of m.stmts) if (s.effect) reach(m.file, s.i)

  // Who writes each statement's variables, for module state.
  const namedBy = new Map<string, [string, number][]>()
  for (const m of mods.values()) {
    for (const s of m.stmts) {
      for (const ref of s.writes) {
        for (const [f, i] of targets(m, ref)) {
          const k = `${f}|${i}`
          const a = namedBy.get(k) ?? []; a.push([m.file, s.i]); namedBy.set(k, a)
        }
      }
    }
  }
  while (queue.length) {
    const [f, i] = queue.pop()!
    const m = mods.get(f)!, s = m.stmts[i]
    for (const ref of s.refs) for (const [g, j] of targets(m, ref)) reach(g, j)
    for (const w of s.whole) for (const t of mods.get(w)!.stmts) reach(w, t.i)
    if (s.mutable) for (const [g, j] of namedBy.get(`${f}|${i}`) ?? []) reach(g, j)
  }

  const h = createHash('sha256')
  const files: Record<string, string> = {}
  const data = new Set<string>()
  let pkg = false
  const publicFiles = new Set(tree.list('web/public'))
  for (const f of [...mods.keys()].sort()) {
    files[f] = sha(tree.read(f)!)
    const m = mods.get(f)!
    const idx = [...(reached.get(f) ?? [])].sort((a, b) => a - b)
    if (idx.length && m.pkg) pkg = true
    for (const i of idx) {
      if (m.stmts[i].envPath) for (const p of publicFiles) data.add(p)
      if (m.stmts[i].typeOnly) continue
      h.update(`${f}\n${m.stmts[i].canon}\n`)
      for (const str of m.stmts[i].strings) {
        const p = normalize(join('web/public', str))
        if (p.startsWith('web/public/') && publicFiles.has(p)) data.add(p)
      }
    }
  }
  if (pkg) data.add('web/package-lock.json')
  for (const d of [...data].sort()) {
    const b = tree.read(d)
    if (!b) continue
    files[d] = sha(b)
    h.update(`data ${d} ${files[d]}\n`)
  }
  return { hash: h.digest('hex'), files, reached: new Map([...reached].map(([f, s]) => [f, [...s].sort((a, b) => a - b)])) }
}

export const relativeToRoot = (p: string) => relative(ROOT, resolve(p))
