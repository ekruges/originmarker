// No name from the code in printed prose. Run: node --experimental-strip-types src/prose.check.ts
//
// Over every source file the page and the report are built from (web/src, checks left out), every
// string literal with a space in it is prose a reader may see. None may hold a camelCase word
// outside ALLOW, or a word ending in .ts or .tsx. Comments are not strings and are not read.
// A template literal's ${...} parts are code, not prose, and are left out.
import assert from 'node:assert/strict'
import { readdirSync, readFileSync } from 'node:fs'

/**
 * Words that look like code and are not: the chromosome names chrX and chrY; the UCSC tables and the
 * tool the report cites by name (chromInfo, cytoBand, refGene, wgEncodeUwRepliSeq, liftOver); and two
 * terms of the field, rsID (a dbSNP identifier) and gDNA (genomic DNA).
 */
export const ALLOW = new Set(['chrX', 'chrY', 'chromInfo', 'cytoBand', 'refGene', 'wgEncodeUwRepliSeq', 'liftOver', 'rsID', 'gDNA',
  // The panel page's sources and units: the gnomAD and dbSNP databases, the deCODE map, cM, dNTP.
  'gnomAD', 'dbSNP', 'deCODE', 'cM', 'dNTP'])
const CAMEL = /\b[a-z]+[A-Z][A-Za-z0-9]*\b/g
const FILE = /\b[\w-]+\.tsx?\b/g

/** The string literals of one source text, template expressions removed. */
export function literals(src: string): { text: string; line: number }[] {
  const out: { text: string; line: number }[] = []
  let i = 0; let line = 1
  const n = src.length
  while (i < n) {
    const c = src[i]
    if (c === '\n') { line += 1; i += 1; continue }
    if (c === '/' && src[i + 1] === '/') { while (i < n && src[i] !== '\n') i += 1; continue }
    if (c === '/' && src[i + 1] === '*') { const e = src.indexOf('*/', i + 2); for (let k = i; k < (e < 0 ? n : e); k += 1) if (src[k] === '\n') line += 1; i = e < 0 ? n : e + 2; continue }
    if (c === "'" || c === '"' || c === '`') {
      const start = line; let text = ''; i += 1
      while (i < n && src[i] !== c) {
        if (src[i] === '\\') { text += src[i + 1] ?? ''; i += 2; continue }
        if (c === '`' && src[i] === '$' && src[i + 1] === '{') {
          let depth = 1; i += 2
          while (i < n && depth) { if (src[i] === '{') depth += 1; else if (src[i] === '}') depth -= 1; else if (src[i] === '\n') line += 1; i += 1 }
          text += ' '
          continue
        }
        if (src[i] === '\n') { if (c !== '`') break; line += 1 }
        text += src[i]; i += 1
      }
      i += 1
      out.push({ text, line: start })
      continue
    }
    i += 1
  }
  return out
}

/** Code words and file names in the prose literals of one source text. */
export function proseFaults(src: string): string[] {
  const out: string[] = []
  for (const { text, line } of literals(src)) {
    // Prose has a space; markup (attribute="value") and a web address are not prose.
    if (!/\S \S/.test(text) || /\w="/.test(text) || /^https?:\/\//.test(text)) continue
    for (const w of text.match(CAMEL) ?? []) if (!ALLOW.has(w)) out.push(`${line}: ${w}`)
    for (const w of text.match(FILE) ?? []) out.push(`${line}: ${w}`)
  }
  return out
}

{
  assert.deepEqual(proseFaults("const a = 'prefer dropoutFromReplicates here'"), ['1: dropoutFromReplicates'])
  assert.deepEqual(proseFaults('const a = `a sentence ${fooBar} on chrX`'), [], 'template code is not prose, chrX is allowed')
  assert.deepEqual(proseFaults("// a comment naming fooBar.ts\nconst b = 'oneWord'"), [], 'comments and one-word strings are not prose')
  assert.deepEqual(proseFaults("const c = 'read in (timing.ts) here'"), ['1: timing.ts'])
}

const dir = new URL('./', import.meta.url).pathname
/** The files the Syngamy page and its report are built from: Syngamy.tsx and everything it imports. */
const closure = new Set<string>()
const stack = ['Syngamy.tsx']
while (stack.length) {
  const f = stack.pop()!
  const file = [f, `${f}.ts`, `${f}.tsx`].find((x) => /\.tsx?$/.test(x) && readdirSync(dir).includes(x))
  if (!file || closure.has(file)) continue
  closure.add(file)
  for (const m of readFileSync(dir + file, 'utf8').matchAll(/from '\.\/([\w.]+?)'/g)) stack.push(m[1])
}
assert.ok(closure.has('syngamyPdf.ts') && closure.has('scoreSample.ts') && closure.size > 50, `closure: ${closure.size}`)
const faults: string[] = []
for (const f of [...closure].sort()) {
  for (const x of proseFaults(readFileSync(dir + f, 'utf8'))) faults.push(`${f}:${x}`)
}
// The marker panel's page: its prose, and its one-word labels, which a code name reads as.
for (const f of ['SearchPanel.tsx', 'App.tsx', 'PanelTable.tsx', 'DocsPage.tsx']) {
  const src = readFileSync(dir + f, 'utf8')
  for (const x of proseFaults(src)) faults.push(`${f}:${x}`)
  for (const m of src.matchAll(/\blabel="([^"\s]+)"/g)) if (/_|[a-z][A-Z]/.test(m[1])) faults.push(`${f}: label ${m[1]}`)
}
assert.deepEqual(faults, [], `code words or file names in printed prose:\n${faults.join('\n')}`)
console.log('prose.check: ok')
