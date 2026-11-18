// Self-check for the file reader both pages use. Run: node src/fileLines.check.ts
// A gzipped table was refused as having no columns, because neither page inflated what was dropped
// on it; the examples were inflated by the code that fetches them and nothing else was.
import assert from 'node:assert/strict'
import { eachLine, isGzip } from './fileLines.ts'

const TABLE = 'probeset_id\tchromosome\tposition\tcall\n'
  + 'AX-1\t1\t1000\tAA\n'
  + 'AX-2\t1\t2000\tAB\n'
  + 'AX-3\t2\t3000\tBB'

const gzip = async (s: string): Promise<Blob> => new Response(
  new Blob([s]).stream().pipeThrough(new CompressionStream('gzip')),
).blob()

const read = async (file: File): Promise<{ lines: string[]; bytes: number }> => {
  const lines: string[] = []
  let bytes = 0
  await eachLine(file, (l) => lines.push(l), (b) => { bytes = b })
  return { lines, bytes }
}

// The magic number, which is what decides: a file saved without its .gz suffix is still gzip, and
// a host serving one with Content-Encoding: gzip hands over bytes that are already inflated.
assert.ok(isGzip(new Uint8Array(await (await gzip('x')).slice(0, 2).arrayBuffer())))
assert.ok(!isGzip(new TextEncoder().encode('probeset_id')))

{
  const plain = await read(new File([TABLE], 'sample.probes'))
  assert.deepEqual(plain.lines, TABLE.split('\n'))
  assert.equal(plain.bytes, TABLE.length, 'a plain file reports its own length')

  const blob = await gzip(TABLE)
  const gz = await read(new File([blob], 'sample.probes.gz'))
  assert.deepEqual(gz.lines, plain.lines, 'a gzipped table reads as the same table')
  // The name decides nothing: the same bytes under a plain name read the same way.
  assert.deepEqual((await read(new File([blob], 'sample.probes'))).lines, plain.lines)
  // Progress is counted in the bytes of the file on disk, because the bar is drawn against
  // file.size. Counted after inflation it would run past the end of the bar.
  assert.ok(gz.bytes <= blob.size, `${gz.bytes} of ${blob.size} bytes read`)
  assert.ok(blob.size < TABLE.length * 4, 'the compressed file is not the inflated one')
}

// A file whose last line has no newline keeps that line; a trailing newline adds no empty one that
// the header parser would then see.
{
  const end = await read(new File(['a\nb\n'], 'x.csv'))
  assert.deepEqual(end.lines, ['a', 'b'])
}

console.log('fileLines.check.ts: a dropped table reads the same gzipped or not, by its first two '
  + 'bytes, and progress stays in the bytes of the file on disk')
