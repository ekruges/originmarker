/**
 * Reading a dropped file line by line, gzipped or not.
 *
 * ONE READER FOR BOTH PAGES. Syngamy and Germline each carried their own copy of this loop, and
 * neither inflated anything: a `.probes.gz` dropped on either was parsed as its compressed bytes,
 * found no header row and was refused as having no columns, although the Germline page's own text
 * says a table may be gzipped or not and the bundled examples are themselves gzipped. Two copies
 * meant the fix would have had to be made twice, which is how they came to differ in the first
 * place.
 *
 * THE BYTES SAY WHETHER IT IS GZIPPED, NOT THE NAME. A file saved without its `.gz` suffix is still
 * gzip, and a host that serves one with Content-Encoding: gzip hands over bytes that are already
 * inflated. The magic number is the only thing that is true in both cases, and it is the same test
 * `answerKeys.panelText` makes on the founder panel.
 *
 * PROGRESS IS COUNTED IN THE BYTES OF THE FILE ON DISK, not in the inflated ones, because the bar
 * it drives is drawn against `file.size`. The count is taken before the inflater rather than after.
 *
 * Self-check: node --experimental-strip-types src/fileLines.check.ts
 */

/** Gzip's magic number: the first two bytes of any gzip member. */
export const isGzip = (head: Uint8Array): boolean => head[0] === 0x1f && head[1] === 0x8b

/** Every line of a file, inflating it where it is gzipped. `onChunk` gets bytes of the file read. */
export async function eachLine(
  file: File, fn: (line: string) => void, onChunk?: (bytesRead: number) => void,
): Promise<void> {
  const gz = isGzip(new Uint8Array(await file.slice(0, 2).arrayBuffer()))
  let read = 0
  const counted = file.stream().pipeThrough(new TransformStream<Uint8Array<ArrayBuffer>, Uint8Array<ArrayBuffer>>({
    transform(chunk, c) { read += chunk.byteLength; c.enqueue(chunk) },
  }))
  const reader = (gz ? counted.pipeThrough(new DecompressionStream('gzip')) : counted).getReader()
  const decoder = new TextDecoder()
  let carry = ''
  for (;;) {
    const { done, value } = await reader.read()
    if (done) break
    carry += decoder.decode(value, { stream: true })
    let i = carry.indexOf('\n')
    while (i >= 0) {
      fn(carry.slice(0, i))
      carry = carry.slice(i + 1)
      i = carry.indexOf('\n')
    }
    onChunk?.(read)
    // Yield, so the log and the bars paint while a large file streams.
    await new Promise((r) => { setTimeout(r, 0) })
  }
  if (carry) fn(carry)
}
