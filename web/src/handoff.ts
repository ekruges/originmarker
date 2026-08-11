/**
 * Hand a reconstructed array straight from Progenitor to Syngamy, in a new tab.
 *
 * The two features are one workflow: Progenitor builds the array of a parent nobody genotyped,
 * and Syngamy is what you point that array at. Doing it by hand means downloading a 16MB file,
 * opening the other page and finding it again in a file picker, which is three steps of
 * clerical work between two halves of the same run.
 *
 * The transfer is `postMessage` between the two tabs, NOT storage. A tab opened with
 * `window.open` keeps a live reference to its opener, so the file can be passed as a structured
 * clone with nothing serialised, nothing persisted and nothing to clean up afterwards. The
 * alternatives are worse in specific ways: `localStorage` caps out around 5MB and this file is
 * three times that; `sessionStorage` is per tab and its inheritance by an opened tab is not
 * something to rely on; IndexedDB works but leaves a copy of an identifiable person's genotype
 * sitting in the browser after the run, which is the one thing this tool has always refused to
 * do.
 *
 * The new tab asks; the opener answers. That ordering matters because the opener cannot know
 * when the new tab's listener is mounted, and a message sent before then is simply lost.
 */

/** Marks a Syngamy tab as one that should expect a file rather than wait for a drop. */
export const HANDOFF_FLAG = 'handoff'

interface ReadyMessage { type: 'om-handoff-ready' }
interface FileMessage { type: 'om-handoff-file'; file: File; role: 'donor'; note: string }

const READY: ReadyMessage['type'] = 'om-handoff-ready'
const FILE: FileMessage['type'] = 'om-handoff-file'

/** True when this page was opened expecting a handoff. Read from the hash, since the app is
 *  hash-routed and a query string before the hash would be dropped by the router. */
export const wantsHandoff = (hash: string): boolean =>
  new URLSearchParams(hash.split('?')[1] ?? '').get(HANDOFF_FLAG) === '1'

/**
 * Opener side. Opens Syngamy in a new tab and answers its request with the file.
 *
 * Returns false when the tab could not be opened, which in practice means a popup blocker: the
 * call has to happen inside the click that caused it or the browser treats it as unsolicited.
 */
export function sendToSyngamy(file: File, note: string): boolean {
  // The blank tab is opened FIRST, with no URL, and only then pointed at Syngamy. Opening it
  // with the URL directly is what a naive version does and it has a bad failure mode: the
  // target differs from the current page only by its hash, so a browser that blocks the popup
  // can fall back to navigating THIS tab to it. That looks like it worked, and it destroys the
  // Progenitor run whose output was being handed over. Measured: with popups blocked,
  // window.open returned null AND the current tab navigated anyway. An empty target has no URL
  // to fall back to, so a blocked open stays blocked.
  const win = window.open('', '_blank')
  if (!win || win === window) return false
  win.location.href = `${window.location.pathname}#/syngamy?${HANDOFF_FLAG}=1`
  const onReady = (e: MessageEvent): void => {
    if (e.origin !== window.location.origin) return
    if ((e.data as ReadyMessage)?.type !== READY) return
    if (e.source !== win) return
    const msg: FileMessage = { type: FILE, file, role: 'donor', note }
    win.postMessage(msg, window.location.origin)
    window.removeEventListener('message', onReady)
  }
  window.addEventListener('message', onReady)
  // Give up rather than leak a listener onto a tab the user closed before it loaded.
  setTimeout(() => window.removeEventListener('message', onReady), 60_000)
  return true
}

/**
 * Receiver side. Announces readiness to the opener and calls back once with the file.
 *
 * Returns a teardown function. Does nothing at all when the page was not opened for a handoff,
 * so the ordinary Syngamy page never listens for messages it will not get.
 */
export function receiveHandoff(
  onFile: (file: File, note: string) => void,
): () => void {
  if (!wantsHandoff(window.location.hash) || !window.opener) return () => {}
  let done = false
  const onMessage = (e: MessageEvent): void => {
    if (done || e.origin !== window.location.origin) return
    const d = e.data as FileMessage
    if (d?.type !== FILE || !(d.file instanceof File)) return
    done = true
    window.removeEventListener('message', onMessage)
    onFile(d.file, d.note)
  }
  window.addEventListener('message', onMessage)
  const ready: ReadyMessage = { type: READY }
  window.opener.postMessage(ready, window.location.origin)
  return () => window.removeEventListener('message', onMessage)
}
