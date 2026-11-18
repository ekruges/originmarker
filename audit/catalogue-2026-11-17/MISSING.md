# What the catalogue could not show

Every case on the spine that no open public array would produce, with the reason. Nothing here
was faked and no array was built to fill a gap. The numbers are the spine numbers in `index.md`.

## Cases with no picture

**12. UNCLEAR: too few chromosomes could be read.** Not reachable from any array. The sentence is
printed when the per-chromosome parent check comes back "unread", and that value is set in one
place in the code: where a loaded parent has no marker tally at all. That happens when no parent
is loaded, and then the card prints the material word instead, or when the shipped founder panel
fails to load. So the wording belongs to a failed panel, not to any array.

**23. No parent: UNKNOWN.** The material word "unknown" is given to an array that reports neither a
heterozygous rate nor a call rate. Every array that can be read at all reports both. An array that
cannot be read is refused earlier, with its own wording.

**35. An extra copy named with one parent, from the run's other cells.** The one public cell whose
extra copy is named with both parents loaded (cell 19, the long arm of chromosome 6, the father's)
reads "extra copy, parent not named" when the same cell is read against the father alone with five
other cells of the same mother in the run. The feature was measured on built extra copies, not on
real ones, and the closing plan's own open ledger says no public cell gives a clean extra copy with
its parent named.

**40. A change in a one-genome cell, the parent inherited from the genome call.** Scored every open
pronucleus of the public study, with the father alone and with both parents. In each one the allele
scan or the Mendelian test named the parent of every change, and the inherited grade is only used
where both of those decline. So the word exists in the code and no public cell reaches it.

**47. One X and no Y (45,X).** No open public array reads one X with no Y in a cell shown to hold
two copies. Checked the whole no-parent sex table (330 readings) and, with a father loaded, the
second study's own chromosome-X experiment, which is the only public place such a cell would sit.

**49. Two X and a Y (47,XXY).** The release of 17 November removed this reading from every public
array it used to appear on: 13 cells and one adult woman's array, all of which were XY or
unresolved once the rule changed. On the current code no open public array prints it.

## What this catalogue found that looks wrong

Written down rather than fixed, as asked.

**A. One card says two different things about one stretch.** On cell 17 read against the father
alone, the change row reads "a parent's copy missing, parent not named" and the figure directly
under it, in the same opened row, reads "Chromosome 6, 31.3 to 168.8 Mb: the father's copy is
missing." The row withholds the parent and the picture under it names him. See
`30-piece-no-parent-named.png`.

**B. A class from the taxonomy prints its code word.** The change row for a copy-neutral stretch
prints "segmental upd" in the column that otherwise says things like "the father's copy missing".
The code already holds a plain phrase for every class ("segment from one parent only") and the row
does not use it. See `39-class-from-the-taxonomy.png`.

**C. A one-genome cell can be headed BIPARENTAL.** One paternal pronucleus read with both donors
prints "BIPARENTAL. Both parents contributed." beside its own badge reading "uniparental
homozygous", which says the cell holds one genome, and directly above its own notice reading "None
of the loaded mother's genome is in this array". Three statements on one card, two of them against
the headline. The same card is also marked a damaged sample at a 55.7% call rate, so nothing on it
should be read at all. See `09-biparental-word.png`.

**D. With only the mother loaded the header ends in "vs" and nothing.** The line above the headline
reads "<file> · vs " with no name after it, because the name printed there is always the sperm
donor's. See `05-mother-and-second-parent.png`.

**E. A gzipped file is refused by the first page and accepted by the second.** The Syngamy drop zone
does not unzip a dropped file, so a `.probes.gz` is refused as having no columns, although the file
picker offers `.gz` and the Germline page's own text says a table may be "gzipped or not". Both
pages read files with the same code; neither unzips. Only the bundled examples are unzipped, by
the code that fetches them.

**F. The founder panel is read in a way that fails on a nearly full disk, and the run carries on
without it.** The panel is fetched with `res.blob()`, which stages a 17.6 MB body through the
browser's blob store and needs room on the volume holding the browser profile. With 6.9 MB free
that call fails outright, where reading the same file as a buffer in the same page succeeds. The
result is then kept for the whole page, so one failed fetch leaves every file afterwards scored
with no chip check, and, measured against the same arrays read with the panel, a run with one
parent loaded prints a different material, a different sex, a different genome-level class and a
different set of calls with different coordinates, and says none of it. The page logs one line.
This is the same failure the owner has seen before on a full disk, and it is why several runs
behind this catalogue had to be redone.

**G. One cell can be headed with either uniparental word, depending on which adult is loaded.** The
same maternal pronucleus reads GYNOGENETIC against the sperm donor and ANDROGENETIC against an egg
donor who is not its mother. Both cards carry the caveat that the reading holds only if the loaded
array is this sample's parent, so the page is not claiming more than it has; but the headline word
itself flips, and a reader who sees only one card has no way to tell. See `07-gynogenetic.png` and
`07b-gynogenetic-mother-loaded.png`.
