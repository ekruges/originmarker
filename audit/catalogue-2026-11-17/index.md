# Every result OriginMarker can print

One full screenshot of each, taken in the app at version 6.14.0, on public arrays only. Every
picture is a whole result card, top to bottom, as the page draws it. Where a case only appears
inside an opened change, that change is open in the picture.

Cells of the first study are named by the number they carry in its published table, cell 1 to
cell 31. Everything else is described: a sperm donor, an egg donor, a pronucleus, a cell of the
second study. "Published" is what that study's own table states for that cell.

Cases on the spine with no picture are in `MISSING.md`, each with its reason. That file also lists
seven things this catalogue found that look wrong.

Words the page prints, in plain terms. A trophectoderm biopsy is a few outer cells of a day-six
embryo. A blastomere is one cell of an embryo a few days old. A pronucleus is one parent's half of
a just-fertilised egg, before the two halves join. Androgenetic means the cell carries the father's
genome and no copy from the mother; gynogenetic is the mirror. The Mendelian test asks whether an
allele that one parent could not have given is there. The karyotype is the published count and
shape of the chromosomes for that cell.

## The spine

**The headline, 01 to 23.** Both parents with and without changes; the father and a second parent,
measured and rebuilt; the mother and a second parent; one genome of the father's, one of the
mother's read two ways; neither parent; both parents named on a one-genome cell; the seven reasons
the page says UNCLEAR; and the material word, one per word, read with no parent.

**The changes, 24 to 40.** A whole chromosome missing a parent's copy; a copy lost with no parent
named; a stretch found two ways; a stretch from the absence scan; a stretch named by the Mendelian
test with one parent; a stretch from the allele scan; a stretch with no parent named; a stretch
read by two scans at once; a second copy lost inside a lost chromosome; a copy lost with no parent
loaded; an extra copy with its parent named, with one parent, and with none; parents unbalanced;
the seen-once group; a class from the taxonomy; a change inheriting the genome's own parent.

**The rest, 41 to 55.** Finding the father among several files; the eight sex-chromosome readings;
a file from a chip the tool is not proven on; a widely disturbed genome; the germline page's own
result; a change that draws its figure and one that cannot be drawn; a sex reading that was not
made, with its reason.

## The headline

**01 Both parents, nothing changed.** Cell 4, a trophectoderm biopsy, with both its parents
loaded. Published: all 22 chromosomes normal. The card prints BOTH PARENTS and "Both parents'
genomes are in this cell. No change found." `01-both-parents-no-change.png`

**02 Both parents, with changes.** Cell 22, a cleavage-stage cell, with both its parents.
Published: ten changes, each with its parent. The card prints BOTH PARENTS, ten changes, and one
more stretch seen once. `02-both-parents-changes.png`

**03 The father and a second parent.** Cell 4 again, read against the sperm donor alone.
Published: all 22 normal. The card prints FATHER + SECOND PARENT and "The father's genome is in
this cell, with a second parent's copy. No change found." `03-father-and-second-parent.png`

**04 The same, with the father rebuilt from his own haploid cells.** Cell 4 read against a father
the tool reconstructed from his pronuclei rather than measured. The card prints the same words and
carries the badge "father rebuilt". `04-father-rebuilt.png`

**05 The mother and a second parent.** Cell 4 read against the egg donor alone. The card prints
MOTHER + SECOND PARENT and "The mother's genome is in this cell, with a second parent's copy. No
change found." `05-mother-and-second-parent.png`

**06 One genome, the father's.** A paternal pronucleus with both donors loaded. Published: all 22
normal. The card prints ANDROGENETIC and "Paternal genome only. No maternal complement."
`06-androgenetic.png`

**07 One genome, the mother's, read against the father.** A maternal pronucleus with the sperm
donor loaded. The card prints GYNOGENETIC in the undecided colour, with "No paternal contribution.
Maternal in origin. That holds if the loaded array is this sample's parent, which the alleles
cannot show." `07-gynogenetic.png`

**07b The same cell read against an egg donor instead.** The same maternal pronucleus, this time
with only an egg donor loaded, who is not its mother. The card prints ANDROGENETIC, the opposite
word, with the same caveat: "Paternal genome only. No maternal complement. That holds if the loaded
array is this sample's parent, which the alleles cannot show." The word follows whose genome is
absent, so it follows which adult was loaded. `07b-gynogenetic-mother-loaded.png`

**08 Neither parent.** A maternal pronucleus with both donors loaded, neither of whom is its
parent. The card prints NEITHER PARENT and "Neither declared parent accounts for this genome."
`08-neither-parent.png`

**09 Both parents named on a one-genome cell.** A paternal pronucleus with both donors loaded. The
card prints BIPARENTAL, "Both parents contributed. No change found; 16 more stretches seen once,
not confirmed", beside its own badge reading "uniparental homozygous", which says the cell holds
one genome. The two disagree; see fault C in `MISSING.md`. `09-biparental-word.png`

**10 UNCLEAR: the array failed.** Cell 21, with both parents. Published: no call on any of the 22
chromosomes. The card prints "The array failed: its calls are noise, so no parent's genome can be
measured in it." `10-unclear-array-failed.png`

**11 UNCLEAR: a parent's genome is not in this cell.** Cell 5, read against the sperm donor and an
egg donor who is not its mother. Published: all 22 normal. The card prints "The mother's genome is
not in this cell, so the loaded array may not be hers."
`11-unclear-genome-not-in-cell.png`

**12 UNCLEAR: too few chromosomes could be read.** No picture; see `MISSING.md`.

**13 UNCLEAR: the genome is not shown.** Cell 2, with both parents. Published: chromosome 6 lost
whole, the father's. The card prints "The mother's genome is not shown in this cell: found on 1 of
21 chromosomes read. 2 changes found; 7 more stretches seen once, not confirmed."
`13-unclear-genome-not-shown.png`

**14 UNCLEAR: two readings disagree.** Cell 10, with both parents. Published: chromosomes 5, 6,
10, 12, 13 and 19 lost whole, all the father's. The card prints "Two readings disagree. The
father's genome is found on 0 of 22 chromosomes read, but his alleles are missing across the
genome at 3.3x the noise ceiling. 6 changes found."
`14-unclear-two-readings-disagree.png`

**15 UNCLEAR: the copy number could not be read.** Cell 26, with both parents. Published: nine
chromosomes with no copy at all and thirteen with one. The card prints "The copy number could not
be read, so one genome cannot be told from two." `15-unclear-copy-number-not-read.png`

**16 UNCLEAR: which parent's genome this is.** A paternal pronucleus read against the sperm donor
alone. The badge reads "uniparental homozygous", so the cell holds one genome, and the card prints
UNCLEAR with "Not separated by the available evidence. No change found."
`16-unclear-which-parent-alone.png`

**17 No parent: TROPHECTODERM.** Cell 4 with nothing loaded. The card prints the material word and
"Read with no parent: the stage, the sex chromosomes and copy number. No change found."
`17-no-parent-trophectoderm.png`

**18 No parent: BULK.** Cell 3 with nothing loaded. Published: a piece of the long arm of
chromosome 6 lost, the father's. `18-no-parent-bulk.png`

**19 No parent: BLASTOMERE.** Cell 30 with nothing loaded. Published: nine changes, each with its
parent. The card prints eight changes and, above them, the banner "DAMAGED SAMPLE: call rate 58.2%
is below the 70% this tool can reference against." `19-no-parent-blastomere.png`

**20 No parent: SINGLE CELL.** Cell 19 with nothing loaded. Published: two chromosomes lost, one in
three copies, and a gained arm. `20-no-parent-single-cell.png`

**21 No parent: HAPLOID.** A paternal pronucleus with nothing loaded. `21-no-parent-haploid.png`

**22 No parent: FAILED.** Cell 21 with nothing loaded. Published: no call on any of the 22
chromosomes. The word is printed in the undecided colour, with the banner "DAMAGED SAMPLE: this
array did not resolve to any material." `22-no-parent-array-failed.png`

**23 No parent: UNKNOWN.** No picture; see `MISSING.md`.

## The changes

**24 A whole chromosome, a parent's copy missing.** Cell 22 with both parents, chromosome 2 open.
Published for that chromosome: lost whole, the father's. The row reads "chr2 / whole / father's
copy missing / consistent", and the opened call reads "The alleles show one parent's copy missing
along the whole chromosome. 89 of 96 calls of this kind agree with the published karyotype."
`24-whole-chromosome-consistent.png`

**25 A copy lost, no parent named.** Cell 2 with both parents, the 6.9 Mb stretch of chromosome 11
open. The row reads "chr11 / 6.9 Mb / copy lost, parent not named / one reading", and the opened
call reads "One reading: copy number from intensity. 3 of 7 calls of this kind sit where the
published karyotype states a change." `25-copy-lost-parent-not-named.png`

**26 A stretch found two ways.** Cell 3 with both parents, chromosome 6 open. Published: a piece of
the long arm of chromosome 6 lost, the father's. The row reads "chr6 / 117.8 Mb / father's copy
missing / found two ways", and the opened call reads "Two readings agree: the absence scan, and
the Mendelian test over the same stretch. 19 of 26 calls of this kind agree with the published
karyotype." `26-piece-found-two-ways.png`

**27 A stretch from the absence scan alone.** Cell 2 with both parents, the 14.2 Mb stretch of
chromosome 3 open: "father's copy missing / one reading", with "The Mendelian test over the same
stretch did not confirm it." `27-piece-absence-scan-only.png`

**28 A stretch named by the Mendelian test, one parent loaded.** Cell 28 read against the sperm
donor, chromosome 6 open. Published: chromosome 6 part lost and part gained, the father's. The row
reads "chr6 / 64.3 Mb / father's copy missing / consistent", and the opened call reads "The
Mendelian test against the loaded parent names whose copy is missing. 13 of 15 calls of this kind
agree with the published karyotype." `28-piece-mendelian-one-parent.png`

**29 A stretch from the allele scan.** Cell 22 with both parents, the 6.1 Mb stretch of chromosome
4 open: "father's copy missing / one reading", with "One reading: the scan of one parent's alleles
along the chromosome. 3 of 17 calls of this kind agree with the published karyotype."
`29-piece-allele-scan.png`

**30 A stretch with no parent named.** Cell 17 read against the sperm donor, chromosome 6 open.
Published: the long arm of chromosome 6 lost, the father's. The row reads "a parent's copy
missing, parent not named / one reading". The figure under it names the father, which is fault A
in `MISSING.md`. `30-piece-no-parent-named.png`

**31 A stretch read by two scans at once.** Cell 2 with both parents, the 23.6 Mb stretch of
chromosome 15 open: "one reading", with "a copy lost here, reported by the copy scan and by an
allele scan over the same stretch." `31-piece-two-scans-one-stretch.png`

**32 A second copy lost inside a lost chromosome.** Cell 2 with both parents, the 50.4 Mb stretch
of chromosome 6 open, on a chromosome the same card calls lost whole. The row reads "second copy
lost too". New in 6.13.0. `32-second-copy-lost-too.png`

**33 A copy lost with no parent loaded.** A cell of the second study, nothing loaded, chromosome 16
open. Published: a piece of the long arm of chromosome 16 lost. The row reads "copy lost, parent
not named / position only". `33-copy-lost-position-only.png`

**34 An extra copy with its parent named.** Cell 19 with both parents, the 74.3 Mb stretch of
chromosome 6 open. Published for that cell: the long arm of chromosome 6 gained, the father's, and
that is what the card names. The row reads "extra copy, the father's / one reading", and the
opened call reads "One reading: an extra copy, with its parent named. 1 of 4 calls of this kind
agree with the published karyotype." `34-extra-copy-parent-named.png`

**35 An extra copy named with one parent, off the run's other cells.** No picture; see
`MISSING.md`.

**36 An extra copy with no parent named.** The same cell and card, the other stretch of chromosome
6 open: "extra copy, parent not named / one reading", read from copy number alone. The two
stretches sit on one chromosome and only one of them carries a parent.
`36-extra-copy-parent-not-named.png`

**37 Parents unbalanced.** Cell 22 with both parents, chromosome 6 open in the group headed "Also
read, not counted as a change": "parents unbalanced / also seen", with "The parents' alleles read
unbalanced along the chromosome. Not counted as a change." `37-parents-unbalanced.png`

**38 The seen-once group.** Cell 31 with both parents. Published: twelve chromosomes lost whole,
with their parents. The card counts ten changes and lists six more stretches under the heading
"Seen once, not confirmed", which the headline does not count.
`38-seen-once-not-confirmed.png`

**39 A class from the taxonomy.** Cell 3 with both parents, the 104.5 Mb stretch of chromosome 6
open: a copy-neutral stretch from one parent only. The row prints the class as a code word,
"segmental upd", which is fault B in `MISSING.md`. `39-class-from-the-taxonomy.png`

**40 A change inheriting the genome's own parent.** No picture; see `MISSING.md`.

## The rest

**41 Finding the father among several files.** The bundled one-parent example: two men offered as
the father of three cells, the study's sperm donor and an unrelated adult. The panel reads each
candidate against every cell and prints "its genome is in 3 of 3 arrays, and no other candidate's
is in any", then scores the cells against the one it named. It also states the limit it cannot
pass: the check cannot tell a parent from a close relative of that parent.
`41-finding-the-father.png`

**42 Sex: XX.** Cell 22 read against the sperm donor, with the XX badge on the headline.
`42-sex-xx.png`

**43 Sex: XY.** Cell 20 with both parents, the XY badge outlined rather than filled because an
ordinary male karyotype is not an abnormality. Published: chromosome 5 lost and the long arm of
chromosome 6 lost, both the father's; the card counts one of the two. `43-sex-xy.png`

**44 Sex: one X.** Cell 22 read with nothing loaded. With no parent the material reads haploid, one
genome, and one X is then the normal complement rather than a chromosome missing. Published for
that cell: ten changes; with no parent the card prints one. `44-sex-one-x.png`

**45 Sex: XX or one X.** Cell 13 with nothing loaded: chromosome X sits at the level of the
autosomes, and whether the cell holds one genome or two is not settled, so both readings are
printed. Published for that cell: four chromosomes lost whole, chromosome 16 gained with a piece
of it lost, and a piece of 8 lost. With no parent the card finds none of those and prints one
change, a 5.7 Mb copy lost on chromosome 12. `45-sex-xx-or-one-x.png`

**46 Sex: XY read with no parent.** Another cell of the second study with nothing loaded:
chromosome X a copy step under the autosomes with a chromosome Y, which one genome cannot hold.
Published for that cell: chromosome 8 in three copies and chromosome 15 lost. The card finds
neither, and prints one other change, a 3.6 Mb extra copy on chromosome 17 with no parent named.
`46-sex-xy-no-parent.png`

**47 Sex: one X and no Y.** No picture; see `MISSING.md`.

**48 Sex: three X.** Cell 8 read against the sperm donor. The badge reads 47,XXX: chromosome X more
than two copies and no chromosome Y, and the badge is filled rather than outlined because the
reading is itself an abnormality. Published for that cell: pieces of chromosomes 12 and 7 lost,
both the father's, and the card finds both. `48-sex-47xxx.png`

**49 Sex: two X and a Y.** No picture; see `MISSING.md`.

**50 A file from a chip the tool is not proven on.** A public array of a third study, on another
manufacturer's chip. The file is refused as it is read, with every cause named, and no result card
is made. `50-chip-not-proven-refused.png`

**51 A widely disturbed genome.** Cell 26 with nothing loaded. Published: nine chromosomes with no
copy at all and thirteen with one. The card prints eighteen changes and, under the list, "10 of 22
chromosomes carry a change: this cell's genome is widely disturbed."
`51-widely-disturbed-genome.png`

**52 The germline page's own result.** Five haploid cells of one man, the bundled example. The page
prints "One parent, enough products to start a reconstruction", with "Every pair agrees at the rate
two products of one parent do", and the counts beneath: five arrays in, five usable, one parent
group, five in the largest group. `52-germline-father-rebuilt.png`

**53 A change that draws its own figure.** Cell 22 with both parents, the 31.3 Mb stretch of
chromosome 5 open, which is long enough and named, so the page draws it.
`53-change-with-its-figure.png`

**54 A change that cannot be drawn on its own.** The same cell and stretch as 45, opened: the page
says a single cell cannot show this stretch by eye, that it is shorter than one band, and draws the
whole chromosome instead. `54-change-without-its-figure.png`

**55 A sex reading that was not made, with its reason.** An adult woman's own cells from the same
study, read with nothing loaded. No sex badge on the headline, and under "How this was read", which
is open in the picture, two reasons in full: "Not read: chrX sits 0.25 log2 from this sample's
autosomes, beyond this array's own noise but short of the 0.40 log2 copy step, so one X or two is
not separable here", and, for the sex taken from heterozygosity, that two copies of the autosomes
are not shown with no parent loaded. `55-sex-not-read-with-reason.png`
