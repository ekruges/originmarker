# Bundled genetic map (GRCh38)

Recombination map used by `genetic_map.py` to annotate each candidate marker's
approximate genetic distance (cM) from the pathogenic variant, its recombination fraction,
and whether a recombination hotspot sits between them.

The autosomes are sex-averaged. **chrX is not** (see Coverage below), and a sex-averaged
value describes neither parent: linkage runs through one carrier, of known sex. The map
source string carried alongside every cM value says so, per chromosome.

## Provenance

- **Source:** deCODE 2019 genetic map (Halldorsson et al., *Science* 2019), sex-averaged
  on the autosomes,
  lifted to GRCh38 and distributed in PLINK format by the Beagle project.
- **Downloaded from:** `https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip`
- **Variant used:** the `no_chr_in_chrom_field/` set (chromosome field is `11`, not
  `chr11`), matching the bare chromosome names `panelbuilder` gets back from gnomAD.
- **Retrieved:** 2026-07-15.

## Why these files look like this

The upstream PLINK map is four columns (`chrom rsid cM pos`) and ~77 MB uncompressed
across 23 chromosomes. Only `pos` and `cM` are ever read, so each file here is reduced to
those two columns and gzipped:

```
198062 0
198510 0.000421
```

That is `<pos> <cM>`, ascending by position, one chromosome per file
(`chr11.pos_cm.gz`). The whole set is 23 MB, small enough to commit, which is the point:
the Docker image needs no network to build, and the golden test stays deterministic and
offline. Nothing is recomputed or interpolated at rest; interpolation happens in
`genetic_map.py` at query time.

## Regenerating

```bash
curl -sSLO https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip -q plink.GRCh38.map.zip 'no_chr_in_chrom_field/*' -d /tmp/gmaps
cd /tmp/gmaps/no_chr_in_chrom_field
for f in plink.chr*.GRCh38.map; do
  chr=$(echo "$f" | sed 's/plink.chr\(.*\).GRCh38.map/\1/')
  awk '{print $4" "$3}' "$f" | gzip -9 > "chr${chr}.pos_cm.gz"
done
```

## Coverage and the fallback

chrY and chrM have no map and are not bundled.

**chrX is bundled, but it is the female map, not a sex-averaged one.** deCODE reports X
distances in female meioses by convention, and the bundled file agrees: chrX reads
1.297 cM/Mb against 1.175 (chr7) and 1.160 (chr8). A sex-averaged X would have to fall
*well below* the autosomes, since paternal meioses are a third of X transmissions and
contribute no crossovers outside the PARs; this one sits 10% above them, so it was never
scaled. Non-PAR X carries 180.84 cM where a 2/3-scaled map would carry ~120. `genetic_map`
therefore returns a distinct source string for chrX: it is correct for a female carrier,
and for a male carrier it overstates recombination. The self-check asserts the rate
comparison above, so swapping in a scaled X map without relabelling fails loudly.

If a chromosome file is absent, `genetic_map.load()` returns a map with `approx=True` that falls back to a uniform
**1 cM/Mb approximation**. The same fallback and the same label apply to a position that
falls *beyond either end* of a bundled map: the maps stop short of the telomeres (chr11 by
~11 kb), and cM there is extrapolated at 1 cM/Mb rather than clamped to the last reading.
Clamping is the dangerous choice and was the original bug: it collapsed everything past the
edge onto one cM, reporting theta=0, which is not merely unlabelled but *ranks as a perfect
marker*.

Every surface that displays a cM value derived that way
(UI, CSV, XLSX, JSON, PDF) must label it as an approximation. That labelling is a
correctness requirement, not a nicety: a fabricated 1 cM/Mb figure presented as a real map
reading would misinform a recombination-risk judgement.

## Caveat for interpretation

This is a **population-averaged** map. It describes expected recombination across meioses
in a reference cohort, not what happened in a particular meiosis in a particular family.
Like expected heterozygosity (2pq), it is a prior. A marker with a low recombination
fraction is *likely* to stay in phase with the locus; it is not guaranteed to have done so
in the embryo being tested. This is exactly why the Layer-B protocol requires ≥2 concordant
markers per side and markers on both sides.
