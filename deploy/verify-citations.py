import json, urllib.request, time

FINAL = {
    "karyomapping":  "10.1136/jmg.2009.069971",
    "eshre_pgt_m":   "10.1093/hropen/hoaa018",
    "gnomad_v4":     "10.1038/s41586-023-06045-0",
    "clinvar":       "10.1093/nar/gkz972",
    "ensembl":       "10.1093/nar/gkad1049",
    "ensembl_rest":  "10.1093/bioinformatics/btu613",
    "decode_map":    "10.1126/science.aau1043",
    "thousand_g":    "10.1038/nature15393",
    "ldlink":        "10.1093/bioinformatics/btv402",
    "hgvs":          "10.1002/humu.22981",
    "dbsnp":         "10.1093/nar/29.1.308",
    "abcc8_chi":     "10.1186/1750-1172-6-63",
    "ado":           "10.1002/pd.109",
    "kosambi":       "10.1111/j.1469-1809.1943.tb02321.x",
}

def get(doi):
    req = urllib.request.Request(f"https://api.crossref.org/works/{doi}",
        headers={"User-Agent": "OriginMarker/0.1 (mailto:kruger.ezra.s@gmail.com)"})
    with urllib.request.urlopen(req, timeout=30) as f:
        return json.load(f)["message"]

def authors(m):
    a = m.get("author") or []
    if not a:
        return None
    fams = [x.get("family") for x in a if x.get("family")]
    if not fams:
        return None
    if len(fams) == 1:
        return fams[0]
    if len(fams) == 2:
        return f"{fams[0]} & {fams[1]}"
    return f"{fams[0]} et al."

out = {}
for k, doi in FINAL.items():
    m = get(doi)
    yr = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    out[k] = {
        "authors": authors(m) or "Consortium",
        "year": yr,
        "title": (m.get("title") or ["?"])[0],
        "journal": (m.get("container-title") or ["?"])[0],
        "volume": m.get("volume") or "",
        "page": m.get("page") or "",
        "doi": doi,
        "url": f"https://doi.org/{doi}",
    }
    print(f"{k:14} {out[k]['authors']} ({yr}) {out[k]['journal']} {out[k]['volume']}")
    print(f"{'':14} {out[k]['title'][:78]}")
    time.sleep(0.4)

# Haldane 1919 predates DOIs. Cited classically rather than with a fabricated identifier.
out["haldane"] = {
    "authors": "Haldane", "year": 1919,
    "title": ("The combination of linkage values, and the calculation of distances "
              "between the loci of linked factors"),
    "journal": "Journal of Genetics", "volume": "8", "page": "299-309",
    "doi": None, "url": None,
    "note": "Predates DOI assignment; no persistent identifier exists for the original.",
}
print(f"\nhaldane        cited classically (no DOI exists - 1919)")

json.dump(out, open("/tmp/citations.json", "w"), indent=1)
print(f"\n{len(out)} citations locked -> /tmp/citations.json")
