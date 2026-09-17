"""Stage and embryo identity for the embryo series, from their own GEO metadata.

No series states either in a structured field. Every one of them writes the developmental stage and
the embryo into free text, in its own spelling, so each series needs its own reading:

    GSE148488   the sample title: "maternal nucleus isolated from 2PN zygote 3", "Cleavage stage
                blastomere_day4_... embryo 1 rep 1", "genomic DNA sperm donor rep 1". The embryo
                label repeats between egg donors, so the key carries the donor from source_name.
    GSE186407   source_name gives the material, the description gives the cell count and the embryo:
                "Cleavage stage blastomere_8 cell embryo day 3_..._embryo Z10_1".
    GSE290961   source_name gives material and injection stage, the description gives the cell count
                and an embryo number that restarts per experiment: "Blastomere #2 from 24-C human
                embryo #3 ... _Experiment BE5", so the key is experiment plus number.

Writes one row per sample: stage, cell count where stated, embryo key, and the parental role for
arrays that are a parent rather than a sample.

Run: python geo_sheet.py <GSE> [<GSE> ...] > sheet.json
"""
import json
import re
import sys
import urllib.request

URL = ('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi'
       '?acc={acc}&targ=gsm&form=text&view=brief')


def records(acc):
    """Every sample of a series, as dicts of the GEO keys that carry text."""
    req = urllib.request.Request(URL.format(acc=acc), headers={'User-Agent': 'OriginMarker/0.1'})
    with urllib.request.urlopen(req, timeout=180) as f:
        body = f.read().decode('utf8', 'replace')
    out, cur = [], None
    for line in body.split('\n'):
        if line.startswith('^SAMPLE'):
            cur = {'gsm': line.split('=', 1)[1].strip(), 'chars': []}
            out.append(cur)
        elif cur is None or '=' not in line:
            continue
        else:
            k, v = line.split('=', 1)
            k, v = k.strip('! ').strip(), v.strip()
            if k == 'Sample_title':
                cur['title'] = v
            elif k == 'Sample_source_name_ch1':
                cur['source'] = v
            elif k == 'Sample_description':
                cur['desc'] = (cur.get('desc', '') + ' ' + v).strip()
            elif k == 'Sample_characteristics_ch1':
                cur['chars'].append(v)
    return out


def cells_of(text):
    """The embryo's cell count where the text states one, as an integer."""
    m = re.search(r'(\d+)[\s-]?(?:cell|C)\b', text, re.I)
    return int(m.group(1)) if m else None


def gse148488(r):
    t = r.get('title', '')
    src = r.get('source', '')
    donor = (re.search(r'egg donor ([A-Z])', src) or re.search(r'egg donor ([A-Z])', t))
    donor = donor.group(1) if donor else ''
    if re.search(r'genomic DNA sperm donor', t, re.I):
        return dict(stage='bulk gDNA', embryo='', role='father')
    if re.search(r'genomic DNA egg donor|somatic cells from egg|cumulus', t, re.I):
        return dict(stage='bulk gDNA', embryo='', role='mother', donor=donor)
    m = re.search(r'(maternal|paternal) nucleus isolated from 2PN zygote\s*(\d+)', t, re.I)
    if m:
        return dict(stage='pronucleus', cells=1, embryo=f'zygote {m.group(2)}',
                    parentOfOrigin=m.group(1).lower())
    if re.match(r'2PN zygote', t, re.I):
        return dict(stage='whole zygote', cells=1, embryo=t)
    emb = re.search(r'embryo\s*([A-Za-z0-9]+)', t)
    key = f"{donor}/{emb.group(1)}" if emb else ''
    if re.search(r'blastomere', t, re.I):
        return dict(stage='cleavage blastomere', embryo=key)
    if re.search(r'trophectoderm', t, re.I):
        return dict(stage='trophectoderm', embryo=key)
    if re.search(r'\bES|pluripotent|stem cell', t, re.I):
        return dict(stage='stem cell line', embryo=key)
    return dict(stage='other', embryo=key)


def gse186407(r):
    src, desc = r.get('source', ''), r.get('desc', '') or r.get('title', '')
    emb = re.search(r'embryo\s+([A-Za-z]\d+)', desc + ' ' + r.get('title', ''))
    key = emb.group(1) if emb else ''
    cells = cells_of(desc)
    if re.search(r'polar body', src, re.I):
        return dict(stage='polar body', embryo=key, cells=cells)
    if re.search(r'trophectoderm', src, re.I):
        return dict(stage='trophectoderm', embryo=key, cells=cells)
    if re.search(r'fragment', src, re.I):
        return dict(stage='cell fragment', embryo=key, cells=cells)
    if re.search(r'2PN embryo from', src, re.I):
        return dict(stage='whole zygote', embryo=key, cells=1)
    return dict(stage='cleavage blastomere', embryo=key, cells=cells)


def gse290961(r):
    src, desc = r.get('source', ''), r.get('desc', '')
    exp = re.search(r'Experiment(?:al Group)?\s+([\w-]+)', desc)
    num = re.search(r'embryo\s*#\s*(\d+)', desc, re.I)
    key = f"{exp.group(1) if exp else '?'}/{num.group(1)}" if num else ''
    cells = cells_of(desc)
    if re.search(r'ESC|stem cell|ICM', src + desc, re.I):
        return dict(stage='stem cell line', embryo=key)
    if re.search(r'trophectoderm', src, re.I):
        return dict(stage='trophectoderm', embryo=key, cells=cells)
    if re.search(r'blastocyst|morula', src + desc, re.I):
        return dict(stage='whole embryo', embryo=key, cells=cells)
    if re.search(r'zygote|\b1-C\b', src + desc, re.I) and not re.search(r'[Bb]lastomere', src):
        return dict(stage='whole zygote', embryo=key, cells=1)
    return dict(stage='cleavage blastomere', embryo=key, cells=cells)


READERS = {'GSE148488': gse148488, 'GSE186407': gse186407, 'GSE290961': gse290961}


def main():
    rows = []
    for acc in sys.argv[1:]:
        read = READERS.get(acc)
        if not read:
            raise SystemExit(f'no reading rule for {acc}')
        for r in records(acc):
            row = {'gsm': r['gsm'], 'series': acc, 'title': r.get('title', ''), 'role': 'sample'}
            row.update(read(r))
            rows.append(row)
        print(f'{acc}: {sum(1 for x in rows if x["series"] == acc)} samples', file=sys.stderr)
    json.dump(rows, sys.stdout, indent=1)


if __name__ == '__main__':
    main()
