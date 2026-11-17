"""Join every stated answer this project holds to the arrays it holds, in one file.

Answers come from the supplementary tables of three published studies and from the laboratory's own
sheets. Nothing here is produced by OriginMarker: each row is what a table states about one array.

  Jerabek 2026   sample overview (material, embryo, cell), cell ploidy, break sites
  Turocy 2026    per-chromosome karyotypes, karyotype text, break points
  Zuccaro 2020   per-chromosome karyotypes with the parent named (audit/zuccaro-table-s3.tsv),
                 and the GEO series' own material and pronucleus labels
  laboratory     zygote sheets (pronucleus labels, karyotype text) and manual calls

Writes, under --out (untracked: these carry laboratory identifiers):
  arrays.tsv   one row per array on disk: id, series, the file chosen, how many copies
  truth.tsv    one row per stated fact: array, source, table, kind, chrom, start, end, value, verbatim
  units.tsv    one row per array with a fact: its unit and whether it is open, sealed or exposed
  unmatched.tsv  table rows that name no array on disk

OPEN AND SEALED. A unit is an experimental group (Jerabek) or an embryo group (Turocy). A unit is
OPEN when any array of it belongs to a series this tool was developed or scored on (--burned, by
pseudonym, and all of GSE148488); every other unit is assigned by the parity of a hash of its name,
so the assignment cannot depend on anything the tool reports. Sealed arrays are scored once per
major version and never used to set a constant.

EXPOSED. An array listed in --exposed (an untracked file, one array per line) is neither open nor
sealed: it is never scored, of either kind. It is for an array of a sealed unit that a fixed test
of the tool runs, such as a bundled example file: the tool has been run on it, so it cannot count
as held out, and its unit was never open, so it cannot be developed on either.

Run: python3 audit/truth/build.py --supplements <dir> --inventory <tsv> --manifest <csv> \
       --lab-sheets <csv> --manual <tsv> --barcode-map <tsv> --exposed <list> --out audit/local/truth
"""
import argparse
import collections
import csv
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xlsx import sheets  # noqa: E402

BARCODE = re.compile(r'^(?:GSM\d+_)?(\d{5}-\d{2})(?!\d)')
ARRAY_EXT = re.compile(r'\.CEL\.probes(\.txt)?(\.gz)?$|\.probes(\.gz)?$', re.I)
# Preference among copies of one array: the laboratory folders first, the GEO re-deposit last.
FOLDER_RANK = ['DIETER', 'TREFF', 'ROBLES', 'JENNA', 'GSE186407', 'GSE315671', 'trios']
# Embryo B2's two polar bodies: Table S2 names them by tube, the files by index.
B2_POLAR_BODY = {'27': '1', '28': '2'}
# One file is named G1 where GEO titles it embryo G2, blastomere 6. The series has no embryo G1.
FILE_ALIAS = {'g1:6': 'g2:6'}


# ---------------------------------------------------------------- identifiers
def turocy_key(name):
    """A sample of the Turocy series, from a table name or a file name, as embryo:cell.

    The tables write `A1-1_16p35Xp` and `Z10_01 16p35cq`; the files put a tube number first or the
    array position last. Both carry the embryo and the cell or tube, and nothing else agrees.
    """
    s = name.strip().replace('isolated nucleus', 'nucleus')
    s = re.sub(r'\.CEL$', '', s, flags=re.I)
    z_file = re.match(r'^(\d+)_(?:chr16)?p(?:cq)?_z(\d+)(?:_|$)', s, re.I)
    if z_file:
        return f'z{int(z_file.group(2))}:{int(z_file.group(1))}'
    z_tab = re.match(r'^Z(\d+)_(\d+)\b', s)
    if z_tab:
        return f'z{int(z_tab.group(1))}:{int(z_tab.group(2))}'
    pb = re.match(r'^(?:(\d+)_)?([A-J]\d)(?:-1)?[-_]PB(?:[_-]?0?([12]))?(?=_|$)', s)
    if pb:
        tube, embryo, index = pb.groups()
        if embryo == 'B2' and not index:
            index = B2_POLAR_BODY.get(tube)
        return f'{embryo.lower()}:pb{index or ""}'
    s = re.sub(r'^\d+_(?=[A-J]\d)', '', s)
    m = re.match(r'^([A-J]\d)[-_](nucleus_\d+|PB[_\w]*|\d+|[fF])(?![\w-]*PB)', s)
    if not m:
        return None
    cell = m.group(2).lower()
    cell = re.sub(r'^nucleus_0*(\d+)$', r'nucleus\1', cell)
    cell = re.sub(r'^0+(\d)', r'\1', cell)
    return f'{m.group(1).lower()}:{cell}'


def load_inventory(path, pseud):
    """{array id: {'paths': [...], 'series': ...}} from a listing of `path<TAB>bytes`."""
    arrays = collections.OrderedDict()
    for line in open(path):
        p = line.rstrip('\n').split('\t')[0]
        if not p or not ARRAY_EXT.search(p):
            continue
        base = os.path.basename(p)
        stem = ARRAY_EXT.sub('', base)
        top = next((f for f in FOLDER_RANK if f'/{f}/' in p), 'other')
        m = BARCODE.match(base)
        gsm = re.match(r'^(GSM\d+)$', stem)
        if m:
            aid, series = m.group(1), pseud.get(m.group(1).split('-')[0], 'lab-unmapped')
        elif gsm:
            aid, series = gsm.group(1), 'GSE148488'
        elif top in ('JENNA', 'GSE186407'):
            key = turocy_key(re.sub(r'^GSM\d+_', '', stem))
            key = FILE_ALIAS.get(key, key)
            aid, series = (f'turocy:{key}' if key else f'jenna:{stem}'), 'GSE186407'
        else:
            aid, series = f'file:{stem}', 'unnamed'
        a = arrays.setdefault(aid, {'paths': [], 'series': series})
        a['paths'].append((FOLDER_RANK.index(top) if top in FOLDER_RANK else 99, p))
    for a in arrays.values():
        a['paths'].sort()
    return arrays


# ---------------------------------------------------------------- what a table cell says
def chrom_state(text):
    """(class, scope) of one per-chromosome karyotype cell. The cell's own words are kept beside it."""
    t = text.lower().strip()
    if not t or t in ('na', 'nd', '/'):
        return 'blank', ''
    if re.fullmatch(r'normal(\s*\d+)?', t):
        return 'normal', ''
    t = re.sub(r'\bnormal\b[^;,(]*', '', t)
    gain = bool(re.search(r'tri?s[io]o?my|tirsomy|tetrasom|dup|gain|\+\s*\d', t))
    loss = bool(re.search(r'monoso?my|nullisom|del|loss|(?<![\w+])-\s*\d|fragment only', t))
    whole = bool(re.search(r'tri?s[io]o?my|tirsomy|tetrasom|monoso?my|nullisom(?:y|ic)(?!\s*\d*[pq])', t))
    part = bool(re.search(r'del|dup|[pq]\b|[pq]\d|segment|telomer|fragment|@', t))
    scope = 'whole' if whole and not part else 'segmental' if part and not whole else 'both' if whole else ''
    if gain and loss:
        return 'gain+loss', scope
    if gain:
        return 'gain', scope
    if loss:
        return 'loss', scope
    return 'other', scope


def interval(text):
    """(chrom, start, end) from `chr1:55,529,828-55,534,271`, or None."""
    m = re.search(r'chr\s*(\d+|X|Y)\s*:\s*([\d,]+)\s*-\s*([\d,]+)', text, re.I)
    if not m:
        return None
    return m.group(1).upper(), int(m.group(2).replace(',', '')), int(m.group(3).replace(',', ''))


def jerabek_description(d):
    """Material, embryo number, cell number and cells in the embryo, from a sample description."""
    out = {}
    e = re.search(r'embryo\s+#(\d+)', d)
    if e:
        out['embryo'] = e.group(1)
    b = re.match(r'Blastomere #(\d+) from (?:(?:\d+ )?(~?\d+)-C|morula-stage) human embryo', d)
    if b:
        out.update(material='blastomere', cell=b.group(1))
        if b.group(2):
            out['n_cells'] = b.group(2).lstrip('~')
    elif re.match(r'1-C human embryo', d):
        out['material'] = 'zygote'
    elif re.match(r'Pronucleus #(\d+)', d):
        out.update(material='pronucleus', cell=re.match(r'Pronucleus #(\d+)', d).group(1))
    elif d.startswith('TE biopsy'):
        out['material'] = 'trophectoderm-biopsy'
    elif d.startswith('TE cells'):
        out['material'] = 'trophectoderm-outgrowth'
    elif d.startswith('ES cells'):
        out['material'] = 'esc-cells'
    elif d.startswith('gDNA isolated'):
        out['material'] = 'bulk-gdna'
    else:
        out['material'] = 'multi-cell'
    return out


def header_row(rows, first):
    return next(i for i, r in enumerate(rows) if r and r[0].strip() == first)


def cols(header):
    return {h.strip(): i for i, h in enumerate(header) if h.strip()}


def get(row, i):
    return row[i].strip() if i is not None and i < len(row) else ''


# ---------------------------------------------------------------- the build
def build(a):
    pseud = {}
    for line in open(a.barcode_map):
        p = line.rstrip('\n').split('\t')
        if len(p) == 2 and p[0] != 'prefix':
            pseud[p[0]] = p[1]
    arrays = load_inventory(a.inventory, pseud)
    truth, unmatched, unit_of = [], [], {}

    def fact(aid, source, table, kind, value, verbatim='', chrom='', start='', end=''):
        truth.append([aid, source, table, kind, chrom, start, end, value, verbatim])

    def have(aid, source, table, name):
        if aid in arrays:
            return True
        unmatched.append([source, table, name])
        return False

    S = a.supplements

    # Jerabek 2026: sample overview
    t = sheets(os.path.join(S, 'jerabek2026-nature/MOESM4_ESM.xlsx'))['Sample Overview_Sanger']
    h = header_row(t, 'Sample Barcode')
    c = cols(t[h])
    for r in t[h + 1:]:
        aid = get(r, 0)
        if not BARCODE.match(aid) or not have(aid, 'jerabek2026', 'T2', aid):
            continue
        d = get(r, c['Sample Description'])
        group = get(r, c['Experimental Group'])
        info = jerabek_description(d)
        fact(aid, 'jerabek2026', 'T2', 'material', info['material'], d)
        fact(aid, 'jerabek2026', 'T2', 'group', group)
        if 'embryo' in info:
            fact(aid, 'jerabek2026', 'T2', 'embryo', f"jerabek:{group}:{info['embryo']}", d)
        for k in ('cell', 'n_cells'):
            if k in info:
                fact(aid, 'jerabek2026', 'T2', k, info[k])
        unit_of[aid] = f'jerabek:{group}'

    # Jerabek 2026: cell ploidy
    t = sheets(os.path.join(S, 'jerabek2026-nature/MOESM5_ESM.xlsx'))['ED Fig 2c_Source Data']
    h = header_row(t, 'Single Cell Sample Barcode')
    c = cols(t[h])
    for r in t[h + 1:]:
        aid, p = get(r, 0), get(r, c['Cell Ploidy']).lower()
        if BARCODE.match(aid) and p and have(aid, 'jerabek2026', 'T3', aid):
            fact(aid, 'jerabek2026', 'T3', 'ploidy', p)

    # Jerabek 2026: break sites
    for sheet in ('Breakages_CRISPR-Cas9', 'Breakages_ABE'):
        t = sheets(os.path.join(S, 'jerabek2026-nature/MOESM7_ESM.xlsx'))[sheet]
        h = header_row(t, 'Blastomere Barcode')
        coord = next(i for i, x in enumerate(t[h]) if x.startswith('Genomic Coordinates Of Chromosome Break'))
        c = cols(t[h])
        for r in t[h + 1:]:
            aid = get(r, 0)
            if not BARCODE.match(aid) or not have(aid, 'jerabek2026', 'T5', aid):
                continue
            where = get(r, coord)
            iv = interval(where)
            chrom = iv[0] if iv else get(r, c['Chromosome']).replace('chr', '')
            note = f"{get(r, c['Number Of Chromosome Breaks'])}; {get(r, c['Method Of Mapping'])}; {where}"
            fact(aid, 'jerabek2026', f'T5:{sheet}', 'break', get(r, c['Breakage Category']), note,
                 chrom, iv[1] if iv else '', iv[2] if iv else '')

    # Turocy 2026: per-chromosome karyotypes
    t = sheets(os.path.join(S, 'turocy2026-natcommun/MOESM10_ESM.xlsx'))['Table S2']
    head = [x.replace('\n', ' ').strip() for x in t[0]]
    chrom_cols = [(i, re.match(r'^chromosome\s*(\d+)$', x, re.I).group(1))
                  for i, x in enumerate(head) if re.match(r'^chromosome\s*(\d+)$', x, re.I)]
    kary = next(i for i, x in enumerate(head) if x.lower().startswith('karotype'))
    tissue = next(i for i, x in enumerate(head) if x.lower().startswith('tissue'))
    for r in t[1:]:
        name = get(r, 0)
        key = turocy_key(name)
        if not key or len(r) <= kary:
            continue
        aid = f'turocy:{key}'
        if not have(aid, 'turocy2026', 'S2', name):
            continue
        embryo = get(r, 2) or key.split(':')[0].upper()
        fact(aid, 'turocy2026', 'S2', 'material', get(r, tissue).lower(), name)
        fact(aid, 'turocy2026', 'S2', 'embryo', f'turocy:{embryo.lower()}')
        fact(aid, 'turocy2026', 'S2', 'karyotype_text', '', get(r, kary))
        unit_of[aid] = f"turocy:{re.match(r'[a-z]+', key).group(0)}"
        if len(r) > chrom_cols[-1][0] - 5:
            for i, ch in chrom_cols:
                cls, scope = chrom_state(get(r, i))
                if cls != 'blank':
                    fact(aid, 'turocy2026', 'S2', 'chrom_state', f'{cls}:{scope}' if scope else cls,
                         get(r, i), ch)

    # Turocy 2026: karyotype text of the remaining samples
    t = sheets(os.path.join(S, 'turocy2026-natcommun/MOESM11_ESM.xlsx'))['Table S3']
    for r in t[1:]:
        name, key = get(r, 0), turocy_key(get(r, 0))
        if key and have(f'turocy:{key}', 'turocy2026', 'S3', name):
            aid = f'turocy:{key}'
            fact(aid, 'turocy2026', 'S3', 'material', get(r, 1).lower(), name)
            fact(aid, 'turocy2026', 'S3', 'karyotype_text', '', get(r, 3))
            unit_of.setdefault(aid, f"turocy:{re.match(r'[a-z]+', key).group(0)}")

    # Turocy 2026: break points. A sample's further breaks are on rows with no name.
    t = sheets(os.path.join(S, 'turocy2026-natcommun/MOESM12_ESM.xlsx'))['Table S4']
    aid = None
    for r in t[1:]:
        name = get(r, 0)
        if name:
            key = turocy_key(name)
            aid = f'turocy:{key}' if key and have(f'turocy:{key}', 'turocy2026', 'S4', name) else None
        if not aid:
            continue
        unit_of.setdefault(aid, f"turocy:{re.match(r'[a-z]+', aid.split(':')[1]).group(0)}")
        where = get(r, 3)
        iv = interval(where)
        if iv:
            fact(aid, 'turocy2026', 'S4', 'break', '', f'{get(r, 5)} {get(r, 6)}'.strip(), *iv)
        elif where.lower() == 'na' and name:
            fact(aid, 'turocy2026', 'S4', 'no_break', '')

    # Zuccaro 2020: the per-chromosome table already matched to arrays, and the series' own labels
    for line in open(a.zuccaro):
        p = line.rstrip('\n').split('\t')
        if not p[0].startswith('GSM') or not have(p[0], 'zuccaro2020', 'S3', p[1]):
            continue
        cls, scope = chrom_state(p[4])
        fact(p[0], 'zuccaro2020', 'S3', 'chrom_state', f'{cls}:{scope}' if scope else cls, p[4], p[3])
        par = re.search(r'\((mat|pat)', p[4].lower())
        if par:
            fact(p[0], 'zuccaro2020', 'S3', 'parent', 'maternal' if par.group(1) == 'mat' else 'paternal',
                 p[4], p[3])
    for r in csv.DictReader(open(a.manifest)):
        if r['gsm'] not in arrays:
            continue
        fact(r['gsm'], 'gse148488', 'GEO', 'material', r['material'], r['title'])
        if r['pronucleus']:
            fact(r['gsm'], 'gse148488', 'GEO', 'pn_parent', r['pronucleus'])
        for role in ('father', 'mother'):
            if r[f'{role}_gsms']:
                fact(r['gsm'], 'gse148488', 'GEO', role, r[f'{role}_gsms'])
        unit_of[r['gsm']] = 'gse148488'

    # The laboratory's sheets
    for r in csv.DictReader(open(a.lab_sheets)):
        aid = r['barcode'].strip()
        if not BARCODE.match(aid) or not have(aid, 'lab', 'zygote-sheets', aid):
            continue
        label = r['sheet_label'].strip()
        if label in ('M', 'P'):
            fact(aid, 'lab', 'zygote-sheets', 'pn_parent', 'maternal' if label == 'M' else 'paternal')
        elif label:
            fact(aid, 'lab', 'zygote-sheets', 'sheet_label', label)
        if r['sheet_karyotype'].strip():
            fact(aid, 'lab', 'zygote-sheets', 'karyotype_text', '', r['sheet_karyotype'].strip())
        fact(aid, 'lab', 'zygote-sheets', 'embryo', f"lab:{r['experiment'].strip()}:{r['zygote'].strip()}")
        unit_of.setdefault(aid, f"lab:{r['experiment'].strip()}")
    for r in csv.DictReader(open(a.manual), delimiter='\t'):
        aid = r['sample'].strip()
        if BARCODE.match(aid) and have(aid, 'lab', 'manual-calls', aid):
            fact(aid, 'lab', 'manual-calls', 'manual_call', r['class'], r['source_text'], r['chrom'])
            unit_of.setdefault(aid, 'lab:manual-calls')

    # Open or sealed, by unit
    burned = set(a.burned.split(',')) | {'GSE148488'}
    open_units = {u for aid, u in unit_of.items() if arrays[aid]['series'] in burned}
    open_units |= {u for u in unit_of.values() if u.startswith(('lab:', 'gse148488'))}

    exposed = {line.strip() for line in open(a.exposed) if line.strip()}

    def split(unit, aid=None):
        if aid in exposed:
            return 'exposed'
        if unit in open_units:
            return 'open'
        return 'sealed' if int(hashlib.sha1((a.salt + unit).encode()).hexdigest(), 16) % 2 else 'open'

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, 'arrays.tsv'), 'w') as f:
        f.write('array\tseries\tcopies\tpath\n')
        for aid, v in arrays.items():
            f.write(f"{aid}\t{v['series']}\t{len(v['paths'])}\t{v['paths'][0][1]}\n")
    with open(os.path.join(a.out, 'truth.tsv'), 'w') as f:
        f.write('array\tsource\ttable\tkind\tchrom\tstart\tend\tvalue\tverbatim\n')
        for row in truth:
            f.write('\t'.join(str(x).replace('\t', ' ').replace('\n', ' ') for x in row) + '\n')
    with open(os.path.join(a.out, 'units.tsv'), 'w') as f:
        f.write('array\tseries\tunit\tsplit\n')
        for aid in sorted(unit_of):
            f.write(f"{aid}\t{arrays[aid]['series']}\t{unit_of[aid]}\t{split(unit_of[aid], aid)}\n")
    with open(os.path.join(a.out, 'unmatched.tsv'), 'w') as f:
        f.write('source\ttable\tname\n')
        for row in unmatched:
            f.write('\t'.join(row) + '\n')
    return arrays, truth, unit_of, split, unmatched


def report(arrays, truth, unit_of, split, unmatched):
    by_kind = collections.defaultdict(lambda: collections.Counter())
    arrays_by_kind = collections.defaultdict(lambda: collections.defaultdict(set))
    for aid, source, table, kind, chrom, start, end, value, verbatim in truth:
        s = split(unit_of[aid], aid)
        label = kind
        if kind == 'chrom_state':
            label = 'chrom_state ' + value.split(':')[0]
        by_kind[label][s] += 1
        arrays_by_kind[label][s].add(aid)
    print(f'arrays on disk: {len(arrays)}; with at least one stated fact: {len(unit_of)}')
    n = collections.Counter(split(u, aid) for aid, u in unit_of.items())
    print(f"  open {n['open']}, sealed {n['sealed']}, exposed {n['exposed']}")
    print(f"{'fact':28s} {'open rows':>10s} {'arrays':>7s} {'sealed rows':>12s} {'arrays':>7s}")
    for k in sorted(by_kind):
        print(f"{k:28s} {by_kind[k]['open']:10d} {len(arrays_by_kind[k]['open']):7d} "
              f"{by_kind[k]['sealed']:12d} {len(arrays_by_kind[k]['sealed']):7d}")
    units = collections.defaultdict(set)
    for aid, u in unit_of.items():
        units[u].add(aid)
    print('units:')
    for u in sorted(units):
        print(f'  {split(u):6s} {u:34s} {len(units[u]):4d} arrays')
    miss = collections.Counter((s, t) for s, t, _ in unmatched)
    print('table rows naming no array on disk:', dict(miss))


def self_check():
    # Table names. File names carry a tube number or an array position and are sample identifiers,
    # so their cases are read from a file kept out of the repository where it is present.
    assert turocy_key('A1-1_16p35Xp') == 'a1:1'
    assert turocy_key('Z10_01 16p35cq') == 'z10:1'
    assert turocy_key('Z16_42 16p35') == 'z16:42'
    assert turocy_key('Z15_38 16p35') == 'z15:38'
    assert turocy_key('D3_3_16qXq') == 'd3:3'
    assert turocy_key('75_D3-f_16qXq') == 'd3:f'
    assert turocy_key('I2_01_XDMD_Chr16pprox') == 'i2:1'
    assert turocy_key('I4_isolated nucleus_01_XDMD_Chr16pprox') == 'i4:nucleus1'
    assert turocy_key('H3_10_Xp') == 'h3:10'
    assert turocy_key('A1-1_PB') == 'a1:pb'
    assert turocy_key('27_B2_PB') == 'b2:pb1'
    assert turocy_key('28_B2_PB') == 'b2:pb2'
    assert turocy_key('C1_PB_02') == 'c1:pb2'
    assert turocy_key('D4_PB01') == 'd4:pb1'
    assert turocy_key('C3_PB') == 'c3:pb'
    assert turocy_key('E1-PB') == 'e1:pb'
    assert turocy_key('I2_PB') == 'i2:pb'
    assert turocy_key('Z10_02') == 'z10:2'
    assert turocy_key('Z12_18 16p35cq') == 'z12:18'
    assert turocy_key('total: 127') is None
    fixtures = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'local', 'truth', 'name-fixtures.tsv')
    if os.path.exists(fixtures):
        for line in open(fixtures):
            name, key = line.rstrip('\n').split('\t')
            assert turocy_key(name) == key, name
    assert chrom_state('normal') == ('normal', '') and chrom_state('') == ('blank', '')
    assert chrom_state('trisomy 8') == ('gain', 'whole') and chrom_state('tirsomy 8') == ('gain', 'whole')
    assert chrom_state('monosomy 16') == ('loss', 'whole') and chrom_state('nullisomy 3') == ('loss', 'whole')
    assert chrom_state('del 16q') == ('loss', 'segmental') and chrom_state('dup 17q') == ('gain', 'segmental')
    assert chrom_state('-16q') == ('loss', 'segmental') and chrom_state('+17q') == ('gain', 'segmental')
    assert chrom_state('dup 16p, monosomy 16')[0] == 'gain+loss'
    assert chrom_state('monosomy 6 (pat)') == ('loss', 'whole')
    assert chrom_state('nullisomic') == ('loss', 'whole') and chrom_state('nullisomy 17p')[1] == 'segmental'
    assert chrom_state('normal; dup 16q') == ('gain', 'segmental')
    assert chrom_state('fragmentation')[0] == 'other'
    assert interval('chr1:55,529,828-55,534,271') == ('1', 55529828, 55534271)
    assert interval('chr16:3417931-3419719') == ('16', 3417931, 3419719) and interval('centromere') is None
    d = jerabek_description('Blastomere #2 from ~12-C human embryo #4')
    assert d == {'embryo': '4', 'material': 'blastomere', 'cell': '2', 'n_cells': '12'}
    assert jerabek_description('1-C human embryo #3') == {'embryo': '3', 'material': 'zygote'}
    assert jerabek_description('TE biopsy from 5AA human embryo')['material'] == 'trophectoderm-biopsy'
    assert jerabek_description('Blastomere #1 from morula-stage human embryo #2')['cell'] == '1'
    print('build self-check ok')


if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == '--check':
        self_check()
        sys.exit(0)
    ap = argparse.ArgumentParser()
    for name in ('supplements', 'inventory', 'manifest', 'zuccaro', 'lab-sheets', 'manual', 'barcode-map',
                 'exposed', 'out'):
        ap.add_argument(f'--{name}', required=True)
    ap.add_argument('--burned', default='S12,S13,S20,S21,S23,S24',
                    help='series pseudonyms the tool was developed or scored on')
    ap.add_argument('--salt', default='originmarker-sealed-v1')
    report(*build(ap.parse_args()))
