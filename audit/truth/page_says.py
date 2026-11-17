"""What the result card says, measured on the 31 published cells of Zuccaro 2020 Table S3.

Two outputs, both from records and the answer key, never from an array:

  the headline each cell gets in each setup, with the inputs that decide it, written as a fixture
  that web/src/pageSays.check.ts replays through headlineOf;

  the measured count behind each chip word: per kind of call and setup, how many printed calls of
  that kind agree with the paper, by the plan's counting rule (per array, chromosome and loss or
  gain, the joined calls and the joined stated region overlap by half of the longer total).

Inputs: run.ts records with tags two, father and recon (the setups), the change lists the page
itself builds with both parents loaded (cli/cardrun.ts output, one JSON card per line), and
audit/local/truth/truth.tsv. The kind of call follows pageSays.kindOf; keep the two in step.

Run: python3 audit/truth/page_says.py --truth audit/local/truth --records h.jsonl --cards cards.jsonl
       [--fixture audit/truth/page-says-cells.json]
"""
import argparse
import collections
import csv
import json
import os
import re
import sys

CEN = {'1': 125, '2': 93.3, '3': 91, '4': 50.4, '5': 48.4, '6': 61, '7': 59.9, '8': 45.6, '9': 49,
       '10': 40.2, '11': 53.7, '12': 35.8, '13': 17.9, '14': 17.6, '15': 19, '16': 36.6, '17': 24,
       '18': 17.2, '19': 26.5, '20': 27.5, '21': 13.2, '22': 14.7}
END = {'1': 249.3, '2': 243.2, '3': 198, '4': 191.2, '5': 180.9, '6': 171.1, '7': 159.1, '8': 146.4,
       '9': 141.2, '10': 135.5, '11': 135, '12': 133.9, '13': 115.2, '14': 107.3, '15': 102.5,
       '16': 90.4, '17': 81.2, '18': 78.1, '19': 59.1, '20': 63, '21': 48.1, '22': 51.3}
# GRCh37 band starts used where the paper names a band to the telomere.
BAND = {'16q23.1': 74.1, '16q21': 64.0, '6q21': 105.0, '5q32': 145.0, '1q25.2': 175.0}
SETUPS = (('two', 'both parents'), ('father', 'measured father'), ('recon', 'rebuilt father'))
OTHER = {'paternal': 'maternal', 'maternal': 'paternal'}


def load(a):
    recs = collections.defaultdict(dict)
    for line in open(a.records):
        if line.startswith('{'):
            r = json.loads(line)
            recs[r['tag']][r['id']] = r
    cards = {}
    for line in open(a.cards):
        if line.startswith('{'):
            c = json.loads(line)
            cards[c['name'].split('.')[0]] = c
    state, parent, verb, split = {}, {}, {}, {}
    for r in csv.DictReader(open(os.path.join(a.truth, 'units.tsv')), delimiter='\t'):
        split[r['array']] = r['split']
    for r in csv.DictReader(open(os.path.join(a.truth, 'truth.tsv')), delimiter='\t'):
        if r['array'] in cards:
            if r['kind'] == 'chrom_state':
                state[(r['array'], r['chrom'])] = r['value']
                verb[(r['array'], r['chrom'])] = r['verbatim']
            elif r['kind'] == 'parent':
                parent[(r['array'], r['chrom'])] = r['value']
    if any(split.get(i) == 'sealed' for i in cards):
        sys.exit('sealed arrays in these cards. Refused.')
    return recs, cards, state, parent, verb


def stated_extent(state, verb, i, ch, cls):
    v, t = state[(i, ch)], verb[(i, ch)]
    if cls not in v and not (cls == 'loss' and v.startswith('other:segmental')):
        return None
    if 'whole' in v or 'only' in t or 'Chr16 gain' in t and cls == 'gain':
        return (0, END[ch])
    part = t
    if '+' in v:
        bits = re.split(r',\s*', t)
        part = next((b for b in bits if (('del' in b) if cls == 'loss' else ('dup' in b or 'gain' in b))), t)
    for b, mb in BAND.items():
        if b in part:
            return (mb, END[ch])
    arm = re.search(r'\d([pq])', part.replace(' ', '').replace('(', ''))
    if arm:
        return (0, CEN[ch]) if arm.group(1) == 'p' else (CEN[ch], END[ch])
    return (0, END[ch])


def where_key(ch, s, e):
    return f'chr{ch} {s / 1e6:.1f}-{e / 1e6:.1f}Mb'


def rows_two(card, rec):
    """The page's own rows with both parents, each with its kind by pageSays.kindOf."""
    ch = card['changes']
    na = len(ch) - len(rec['findings']) - len(rec['segments'])
    one = {o['where']: o for o in rec['origins'] if o['by'] == 'one-parent'}
    out = []
    for k, d in enumerate(ch):
        named = d['origin'] in ('paternal', 'maternal')
        if k < na:
            c = next(x for x in rec['chroms'] if x['chrom'] == d['chrom'])
            if d['kind'] == 'whole-chromosome' or not (c['aneuploidy'] or c['lostParent'] or c['gain'] == 'gain'):
                kind = 'imbalance'  # an allele imbalance, or a mixed signal (not in the records)
            elif d['kind'] == 'copy-gain' and named:
                kind = 'gain-parent'
            elif c['lostParent'] or (named and one.get(f"chr{d['chrom']}")):
                kind = 'whole' if named else 'intensity'
            else:
                kind = 'intensity'
        elif k < na + len(rec['segments']):
            s = rec['segments'][k - na]
            o = one.get(where_key(d['chrom'], d['startBp'], d['endBp']))
            if s['kind'] != 'parental-absence':
                kind = 'gain-parent' if s['kind'] == 'copy-gain' and named else 'intensity'
            elif s.get('lrrShift') is not None:
                kind = 'allele-scan'
            elif not named:
                kind = 'unnamed'
            elif s.get('absentParent'):
                kind = 'two-ways' if (o and o.get('parent') == d['origin'] and o['verdict'] == 'parent-lost') else 'run-only'
            else:
                kind = 'second-test'
        else:
            kind = 'finding'
        out.append(dict(chrom=d['chrom'], start=d['startBp'] / 1e6, end=d['endBp'] / 1e6, kind=kind,
                        origin=d['origin'], gain=d['kind'] == 'copy-gain'))
    return out


def rows_one(rec, role='paternal'):
    """The rows the page builds with one parent, from a record, by the defectsForResult rules."""
    one = [o for o in rec['origins'] if o['by'] == 'one-parent']
    dos = [o for o in rec['origins'] if o['by'] == 'dosage']

    def span(o):
        m = re.match(r'chr(\w+) ([\d.]+)-([\d.]+)Mb', o['where'])
        return (m.group(1), float(m.group(2)), float(m.group(3))) if m else (None, None, None)

    def named(v):
        return role if v in ('loaded-parent', 'known-parent-lost') else OTHER[role] if v in ('other-parent', 'other-parent-lost') else None

    def ov(a, b, c, d):
        return max(0, min(b, d) - max(a, c))
    out = []
    for c in rec['chroms']:
        gain = c['aneuploidy'] == 'gain'
        if not (c['aneuploidy'] or c['lostParent'] or c['gain'] in ('gain', 'imbalance')):
            continue
        mendel = next((o for o in one if o['where'] == f"chr{c['chrom']}"), None)
        by_mendel = named(mendel['verdict']) if mendel and not gain else None
        letters = (OTHER[c['lostParent']] if gain else c['lostParent']) if c['lostParent'] else None
        gp = c['gainParent'] if c['gain'] == 'gain' else None
        scored = next((o for o in dos if o['where'] == f"chr{c['chrom']}"), None)
        dp = named(scored['verdict']) if scored else None
        if (by_mendel and letters and by_mendel != letters) or (gp and gain and dp and dp != gp) \
                or (c['gain'] == 'imbalance' and not c['aneuploidy'] and not c['lostParent']):
            origin = 'unclear'
        elif gp and not letters:
            origin = gp
        else:
            origin = by_mendel or letters or (role if c['aneuploidyParent'] == 'this' else OTHER[role]
                                               if c['aneuploidyParent'] == 'other' else dp) or 'unclear'
        is_named = origin != 'unclear'
        if not c['aneuploidy'] and not c['lostParent'] and c['gain'] == 'imbalance':
            kind = 'imbalance'
        elif (gain or (not c['aneuploidy'] and not c['lostParent'] and gp)) and is_named:
            kind = 'gain-parent'
        elif (c['lostParent'] or by_mendel) and is_named:
            kind = 'whole'
        else:
            kind = 'intensity'
        out.append(dict(chrom=c['chrom'], start=0, end=END[c['chrom']], kind=kind, origin=origin, gain=gain))
    for s in rec['segments']:
        a, b = s['start'] / 1e6, s['end'] / 1e6
        origin = 'unclear'
        if s['kind'] != 'parental-absence':
            kind = 'intensity'
        elif s.get('lrrShift') is not None:
            kind, origin = 'allele-scan', s.get('absentParent') or role
        else:
            o = max((o for o in one if span(o)[0] == s['chrom'] and span(o)[1] is not None
                     and ov(a, b, span(o)[1], span(o)[2]) > 0), key=lambda o: ov(a, b, span(o)[1], span(o)[2]), default=None)
            p = named(o['verdict']) if o else None
            kind, origin = ('second-test', p) if p else ('unnamed', 'unclear')
        if kind == 'intensity':
            d = max((o for o in dos if span(o)[0] == s['chrom'] and span(o)[1] is not None
                     and ov(a, b, span(o)[1], span(o)[2]) > 0), key=lambda o: ov(a, b, span(o)[1], span(o)[2]), default=None)
            if d and named(d['verdict']):
                origin = named(d['verdict'])
                if s['kind'] == 'copy-gain':
                    kind = 'gain-parent'
        out.append(dict(chrom=s['chrom'], start=a, end=b, kind=kind, origin=origin, gain=s['kind'] == 'copy-gain'))
    for f in rec['findings']:
        out.append(dict(chrom=f['chrom'], start=(f['start'] or 0) / 1e6, end=(f['end'] or 0) / 1e6, kind='finding',
                        origin='unclear', gain=False))
    return out


def judge(i, rows, state, parent, verb):
    """Place and parent for each row against the paper."""
    groups = collections.defaultdict(list)
    for r in rows:
        if r['chrom'] in CEN and r['kind'] != 'imbalance':
            groups[(r['chrom'], 'gain' if r['gain'] else 'loss')].append((r['start'], r['end']))
    for r in rows:
        ch, cls = r['chrom'], 'gain' if r['gain'] else 'loss'
        if ch not in CEN:
            r['place'], r['par'] = 'unscored', '-'
            continue
        st = state[(i, ch)]
        if st == 'other':
            r['place'] = 'unscored'
        elif r['kind'] == 'imbalance':
            r['place'] = 'wrong' if st == 'normal' else 'right'
        elif st == 'normal':
            r['place'] = 'wrong'
        else:
            ext = stated_extent(state, verb, i, ch, cls)
            spans = groups[(ch, cls)]
            tot = sum(b - a for a, b in spans)
            ovl = sum(max(0, min(b, ext[1]) - max(a, ext[0])) for a, b in spans) if ext else 0
            r['place'] = 'right' if ext and ovl >= 0.5 * max(tot, ext[1] - ext[0]) else 'wrong'
        tp = parent.get((i, ch))
        if 'only' in verb.get((i, ch), ''):
            tp = None
        elif tp and '+' in st:
            part = [b for b in re.split(r',\s*', verb[(i, ch)]) if (('dup' in b or 'gain' in b) if cls == 'gain' else 'del' in b)]
            if part:
                tp = 'paternal' if '(pat)' in part[0] else 'maternal' if '(mat)' in part[0] else tp
        if r['origin'] not in ('paternal', 'maternal') or r['place'] != 'right' or tp is None:
            r['par'] = '-'
        else:
            r['par'] = 'right' if r['origin'] == tp else 'wrong'
    return rows


def headline_inputs(i, setup, rec, card):
    """What headlineOf reads, from a record. The mother's genome-wide verdict is not in these
    records: taken as not absent until the harness writes it."""
    return {
        'accession': i, 'setup': setup,
        'failed': rec['stage']['stage'] == 'failed',
        'zygosity': rec['zygosity'],
        'cls': card['origin'] if setup == 'two' else rec['originClass'],
        'loaded': [{'who': 'father', 'check': rec['parentGenome']['verdict'],
                    'present': rec['parentGenome']['present'], 'read': rec['parentGenome']['read'],
                    'absent': rec['originClass'] == 'gynogenetic'}]
        + ([{'who': 'mother', 'check': rec['parentGenomeOther']['verdict'],
             'present': rec['parentGenomeOther']['present'], 'read': rec['parentGenomeOther']['read'],
             'absent': False}] if setup == 'two' else []),
    }


def expected_word(h):
    """The rule of plan item 3.14 P1, written here once so the check has an independent oracle."""
    if h['failed']:
        return 'UNCLEAR'
    if h['zygosity'] == 'uniparental_homozygous':
        return h['cls'].replace('_', ' ').upper()
    if h['cls'] == 'neither_parent':
        return 'NEITHER PARENT'
    if any(p['check'] != 'carries' for p in h['loaded']) or any(p['absent'] for p in h['loaded']):
        return 'UNCLEAR'
    if h['zygosity'] != 'diploid':
        return 'UNCLEAR'
    return 'BOTH PARENTS' if len(h['loaded']) == 2 else 'FATHER + SECOND PARENT'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth', required=True)
    ap.add_argument('--records', required=True)
    ap.add_argument('--cards', required=True)
    ap.add_argument('--fixture')
    a = ap.parse_args()
    recs, cards, state, parent, verb = load(a)
    ids = sorted(cards)
    fixture = []
    words = collections.defaultdict(collections.Counter)
    tally = collections.defaultdict(collections.Counter)
    parent_right = collections.Counter()
    for i in ids:
        for setup, _ in SETUPS:
            rec = recs[setup][i]
            h = headline_inputs(i, setup, rec, cards[i])
            h['expected'] = expected_word(h)
            fixture.append(h)
            words[setup][h['expected']] += 1
            rows = rows_two(cards[i], rec) if setup == 'two' else rows_one(rec)
            for r in judge(i, rows, state, parent, verb):
                c = tally[(r['kind'], setup)]
                if r['place'] == 'unscored':
                    c['unscored'] += 1
                    continue
                c['scored'] += 1
                c['right'] += r['place'] == 'right' and r['par'] != 'wrong'
                c['in place'] += r['place'] == 'right'
                c['wrong parent'] += r['par'] == 'wrong'
                if r['kind'] not in ('whole', 'imbalance', 'finding') and r['par'] in ('right', 'wrong'):
                    parent_right[r['par']] += 1
    if a.fixture:
        with open(a.fixture, 'w') as f:
            json.dump(fixture, f, indent=0)
            f.write('\n')
    print(f'{len(ids)} published cells, {len(fixture)} cell-setups')
    print('\nHEADLINE, expected word by setup')
    for setup, name in SETUPS:
        print(f'  {name:16s} ' + ', '.join(f'{v} {k}' for k, v in sorted(words[setup].items())))
    print('\n  UNCLEAR, with the rule that keeps it')
    for h in fixture:
        if h['expected'] != 'UNCLEAR':
            continue
        why = ('failed array' if h['failed'] else
               'two readings disagree (parent check carries, genome-wide absent)' if any(p['absent'] for p in h['loaded']) else
               'a parent not shown' if any(p['check'] != 'carries' for p in h['loaded']) else
               'copy number not read')
        print(f"    {h['accession']}  {h['setup']:6s}  {why}")
    print('\nMEASURED COUNT BEHIND EACH CHIP WORD (right / scored; wrong parent; unscored)')
    for kind in ('whole', 'two-ways', 'run-only', 'second-test', 'allele-scan', 'unnamed', 'intensity',
                 'gain-parent', 'imbalance', 'finding'):
        cells = []
        for setup, _ in SETUPS:
            c = tally.get((kind, setup))
            cells.append(f"{c['right']:3d} of {c['scored']:<3d} (parent wrong {c['wrong parent']}, unscored {c['unscored']})"
                         if c else f"{'-':>40s}")
        print(f'  {kind:12s} ' + '   '.join(cells))
    print('\n  Machine-read by web/src/pageSays.check.ts (kind, setup, right, scored):')
    for kind in ('whole', 'two-ways', 'run-only', 'second-test', 'allele-scan', 'unnamed', 'intensity',
                 'gain-parent', 'imbalance'):
        for setup, _ in SETUPS:
            c = tally.get((kind, setup))
            if c and c['scored']:
                print(f"  MEASURED {kind} {dict(two='two', father='measured', recon='rebuilt')[setup]} {c['right']} {c['scored']}")
    print('\n  "Right" is right place and, where the call and the paper both name one, right parent.')
    print('  For unnamed, intensity and imbalance rows, right is right place.')
    print(f"\nPARENT ON A STRETCH IN THE RIGHT PLACE, all three setups: right {parent_right['right']}, "
          f"wrong {parent_right['wrong']}")


if __name__ == '__main__':
    main()
