#!/usr/bin/env python3
"""
The validation scorecard: every check of audit/validate in one report, as numbers that can be set
against the previous release's.

Sections, each run when its input is given:
  --invariance <jsonl...>   audit/validate/invariance.ts output. shuffle and flip must reproduce the
                            original exactly; half may lose a claim, never contradict one.
  --census <prefix...>      audit/refusal-census.ts TSV shards per configuration (prefix-0.tsv ...),
                            the first being the reference: claims against Table S3, and parity of
                            stage, zygosity and copy number with the reference.
  --constructed <jsonl...>  audit/validate/constructed-gains.ts output: extra chromosomes of known
                            parental origin, per configuration: named, right, wrong, and agreement
                            between the one-parent and two-parent readings.
  --heldout <jsonl...>      the same, on chromosomes no threshold was set on; reported apart.
  --previous <json>         a previous scorecard; every number that moved is listed.

Writes the scorecard as JSON (--out) and prints the report. Stdlib only.
"""
import argparse, collections, glob, json, sys


def load_jsonl(paths):
    out = []
    for p in paths:
        for line in open(p):
            line = line.strip()
            if line and line != 'DONE':
                out.append(json.loads(line))
    return out


CLAIM_FIELDS = ('aneuploidy', 'lostParent', 'aneuploidyParent', 'alleleGain')
# Verdicts that make no claim: the allele-gain test ran and found nothing, or could not run.
NOT_A_CLAIM = ('none', 'unscored')


def claim(v):
    return '' if v in NOT_A_CLAIM else v


def invariance(rows):
    by = collections.defaultdict(dict)
    for r in rows:
        by[r['gsm']][r['variant']] = r
    card = {'cells': len(by)}
    detail = []
    for v in ('shuffle', 'flip'):
        diffs = []
        for gsm, d in sorted(by.items()):
            if 'orig' not in d or v not in d:
                continue
            a, b = dict(d['orig']), dict(d[v])
            a.pop('variant'); b.pop('variant')
            bad = [k for k in a if a[k] != b.get(k)]
            if bad:
                diffs.append((gsm, bad))
                for k in bad:
                    detail.append(f'{v} {gsm} {k}: {summarise(a[k], b.get(k))}')
        card[f'{v}_cells_differing'] = len(diffs)
    contradictions, new_claims, stage_moves = [], [], []
    for gsm, d in sorted(by.items()):
        if 'orig' not in d or 'half' not in d:
            continue
        a, b = d['orig'], d['half']
        if a['stage'] != b['stage'] and 'failed' not in (a['stage'], b['stage']):
            stage_moves.append(gsm); detail.append(f'half {gsm} stage {a["stage"]}/{a["ploidy"]} -> {b["stage"]}/{b["ploidy"]}')
        if a['zygosity'] != b['zygosity'] and 'unknown' not in (a['zygosity'], b['zygosity']):
            contradictions.append((gsm, 'zygosity')); detail.append(f'half {gsm} zygosity {a["zygosity"]} -> {b["zygosity"]}')
        ca = {c[0]: c for c in a['chroms']}
        for c in b['chroms']:
            o = ca.get(c[0])
            if not o:
                continue
            for i, name in enumerate(CLAIM_FIELDS, start=2):
                o_i, c_i = claim(o[i]), claim(c[i])
                if o_i and c_i and o_i != c_i:
                    contradictions.append((gsm, f'chr{c[0]} {name}'))
                    detail.append(f'half {gsm} chr{c[0]} {name}: {o_i} -> {c_i}')
                elif c_i and not o_i:
                    new_claims.append((gsm, f'chr{c[0]} {name}'))
                    detail.append(f'half {gsm} chr{c[0]} {name}: new {c_i}')
        for kind in ('segments', 'findings'):
            for x in sorted(set(b[kind]) - set(a[kind])):
                new_claims.append((gsm, x)); detail.append(f'half {gsm} new {kind[:-1]} {x}')
    card['half_contradictions'] = len(contradictions)
    card['half_new_claims'] = len(new_claims)
    card['half_stage_moves'] = len(stage_moves)
    return card, detail


def summarise(a, b):
    if isinstance(a, list) and isinstance(b, list):
        sa, sb = set(map(json.dumps, a)), set(map(json.dumps, b))
        return f'{len(sa - sb)} only in orig, {len(sb - sa)} only in variant, e.g. {sorted(sa - sb)[:2]} / {sorted(sb - sa)[:2]}'
    return f'{a!r} -> {b!r}'


def census_rows(prefix):
    rows = []
    for f in sorted(glob.glob(prefix + '-*.tsv')):
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            if p and p[0] != 'kind':
                rows.append(p)
    return rows


def census(prefixes):
    card, detail = {}, []
    data = {p: census_rows(p) for p in prefixes}
    cell = lambda rows: {r[1]: (r[2], r[3], r[4]) for r in rows if r[0] == 'cell'}
    def copy_rows(rows):
        out = collections.defaultdict(list)
        for r in rows:
            if r[0] in ('made', 'withheld') and (r[5] == 'aneuploidy' or (r[5] == 'segment' and r[7].startswith('copy'))):
                out[r[1]].append((r[0], r[4], r[5], r[6], r[7]))
        return {k: sorted(v) for k, v in out.items()}
    for p, rows in data.items():
        name = p.split('/')[-1]
        made = [r for r in rows if r[0] == 'made']
        lp = collections.Counter(r[11] for r in made if r[5] == 'lostParent')
        card[f'{name}.lostParent'] = dict(lp)
        card[f'{name}.claims_on_published_normal'] = sum(1 for r in made if r[8] == 'normal' and r[11] == 'contradicted')
        card[f'{name}.withheld_by_gate'] = dict(collections.Counter(r[4] for r in rows if r[0] == 'withheld'))
    ref = prefixes[0]
    a, ca = cell(data[ref]), copy_rows(data[ref])
    for p in prefixes[1:]:
        b, cb = cell(data[p]), copy_rows(data[p])
        name = p.split('/')[-1]
        same = sum(1 for g in a if a[g] == b.get(g))
        card[f'parity.{name}.stage_zygosity'] = f'{same}/{len(a)}'
        for g in a:
            if a[g] != b.get(g):
                detail.append(f'parity {name} {g}: {a[g]} vs {b.get(g)}')
        cells = set(ca) | set(cb)
        card[f'parity.{name}.copy_rows'] = f'{sum(1 for g in cells if ca.get(g) == cb.get(g))}/{len(cells)}'
    return card, detail


def constructed(rows):
    card, detail = {}, []
    for cfg in ('one', 'two', 'meas1'):
        if not any(r['config'] == cfg for r in rows):
            continue
        built = [r for r in rows if r['config'] == cfg and r['constructed'] != 'none']
        named = [r for r in built if r['gain'] == 'gain']
        right = [r for r in named if r['parent'] == r['constructed']]
        card[f'{cfg}.constructed'] = len(built)
        card[f'{cfg}.named'] = len(named)
        card[f'{cfg}.right'] = len(right)
        card[f'{cfg}.wrong'] = len(named) - len(right)
        card[f'{cfg}.intensity_gain_called'] = sum(1 for r in built if r['aneuploidy'] == 'gain')
        base = [r for r in rows if r['config'] == cfg and r['constructed'] == 'none']
        card[f'{cfg}.false_gain_unmodified'] = sum(1 for r in base if r['gain'] == 'gain')
        card[f'{cfg}.other_chromosome_gains'] = sum(len(r['otherGains']) for r in built)
        if any('lostParent' in r for r in built):
            card[f'{cfg}.lost_parent_on_gain'] = sum(1 for r in built if r.get('lostParent'))
        for r in named:
            if r['parent'] != r['constructed']:
                detail.append(f"constructed wrong {cfg} {r['gsm']} chr{r['chrom']} built {r['constructed']} named {r['parent']} z {r['z']:.2f}")
    key = lambda r: (r['gsm'], r['chrom'], r['constructed'])
    one = {key(r): r for r in rows if r['config'] == 'one' and r['constructed'] != 'none'}
    two = {key(r): r for r in rows if r['config'] == 'two' and r['constructed'] != 'none'}
    both = [k for k in one if k in two and one[k]['gain'] == 'gain' and two[k]['gain'] == 'gain']
    card['both_named'] = len(both)
    card['both_named_same_parent'] = sum(1 for k in both if one[k]['parent'] == two[k]['parent'])
    return card, detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--invariance', nargs='*', default=[])
    ap.add_argument('--census', nargs='*', default=[])
    ap.add_argument('--constructed', nargs='*', default=[])
    ap.add_argument('--heldout', nargs='*', default=[])
    ap.add_argument('--previous')
    ap.add_argument('--out')
    a = ap.parse_args()
    card, detail = {}, []
    if a.invariance:
        c, d = invariance(load_jsonl(a.invariance)); card['invariance'] = c; detail += d
    if a.census:
        c, d = census(a.census); card['census'] = c; detail += d
    if a.constructed:
        c, d = constructed(load_jsonl(a.constructed)); card['constructed'] = c; detail += d
    if a.heldout:
        c, d = constructed(load_jsonl(a.heldout)); card['constructed_heldout'] = c; detail += d
    print(json.dumps(card, indent=1))
    if detail:
        print('\nDETAIL'); print('\n'.join(detail[:200]))
    if a.previous:
        prev = json.load(open(a.previous))
        moved = [(s, k, prev.get(s, {}).get(k), v) for s, sec in card.items() for k, v in sec.items()
                 if prev.get(s, {}).get(k) != v]
        print('\nMOVED SINCE PREVIOUS' if moved else '\nNOTHING MOVED SINCE PREVIOUS')
        for s, k, o, n in moved:
            print(f'  {s}.{k}: {o} -> {n}')
    if a.out:
        json.dump(card, open(a.out, 'w'), indent=1, sort_keys=True)


if __name__ == '__main__':
    main()
