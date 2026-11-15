"""What one run.ts record says about each chromosome's copy number.

`intensity` is what the parent-free channels compute: the whole-chromosome call from call rate and
intensity, and the copy segments, printed or withheld. `alleles` adds what needs a parent: a
chromosome the alleles name as missing a parent's copy.
"""
import collections


def chrom_calls(rec, alleles=False, printed_only=False):
    """{chrom: {'loss', 'gain'}} for one record."""
    out = collections.defaultdict(set)
    for c in rec['chroms']:
        if c['aneuploidy']:
            out[c['chrom']].add(c['aneuploidy'])
        if alleles and c['lostParent']:
            out[c['chrom']].add('loss')
    for s in rec['segments']:
        if s['kind'] in ('copy-loss', 'copy-gain'):
            out[s['chrom']].add(s['kind'][5:])
        elif alleles and s['kind'] == 'parental-absence':
            out[s['chrom']].add('loss')
    if not printed_only:
        for x in rec['refusals']:
            w = x.get('withheld')
            if x['claim'] == 'aneuploidy' and w in ('gain', 'loss'):
                out[x['chrom']].add(w)
            elif x['claim'] == 'segment' and isinstance(w, dict) and w['kind'] in ('copy-loss', 'copy-gain'):
                out[x['chrom']].add(w['kind'][5:])
            elif alleles and (x['claim'] == 'lostParent' or (x['claim'] == 'segment' and isinstance(w, dict)
                              and w['kind'] == 'parental-absence')):
                out[x['chrom']].add('loss')
    return {k: frozenset(v) for k, v in out.items()}


def copy_segments(rec, printed_only=False):
    """[(chrom, kind, start, end)] of the copy segments, printed or withheld."""
    segs = [(s['chrom'], s['kind'][5:], s['start'], s['end']) for s in rec['segments']
            if s['kind'] in ('copy-loss', 'copy-gain')]
    if not printed_only:
        for x in rec['refusals']:
            w = x.get('withheld')
            if x['claim'] == 'segment' and isinstance(w, dict) and w['kind'] in ('copy-loss', 'copy-gain'):
                segs.append((w['chrom'], w['kind'][5:], w['start'], w['end']))
    return segs


if __name__ == '__main__':
    r = {'chroms': [{'chrom': '3', 'aneuploidy': 'gain', 'lostParent': ''}, {'chrom': '5', 'aneuploidy': '', 'lostParent': 'maternal'}],
         'segments': [{'chrom': '7', 'kind': 'copy-loss', 'start': 1, 'end': 9}, {'chrom': '8', 'kind': 'absence', 'start': 1, 'end': 2}],
         'refusals': [{'claim': 'aneuploidy', 'chrom': '9', 'withheld': 'loss'}, {'claim': 'segment', 'chrom': '1'},
                      {'claim': 'segment', 'chrom': '2', 'withheld': {'chrom': '2', 'kind': 'copy-gain', 'start': 5, 'end': 6}}]}
    assert chrom_calls(r) == {'3': {'gain'}, '7': {'loss'}, '9': {'loss'}, '2': {'gain'}}
    assert chrom_calls(r, alleles=True)['5'] == {'loss'} and '9' not in chrom_calls(r, printed_only=True)
    r['segments'][1]['kind'] = 'parental-absence'
    assert chrom_calls(r, alleles=True)['8'] == {'loss'} and '8' not in chrom_calls(r)
    assert copy_segments(r) == [('7', 'loss', 1, 9), ('2', 'gain', 5, 6)]
    print('calls ok')
