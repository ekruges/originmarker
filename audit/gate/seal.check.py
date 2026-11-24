"""The open list refuses sealed, exposed and unknown arrays and allows open and laboratory ones.

Run: python3 audit/gate/seal.check.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import seal  # noqa: E402

# Barcodes are made here, not written, so no tracked file holds one (release-check.sh).
bc = lambda d, n: f'{d * 5}-{n:02d}'  # noqa: E731
split = {bc('1', 1): 'open', bc('1', 2): 'sealed', bc('1', 3): 'exposed', 'GSM1': 'open', 'GSM2': 'sealed'}
paths = {
    f"/d/GSE288845/{bc('1', 1)}_9.probes.gz": 'open',
    f"/d/GSE288845/{bc('1', 2)}_9.probes.gz": 'refused',
    f"/d/DIETER/{bc('1', 3)}_9.CEL.probes": 'refused',
    f"/d/DIETER/{bc('2', 1)}_9.CEL.probes": 'open',
    f"/d/GSE315671/{bc('2', 2)}_9.probes.gz": 'refused',
    '/d/trios/GSM1.probes.gz': 'open',
    '/d/trios/GSM2.probes.gz': 'refused',
    '/d/trios/GSM3.probes.gz': 'refused',
    '/d/GSE288845/A1-1_9.CEL.probes': 'refused',
    '/r/rebuilt-father.probes': 'open',
}
got = seal.classify(list(paths), split, {'/r/rebuilt-father.probes'})
for p, want in paths.items():
    assert got[p][1] == want, (p, got[p], want)
print('seal checks pass')
