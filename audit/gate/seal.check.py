"""The open list refuses sealed, exposed and unknown arrays and allows open and laboratory ones.

Run: python3 audit/gate/seal.check.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import seal  # noqa: E402

split = {'11111-01': 'open', '11111-02': 'sealed', '11111-03': 'exposed', 'GSM1': 'open', 'GSM2': 'sealed'}
paths = {
    '/d/GSE288845/11111-01_9.probes.gz': 'open',
    '/d/GSE288845/11111-02_9.probes.gz': 'refused',
    '/d/DIETER/11111-03_9.CEL.probes': 'refused',
    '/d/DIETER/22222-01_9.CEL.probes': 'open',
    '/d/GSE315671/22222-02_9.probes.gz': 'refused',
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
