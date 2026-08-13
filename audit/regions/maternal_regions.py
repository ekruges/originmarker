"""Parent-labelled REGIONS on experiment 14801, against a confirmed maternal reference.

The reference is 14801-05, the single cumulus sample in the whole set: maternal somatic DNA, bulk
quality on this platform (call 0.950, het 0.145 against 0.954 / 0.168 for known bulk gDNA). Every
other usable array in the experiment reads as her material at 0.0012 to 0.114 opposite-homozygote
rate, against ~0.30 for unrelated adults, so the assignment is measured rather than assumed.

THE STATISTIC. At a marker where the MOTHER is homozygous, whatever else the cell carries came
from the father. Orient every marker to MATERNAL SHARE, the fraction of the cell's signal that
matches her allele:

    mother AA -> share = 1 - baf        mother BB -> share = baf

    both parents present      the cell is AB wherever the father differs, so share sits mid-range
    PATERNAL copy missing     only her allele remains, share drives toward 1
    MATERNAL copy missing     only his allele remains, share drives toward 0

Direction is therefore read off the sign, and it is the opposite of the paternal-share convention
used when the father is the reference. Centring is on the SAMPLE'S OWN median, never a theoretical
value, because amplification and marker ascertainment move each array's baseline.
"""
import glob, sys
import numpy as np
from regions import read_probes

USB = "/Volumes/SANDISK USB/SNP array data"
AUTOSOMES = [str(i) for i in range(1, 23)]

WIN = 300          # markers per window
STEP = 150
MARGIN = 0.12      # distance from the sample's own centre before a window counts
MIN_WIN = 2        # consecutive windows, so one noisy window is not a region
MAX_OPP = 0.15     # above this the array is not this mother's material, or has failed
MIN_MARKERS = 20_000


def main():
    mumf = glob.glob(f'{USB}/ROBLES/**/14801-05_*.probes', recursive=True)[0]
    mum = read_probes(mumf)
    hom = {k: v for k, v in mum.items() if v[3] in ('0', '2') and v[0] in AUTOSOMES}
    print(f'maternal reference 14801-05: {len(hom)} autosomal homozygous markers\n')

    cells = sorted(f for f in glob.glob(f'{USB}/ROBLES/**/14801-*.probes', recursive=True)
                   if '14801-05_' not in f)
    print(f'{"cell":<11} {"chr":>3} {"start Mb":>9} {"end Mb":>9} {"span Mb":>8} {"markers":>8} '
          f'{"share":>6} {"dev":>7}  parental call')
    print('-' * 104)
    total = 0
    scanned = 0
    for cf in cells:
        cell = cf.split('/')[-1].split('_')[0]
        s = read_probes(cf)
        opp = tot = 0
        for k, v in hom.items():
            w = s.get(k)
            if w and w[3] in ('0', '2'):
                tot += 1
                if w[3] != v[3]:
                    opp += 1
        if not tot or opp / tot > MAX_OPP:
            continue
        rows = []
        for k, (c, pos, _, gm) in hom.items():
            w = s.get(k)
            if not w or w[3] == '-1' or not np.isfinite(w[2]):
                continue
            rows.append((c, pos, (1.0 - w[2]) if gm == '0' else w[2]))
        if len(rows) < MIN_MARKERS:
            continue
        scanned += 1
        centre = float(np.median([r[2] for r in rows]))
        by = {}
        for c, pos, sh in rows:
            by.setdefault(c, []).append((pos, sh))
        for c in AUTOSOMES:
            arr = sorted(by.get(c, []))
            if len(arr) < WIN:
                continue
            pos = np.array([a[0] for a in arr], dtype=np.int64)
            sh = np.array([a[1] for a in arr])
            flagged = []
            for i in range(0, len(arr) - WIN + 1, STEP):
                if abs(float(np.median(sh[i:i + WIN])) - centre) >= MARGIN:
                    flagged.append(i)
            # merge consecutive flagged windows into one region
            runs, cur = [], []
            for i in flagged:
                if cur and i - cur[-1] > STEP:
                    runs.append(cur); cur = []
                cur.append(i)
            if cur:
                runs.append(cur)
            for run in runs:
                if len(run) < MIN_WIN:
                    continue
                i0, i1 = run[0], min(run[-1] + WIN, len(sh) - 1)
                seg = sh[i0:i1]
                dev = float(np.median(seg)) - centre
                call = ('PATERNAL copy missing' if dev > 0 else 'MATERNAL copy missing')
                print(f'{cell:<11} {c:>3} {pos[i0]/1e6:>9.2f} {pos[i1]/1e6:>9.2f} '
                      f'{(pos[i1]-pos[i0])/1e6:>8.1f} {len(seg):>8} '
                      f'{float(np.median(seg)):>6.3f} {dev:>+7.3f}  {call}')
                total += 1
    print('-' * 104)
    print(f'{total} regions across {scanned} cells scored against the maternal reference')


if __name__ == '__main__':
    main()
