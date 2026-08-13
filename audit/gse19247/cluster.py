"""Estimate per-marker genotype cluster positions for GPL6985 from the series itself.

WHY. Illumina's BAF and genotype calls come from a cluster file (.egt) holding, for every marker,
where the AA / AB / BB clouds sit in theta. Without it a global threshold has to stand in, and
that fails measurably: on this series it called only 31-43% of markers homozygous in bulk DNA
against an expected 67-70%, and could not tell a father from a mother using the father's own
sperm. The cluster file is not public, but it does not have to be: cluster positions are a
property of the CHIP, and 255 arrays on this chip are in this same series. That is enough to
estimate them.

Allele dropout in single-cell MDA moves a cell BETWEEN clusters; it does not move where the
clusters ARE. So a series that is mostly single cells still locates the clouds, which is the only
thing being estimated here.

Output is an .npz holding, per marker, the three cluster means plus the address and locus.
"""
import gzip, json, os, struct, sys
import numpy as np

MIN_R_FRAC = 0.30      # per-sample intensity floor, as a fraction of that sample's own median
AB_MIN_FRAC = 0.02     # below this share of samples the AB cloud is not located, only inferred
ITERS = 6


def read_idat_np(path):
    with gzip.open(path, 'rb') as f:
        buf = f.read()
    assert buf[:4] == b'IDAT', path
    (ver,) = struct.unpack_from('<q', buf, 4)
    assert ver == 3, ver
    (nf,) = struct.unpack_from('<i', buf, 12)
    fields = {}
    for i in range(nf):
        code, off = struct.unpack_from('<Hq', buf, 16 + i * 10)
        fields[code] = off
    (n,) = struct.unpack_from('<i', buf, fields[1000])
    addr = np.frombuffer(buf, '<i4', n, fields[102])
    mean = np.frombuffer(buf, '<u2', n, fields[104])
    return addr, mean


def manifest(path='gpl6985.txt'):
    addr, chrom, pos = [], [], []
    with open(path) as f:
        hdr = None
        for line in f:
            if line.startswith('!platform_table_begin'):
                cols = next(f).rstrip('\n').split('\t')
                hdr = {n: k for k, n in enumerate(cols)}
                continue
            if hdr is None or line.startswith('!'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) <= hdr['CNV_Probe'] or p[hdr['AddressB_ID']].strip():
                continue
            a, c, m = p[hdr['AddressA_ID']].strip(), p[hdr['Chr']].strip(), p[hdr['MapInfo']].strip()
            if not (a and c and m) or p[hdr['Intensity_Only']].strip() == '1':
                continue
            addr.append(int(a)); chrom.append(c); pos.append(int(m))
    return np.array(addr, '<i4'), np.array(chrom), np.array(pos, '<i8')


def _gather(addrs, means, want):
    """means, reordered onto `want`. Zero where the chip has no bead for that address."""
    order = np.argsort(addrs)
    sorted_addr = addrs[order]
    at = np.clip(np.searchsorted(sorted_addr, want), 0, len(sorted_addr) - 1)
    return np.where(sorted_addr[at] == want, means[order][at], 0).astype(np.float32)


BG_PCT, SCALE_PCT = 2, 98


def _normalise(x):
    bg, top = np.percentile(x, BG_PCT), np.percentile(x, SCALE_PCT)
    return np.clip(x - bg, 0, None) / max(top - bg, 1e-6)


def sample_theta(grn, red, addr_idx):
    """theta and total intensity for one sample, aligned to the manifest address order.

    Each channel is aligned through the manifest explicitly rather than assuming the two IDATs
    are written in the same bead order.
    """
    A = _gather(*read_idat_np(grn), addr_idx)
    B = _gather(*read_idat_np(red), addr_idx)
    # Normalise the two channels against each other before taking any angle. The green and red
    # dyes differ in brightness and background, and by a different amount on every array: on one
    # bulk sample here the 98th percentiles are 32,170 and 20,813. Left raw, that skew puts the
    # mode of theta in the MIDDLE of the range, which reads as a 50% heterozygous sample and is
    # impossible in a real person. Background from a low percentile, scale from a high one, both
    # per channel and per sample, which is the affine part of what a cluster file's normalisation
    # does and the part that matters here.
    A, B = _normalise(A), _normalise(B)
    tot = A + B
    theta = (2.0 / np.pi) * np.arctan2(B, np.maximum(A, 1e-6))
    return theta.astype(np.float32), tot


def main():
    addr, chrom, pos = manifest()
    print(f'{len(addr)} Infinium II markers')
    samples = json.load(open('gpl6985_samples.json'))
    usable = []
    for s in samples:
        fs = ['idat/' + f for f in s['idats']]
        g = [f for f in fs if 'Grn' in f]
        r = [f for f in fs if 'Red' in f]
        if g and r and os.path.exists(g[0]) and os.path.exists(r[0]):
            usable.append((s, g[0], r[0]))
    print(f'{len(usable)} samples with both channels on disk')

    T = np.zeros((len(usable), len(addr)), np.float32)
    ok = np.zeros((len(usable), len(addr)), bool)
    for i, (s, g, r) in enumerate(usable):
        th, tot = sample_theta(g, r, addr)
        med = np.median(tot[tot > 0]) if (tot > 0).any() else 1.0
        T[i] = th
        ok[i] = tot >= MIN_R_FRAC * med
        if i % 40 == 0:
            print(f'  read {i}/{len(usable)}', flush=True)

    # 1-D three-means per marker, initialised at the canonical positions. Vectorised across all
    # markers at once: each iteration is three masked means over the sample axis.
    # Chunked over markers: the distance array is samples x chunk x 3, which would be a gigabyte
    # if built over all 343k markers at once.
    mu = np.tile(np.array([0.10, 0.50, 0.90], np.float32), (len(addr), 1))
    counts = np.zeros((len(addr), 3), np.int32)
    CHUNK = 20_000
    for lo in range(0, len(addr), CHUNK):
        hi = min(lo + CHUNK, len(addr))
        Tc, okc, muc = T[:, lo:hi], ok[:, lo:hi], mu[lo:hi]
        for _ in range(ITERS):
            d = np.abs(Tc[:, :, None] - muc[None, :, :])
            d[~okc] = np.inf
            lab = np.argmin(d, axis=2)
            lab[~okc] = -1
            for k in range(3):
                m = (lab == k)
                cnt = m.sum(0)
                hit = cnt > 0
                muc[hit, k] = np.where(m, Tc, 0).sum(0)[hit] / cnt[hit]
        mu[lo:hi] = muc
        counts[lo:hi] = np.stack([(lab == k).sum(0) for k in range(3)], 1)

    # Where the AB cloud is too thin to be located, it is inferred as the midpoint of the two
    # homozygote clouds rather than left at wherever the iteration drifted. This is the same
    # fallback a cluster file uses for a marker with no observed heterozygote.
    n_ok = ok.sum(0).astype(np.float32)
    thin = counts[:, 1] < np.maximum(AB_MIN_FRAC * n_ok, 3)
    mu[thin, 1] = 0.5 * (mu[thin, 0] + mu[thin, 2])
    # Clusters must stay ordered, or BAF interpolation is meaningless.
    mu = np.sort(mu, axis=1)
    print(f'AB cloud inferred rather than observed at {thin.sum()} of {len(addr)} markers')

    np.savez_compressed('clusters_gpl6985.npz', addr=addr, chrom=chrom, pos=pos,
                        mu=mu, counts=counts, n_ok=n_ok,
                        gsms=np.array([s['gsm'] for s, _, _ in usable]))
    print('wrote clusters_gpl6985.npz')


if __name__ == '__main__':
    main()
