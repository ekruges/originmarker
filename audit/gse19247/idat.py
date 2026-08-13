"""Minimal Illumina IDAT reader. Enough to get address -> mean intensity, nothing more.

The genotyping IDAT is a flat tagged-field file: a header of (code, offset) pairs, then the
arrays. Only three fields matter here: how many beads were read, their addresses, and their mean
intensities. Everything else in the format is scanner provenance.
"""
import gzip, struct

N_SNPS, ADDRESS, MEAN, BARCODE, CHIP, POSITION = 1000, 102, 104, 402, 403, 404


def _read(f, fmt):
    n = struct.calcsize(fmt)
    return struct.unpack(fmt, f.read(n))[0]


def _string(f):
    # Illumina's 7-bit continuation length prefix, same encoding .NET uses.
    n, shift = 0, 0
    while True:
        b = _read(f, '<B')
        n |= (b & 0x7F) << shift
        if not b & 0x80:
            break
        shift += 7
    return f.read(n).decode('latin-1')


def read_idat(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rb') as f:
        assert f.read(4) == b'IDAT', f'not an IDAT: {path}'
        ver = _read(f, '<q')
        assert ver == 3, f'unsupported IDAT version {ver}'
        fields = {}
        for _ in range(_read(f, '<i')):
            code = _read(f, '<H')
            fields[code] = _read(f, '<q')

        f.seek(fields[N_SNPS])
        n = _read(f, '<i')

        f.seek(fields[ADDRESS])
        addr = struct.unpack(f'<{n}i', f.read(4 * n))

        f.seek(fields[MEAN])
        mean = struct.unpack(f'<{n}H', f.read(2 * n))

        meta = {}
        for name, code in (('barcode', BARCODE), ('chip', CHIP), ('position', POSITION)):
            if code in fields:
                f.seek(fields[code])
                meta[name] = _string(f)
        return addr, mean, meta


if __name__ == '__main__':
    import sys
    a, m, meta = read_idat(sys.argv[1])
    assert len(a) == len(m) and len(a) > 1000, 'address and mean arrays must agree and be non-trivial'
    assert len(set(a)) == len(a), 'addresses must be unique'
    print(f'{len(a)} beads  {meta}  mean intensity min/median/max = '
          f'{min(m)}/{sorted(m)[len(m) // 2]}/{max(m)}')
