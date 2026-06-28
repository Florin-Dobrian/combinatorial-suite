#!/usr/bin/env python3
"""
mtx2txt.py - convert a SuiteSparse MatrixMarket (.mtx or .mtx.gz) file to the
plain-text edge-list format used by the matching codes.

  bipartite (default): first line  |S| |T| |E|   then |E| lines of  s t
  general  (--mode general): first line  |V| |E|  then |E| lines of  u v

Indices in the output are 0-based (MatrixMarket is 1-based; we subtract 1).
Rows map to the S side, columns to the T side.

Common cases (the large directed graphs: LAW/*, SNAP/*) are stored 'general'
with no explicit zeros, so the default single-pass fast path applies and the
emitted edge count equals the header nnz.

Usage:
  python3 mtx2txt.py in.mtx[.gz] out.txt
  python3 mtx2txt.py in.mtx.gz out.txt --mode general
  python3 mtx2txt.py in.mtx out.txt --drop-zeros --no-selfloops
"""
import sys, gzip, io, argparse

BATCH = 1 << 20  # edges per write flush


def open_text(path):
    raw = open(path, 'rb')
    if raw.read(2) == b'\x1f\x8b':
        raw.seek(0)
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding='ascii', newline='')
    raw.seek(0)
    return io.TextIOWrapper(raw, encoding='ascii', newline='')


def read_header(fh):
    banner = fh.readline()
    if not banner.startswith('%%MatrixMarket'):
        sys.exit('error: not a MatrixMarket file (missing %%MatrixMarket banner)')
    p = banner.split()
    obj = p[2].lower() if len(p) > 2 else ''
    field = p[3].lower() if len(p) > 3 else 'real'
    sym = p[4].lower() if len(p) > 4 else 'general'
    if obj != 'coordinate':
        sys.exit('error: only sparse "coordinate" MatrixMarket is supported (got %r)' % obj)
    line = fh.readline()
    while line and (line.startswith('%') or not line.strip()):
        line = fh.readline()
    d = line.split()
    return field, sym, int(d[0]), int(d[1]), int(d[2])


def edge_stream(fh, nnz, field, drop_zeros, no_selfloops, mirror):
    """Yield 0-based (a,b) pairs, reading exactly nnz data records."""
    has_val = field in ('real', 'integer', 'complex')
    seen = 0
    for line in fh:
        if seen >= nnz:
            break
        if not line or line[0] == '%':
            continue
        t = line.split()
        if len(t) < 2:
            continue
        seen += 1
        if drop_zeros and has_val:
            if field == 'complex':
                if float(t[2]) == 0.0 and float(t[3]) == 0.0:
                    continue
            elif float(t[2]) == 0.0:
                continue
        a = int(t[0]) - 1
        b = int(t[1]) - 1
        if no_selfloops and a == b:
            continue
        yield a, b
        if mirror and a != b:
            yield b, a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input')
    ap.add_argument('output')
    ap.add_argument('--mode', choices=['bipartite', 'general'], default='bipartite')
    ap.add_argument('--drop-zeros', action='store_true',
                    help='skip entries whose stored value is exactly 0')
    ap.add_argument('--no-selfloops', action='store_true',
                    help='skip i==j entries')
    ap.add_argument('--mirror-symmetric', action='store_true',
                    help='force emitting both (i,j) and (j,i); auto-on for symmetric storage')
    args = ap.parse_args()

    fh = open_text(args.input)
    field, sym, nrows, ncols, nnz = read_header(fh)
    sym_storage = sym in ('symmetric', 'skew-symmetric', 'hermitian')
    mirror = args.mirror_symmetric or sym_storage

    sys.stderr.write('header: %s %s  %d x %d  nnz=%d\n' % (field, sym, nrows, ncols, nnz))
    if sym_storage:
        sys.stderr.write('note: symmetric storage detected; mirroring upper triangle\n')

    fast = not (args.drop_zeros or args.no_selfloops or mirror)

    out = open(args.output, 'w', buffering=1 << 22)

    def write_edges(src):
        buf = []
        n = 0
        ap_ = buf.append
        for a, b in src:
            ap_('%d %d\n' % (a, b))
            n += 1
            if len(buf) >= BATCH:
                out.write(''.join(buf)); buf.clear()
                sys.stderr.write('\r  %d edges' % n); sys.stderr.flush()
        if buf:
            out.write(''.join(buf))
        return n

    if fast:
        # emitted count == nnz, write header immediately, single streaming pass
        header = ('%d %d %d\n' % (nrows, ncols, nnz)) if args.mode == 'bipartite' \
                 else ('%d %d\n' % (nrows, nnz))
        out.write(header)
        n = write_edges(edge_stream(fh, nnz, field, False, False, False))
        if n != nnz:
            sys.stderr.write('\nwarning: emitted %d edges but header said %d\n' % (n, nnz))
    else:
        # filtering/mirroring changes the count: pass 1 counts, pass 2 writes
        n = sum(1 for _ in edge_stream(fh, nnz, field,
                                       args.drop_zeros, args.no_selfloops, mirror))
        fh.close()
        fh = open_text(args.input)
        read_header(fh)
        header = ('%d %d %d\n' % (nrows, ncols, n)) if args.mode == 'bipartite' \
                 else ('%d %d\n' % (nrows, n))
        out.write(header)
        write_edges(edge_stream(fh, nnz, field,
                                args.drop_zeros, args.no_selfloops, mirror))

    out.close()
    sys.stderr.write('\rwrote %d edges to %s\n' % (n, args.output))


if __name__ == '__main__':
    main()
