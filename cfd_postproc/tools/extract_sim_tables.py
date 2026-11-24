"""Extract (time, value) tables embedded in a STAR-CCM+ .sim file, without STAR.

The .sim stores its object tree as text dicts. A FileTable's columns live in binary
GlobalArray<Float8> blocks: descriptor {'ClassName': 'GlobalArray<Float8>', 'dataKey': K, 'size': N}
points at file offset K, followed by a 76-byte header and N little-endian doubles.
Each table is stored as [value column, time column]: two consecutive same-size arrays whose
SECOND one is a monotonic time axis. Tables appear in the same order as the FileTable list. Usage: python3 extract_sim_tables.py <file.sim> <outdir>
"""
import mmap, os, re, sys
import numpy as np

HDR = 76
sim, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
f = open(sim, 'rb'); m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
arrs = [(int(a.group(1)), int(a.group(2)), int(a.group(3))) for a in
        re.finditer(rb"'ClassName': 'GlobalArray<Float8>', 'ObjectId': (\d+)L, 'dataKey': (\d+)L, 'size': (\d+)\}", m)]
seen = set(); tables = []
for oid, key, n in arrs:
    if key in seen or n < 20 or n > 100000:
        continue
    seen.add(key)
    v = np.frombuffer(m[key + HDR: key + HDR + 8 * n], dtype='<f8')
    tables.append((key, n, v))
tables.sort()
# table names in file order, to label the pairs
names = re.findall(rb"'ClassName': 'star.common.FileTable', [^}]*?'FileName': '([^']*)'[^}]*?'PresentationName': '([^']*)'", m)
print('FileTables in sim:'); [print('   %-28s <- %s' % (p.decode(), fn.decode())) for fn, p in names if b'Tables\\\\table_' not in fn and b'Tables/table_' not in fn]
k = 0; found = 0
while k < len(tables) - 1:
    (k1, n1, t), (k2, n2, y) = tables[k], tables[k + 1]
    y, t = t, y                                    # stored as [values, time]
    if n1 == n2 and np.all(np.diff(t) > 0) and 0.5 < t[-1] < 10 and np.all(np.isfinite(y)):
        out = os.path.join(outdir, 'table_%d_n%d.csv' % (found, n1))
        np.savetxt(out, np.c_[t, y], delimiter=',', header='"time (s)","value"', comments='', fmt='%.17g')
        print('pair: n=%d  t %.4f..%.4f s  value %.4g..%.4g   -> %s' % (n1, t[0], t[-1], y.min(), y.max(), out))
        found += 1; k += 2
    else:
        k += 1
print('%d time-tables found' % found)
