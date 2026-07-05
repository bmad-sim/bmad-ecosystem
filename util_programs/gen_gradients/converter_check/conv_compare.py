import re, sys
newblock, origfile = sys.argv[1], sys.argv[2]
# Parse new block: for each curve (kind,n) get derivs at z=0.1 (2nd plane).
text = open(newblock).read()
got = {}
for cm in re.finditer(r'curve\s*=\s*\{(.*?)\}\s*\}', text, re.S):
    body = cm.group(1)
    kind = re.search(r'kind\s*=\s*(\w+)', body).group(1)
    n = int(re.search(r'\bn\s*=\s*(-?\d+)', body).group(1))
    dm = re.search(r'derivs\s*=\s*\{(.*)$', body, re.S).group(1)
    for line in dm.splitlines():
        line = line.strip().rstrip(',')
        if ':' not in line: continue
        z, vstr = line.split(':', 1)
        if abs(float(z) - 0.1) < 1e-9:
            for j, v in enumerate(float(x) for x in vstr.split()):
                got[(kind, n, j)] = v
# Original values.
maxerr = 0.0
for line in open(origfile):
    kind, n, j, v = line.split()
    key = (kind, int(n), int(j)); v = float(v)
    g = got.get(key)
    if g is None:
        print("MISSING", key); continue
    e = abs(g - v) / max(abs(v), 1e-12)
    maxerr = max(maxerr, e)
    if e > 1e-6:
        print(f"MISMATCH {key}: got {g:.10e} orig {v:.10e} relerr {e:.2e}")
print(f"max rel err (a/b/bs recovery) = {maxerr:.3e}")
print("CONVERTER RECOVERY PASSED" if maxerr < 1e-6 else "CONVERTER RECOVERY FAILED")
