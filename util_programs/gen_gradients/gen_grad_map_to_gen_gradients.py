#!/usr/bin/env python3
"""
gen_grad_map_to_gen_gradients.py

Convert an obsolete Bmad `gen_grad_map = {...}` block (the straight-frame,
azimuthal-harmonic C_{m,sin/cos} generalized-gradient format) to the new
curved-coordinate `gen_gradients = {...}` block (midplane a_n / b_n / b_s
functions), used by the formalism from GeneralizedGradients.jl (issue #2086).

Usage:
    gen_grad_map_to_gen_gradients.py <old_gg_map_block> [--g-ref <value>] [-o <out>]

`<old_gg_map_block>` is a file containing a single `{ ... }` gen_grad_map body
(as written by Bmad's write_bmad_lattice_file, or GeneralizedGradients.jl's
write_bmad_gg_fit into `<base>_gg.bmad`). Writes the equivalent gen_gradients
block to stdout (or `-o`).

Mapping (inverse of gg_to_bmad_curves in GeneralizedGradients.jl):
    b_m^{[j]} = m! C^{[j]}_{m,s} + (m-1)! Σ_{n≥1, m-2n≥1} Wn(m,n) C^{[j+2n]}_{m-2n,s}
    a_m^{[j]} = m! C^{[j]}_{m,c} + (m-1)! Σ_{n≥1, m-2n≥1} Wc(m,n) C^{[j+2n]}_{m-2n,c}
                                 + [m even] (m-1)! Us(m) b_s^{[m+j-1]}
    b_s^{[i]} = C^{[i+1]}_{0,c}
with Wn(m,n) = (-1)^n (m-2n)!(m-2n)/(4^n n!(m-n)!),
     Wc(m,n) = (-1)^n (m-2n)! m   /(4^n n!(m-n)!),
     Us(m)   = (-1)^{m/2} m/(4^{m/2} ((m/2)!)^2).

Curve kinds: old (m, sin) -> new (b, n=m) [normal]; old (m, cos) -> new (a, n=m)
[skew]; old (m=0, cos) -> new (bs, n=0) [solenoid].

The midplane gradients a_n, b_n, b_s are curvature-independent, so the field
data carries over unchanged; `g_ref` (the new reference-frame curvature, 0 for a
straight map) is supplied with --g-ref (default 0). For a straight map the new
block reproduces the old field exactly; for a curved map it reproduces the
*correct* curved field (the old straight-frame expansion was only approximate).
"""

import sys
import re
import argparse
from math import factorial


def parse_block(text):
    """Parse a gen_grad_map { ... } body. Returns dict with top-level scalars and
    a list of curves: {m, kind, dz, planes:[z...], derivs:[[vals per order]...]}."""
    top = {}
    # r0 is a parenthesized list; capture the whole (...) first.
    m = re.search(r'\br0\s*=\s*(\([^)]*\))', text)
    if m:
        top['r0'] = m.group(1).strip()
    # Other top-level scalars (stop at comma/brace/newline).
    for key in ('field_type', 'ele_anchor_pt', 'dz', 'field_scale',
                'master_parameter', 'curved_ref_frame'):
        m = re.search(r'\b' + key + r'\s*=\s*([^,\}\n]+)', text)
        if m:
            top[key] = m.group(1).strip()

    curves = []
    # Each curve = { m = .., kind = .., derivs = { z: v v v, ... } }
    for cm in re.finditer(r'curve\s*=\s*\{(.*?)\}\s*\}', text, re.S):
        body = cm.group(1)
        mm = re.search(r'\bm\s*=\s*(-?\d+)', body)
        km = re.search(r'\bkind\s*=\s*(\w+)', body)
        dm = re.search(r'derivs\s*=\s*\{(.*)$', body, re.S)
        if not (mm and km and dm):
            continue
        planes, derivs = [], []
        for line in dm.group(1).splitlines():
            line = line.strip().rstrip(',')
            if ':' not in line:
                continue
            zstr, vstr = line.split(':', 1)
            try:
                z = float(zstr.strip())
            except ValueError:
                continue
            vals = [float(x) for x in vstr.split()]
            planes.append(z)
            derivs.append(vals)
        curves.append(dict(m=int(mm.group(1)), kind=km.group(1).lower().strip(),
                           planes=planes, derivs=derivs))
    return top, curves


def build_tower(curve, npl, nder):
    """Return C[j][i] for j=0..nder, i=0..npl-1 from a parsed curve (0 if absent)."""
    T = [[0.0] * npl for _ in range(nder + 1)]
    for i, row in enumerate(curve['derivs']):
        if i >= npl:
            break
        for j, v in enumerate(row):
            if j <= nder:
                T[j][i] = v
    return T


def convert(text, g_ref):
    top, curves = parse_block(text)
    if not curves:
        sys.exit('error: no curves found in gen_grad_map block')

    npl = max(len(c['planes']) for c in curves)
    planes = next(c['planes'] for c in curves if len(c['planes']) == npl)
    dz = float(top.get('dz', (planes[1] - planes[0]) if npl > 1 else 0.0))

    # Sort old curves into normal (sin) cs[m], skew (cos) cc[m], solenoid c0c.
    cs, cc, c0c = {}, {}, None
    mmax = 0
    for c in curves:
        nder = max((len(r) for r in c['derivs']), default=1) - 1
        T = build_tower(c, npl, nder)
        if c['m'] == 0 and c['kind'] == 'cos':
            c0c = T                                  # C^{[j]}_{0,cos}, j = 0..nder (nder = m_max+1)
        elif c['kind'] == 'sin':
            cs[c['m']] = T; mmax = max(mmax, nder)
        elif c['kind'] == 'cos':
            cc[c['m']] = T; mmax = max(mmax, nder)

    kmax = max([m for m in list(cs) + list(cc)], default=0)

    def Wn(k, n): return (-1.0) ** n * factorial(k - 2 * n) * (k - 2 * n) / (4.0 ** n * factorial(n) * factorial(k - n))
    def Wc(k, n): return (-1.0) ** n * factorial(k - 2 * n) * k / (4.0 ** n * factorial(n) * factorial(k - n))
    def Us(k):    return (-1.0) ** (k // 2) * k / (4.0 ** (k // 2) * factorial(k // 2) ** 2)

    def bs(order, i):   # b_s^{[order]} = C^{[order+1]}_{0,cos}
        if c0c is None or order + 1 >= len(c0c):
            return 0.0
        return c0c[order + 1][i]

    # New towers.
    b_new, a_new = {}, {}
    for k in range(1, kmax + 1):
        if k in cs:
            Tb = [[0.0] * npl for _ in range(mmax + 1)]
            for j in range(mmax + 1):
                for i in range(npl):
                    val = factorial(k) * cs[k][j][i]
                    n = 1
                    while k - 2 * n >= 1:
                        lo = cs.get(k - 2 * n)
                        if lo is not None and j + 2 * n < len(lo):
                            val += factorial(k - 1) * Wn(k, n) * lo[j + 2 * n][i]
                        n += 1
                    Tb[j][i] = val
            b_new[k] = Tb
        if k in cc:
            Ta = [[0.0] * npl for _ in range(mmax + 1)]
            for j in range(mmax + 1):
                for i in range(npl):
                    val = factorial(k) * cc[k][j][i]
                    n = 1
                    while k - 2 * n >= 1:
                        lo = cc.get(k - 2 * n)
                        if lo is not None and j + 2 * n < len(lo):
                            val += factorial(k - 1) * Wc(k, n) * lo[j + 2 * n][i]
                        n += 1
                    if k % 2 == 0:
                        val += factorial(k - 1) * Us(k) * bs(k + j - 1, i)
                    Ta[j][i] = val
            a_new[k] = Ta

    bs_new = None
    if c0c is not None:
        bs_new = [[bs(j, i) for i in range(npl)] for j in range(mmax + 1)]

    return top, planes, dz, a_new, b_new, bs_new, mmax, g_ref


def fmt(v):
    return f'{v:.12g}'


def emit(top, planes, dz, a_new, b_new, bs_new, g_ref, out):
    npl = len(planes)
    p = out.write
    p('{\n')
    p('  field_type    = %s,\n' % top.get('field_type', 'magnetic'))
    p('  ele_anchor_pt = %s,\n' % top.get('ele_anchor_pt', 'beginning'))
    if 'master_parameter' in top:
        p('  master_parameter = %s,\n' % top['master_parameter'])
    if 'field_scale' in top:
        p('  field_scale   = %s,\n' % top['field_scale'])
    p('  r0            = %s,\n' % top.get('r0', '(0, 0, 0)'))
    p('  dz            = %s,\n' % fmt(dz))
    p('  g_ref         = %s,\n' % fmt(g_ref))

    # Assemble the curve list, then emit with a separating comma after all but the last.
    clist = []
    if bs_new is not None:
        clist.append(('bs', 0, bs_new))
    clist += [('b', n, b_new[n]) for n in sorted(b_new)]
    clist += [('a', n, a_new[n]) for n in sorted(a_new)]

    for ci, (kind, n, tower) in enumerate(clist):
        nder = len(tower) - 1
        p('  curve = {\n')
        p('    kind = %s,\n' % kind)
        p('    n = %d,\n' % n)
        p('    derivs = {\n')
        for i in range(npl):
            vals = ' '.join(fmt(tower[j][i]) for j in range(nder + 1))
            comma = ',' if i < npl - 1 else ''
            p('      %s: %s%s\n' % (fmt(planes[i]), vals, comma))
        p('    }\n')
        p('  }%s\n' % (',' if ci < len(clist) - 1 else ''))
    p('}\n')


def main():
    ap = argparse.ArgumentParser(description='Convert a gen_grad_map block to gen_gradients.')
    ap.add_argument('infile')
    ap.add_argument('--g-ref', type=float, default=0.0,
                    help='reference-frame curvature 1/rho for the new map (default 0)')
    ap.add_argument('-o', '--out', default=None)
    args = ap.parse_args()

    with open(args.infile) as f:
        text = f.read()
    top, planes, dz, a_new, b_new, bs_new, mmax, g_ref = convert(text, args.g_ref)
    out = open(args.out, 'w') if args.out else sys.stdout
    emit(top, planes, dz, a_new, b_new, bs_new, g_ref, out)
    if args.out:
        out.close()


if __name__ == '__main__':
    main()
