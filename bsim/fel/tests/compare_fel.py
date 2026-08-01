#!/usr/bin/env python3
"""
Compare the Bmad FEL tracker (fel_ss_test) against Genesis 1.3 Version 4, over the
steady-state Benchmark1-SASE configuration. Run through run_fel_benchmark.sh, which
produces all the inputs; this script only reads and compares.

Three tiers, each with its own tolerance sized to what it measures:

  tier1         One undulator segment, no interludes. Isolates the transcribed FEL core:
                push, deposition, field solve. Observed agreement ~2e-13 on the power
                curve; anything above the tolerance is a transcription defect.

  tier2_genesis The full 6-FODO line with the transcribed Genesis interlude model.
                Everything is transcription, so this should agree at transcription level
                too; the residual (~1e-8 observed) is rounding-difference growth through
                exponential gain, chiefly the interlude phase advance evaluated in a
                different operation order.

  tier2_bmad    The full line with the seam: Bmad tracks the interludes, wavefront_drift
                moves the field, theta advances by the exact mapping from Bmad's z. This
                differs from Genesis by a real transport model difference -- Genesis
                samples the path length term px^2+py^2 once at quad mid-element
                (TrackBeam half step, then the theta step, then the second half), Bmad
                integrates it exactly through the quad map -- measured at ~2e-3 rad of
                bunching phase per quadrupole and ~1.3e-2 of power at saturation. The
                tier exists to demonstrate the seam works; its tolerance brackets the
                measured model difference to catch regressions, and tier2_genesis is
                what prices the difference.

Each tier compares the per-record power and bunching curves, and the final field and
particle dumps element by element. Particle ordering is preserved by both codes (no
sorting happens in steady state), so the dumps compare particle by particle.

A fourth check compares Fortran against itself: tier1 rerun with every particle split
into two coincident copies of weights w/3 and 2w/3 (split_weights = T). Collective
observables are linear in the weights, so the curves and the final field must be
unchanged to round-off. This is the only test of the weighted code paths with nonuniform
weights -- the Genesis dump format carries no weights, so every Genesis comparison sees
the uniform case, where a bug like using one particle's weight for all is invisible.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np


def load_genesis_out(fn):
    with h5py.File(fn) as h5:
        return {
            "z": h5["Lattice/zplot"][:],
            "power": h5["Field/power"][:, 0],
            "bunching": h5["Beam/bunching"][:, 0],
            "bunchingphase": h5["Beam/bunchingphase"][:, 0],
            "energy": h5["Beam/energy"][:, 0],
            "xsize": h5["Beam/xsize"][:, 0],
            "ysize": h5["Beam/ysize"][:, 0],
        }


def load_fortran_diag(fn):
    d = np.loadtxt(fn)
    return {
        "z": d[:, 0], "power": d[:, 1], "bunching": d[:, 3],
        "bunchingphase": d[:, 4], "energy": d[:, 5],
        "xsize": d[:, 7], "ysize": d[:, 8],
    }


def load_fld(fn):
    with h5py.File(fn) as h5:
        n = int(h5["gridpoints"][0])
        return (h5["slice000001/field-real"][:].reshape(n, n)
                + 1j * h5["slice000001/field-imag"][:].reshape(n, n))


def load_par(fn):
    with h5py.File(fn) as h5:
        s = h5["slice000001"]
        return {k: s[k][:] for k in ("gamma", "theta", "x", "y", "px", "py")}


def compare_tier(name, fortran_diag, genesis_out, fortran_fld, genesis_fld,
                 fortran_par, genesis_par, tolerance):
    """
    Compare one tier. Returns (worst_relative_difference, ok).
    """
    print(f"--- {name} " + "-" * (74 - len(name)))

    f = load_fortran_diag(fortran_diag)
    g = load_genesis_out(genesis_out)

    n = min(len(f["z"]), len(g["z"]))
    if len(f["z"]) != len(g["z"]):
        print(f"  FAIL: record counts differ: fortran {len(f['z'])}, genesis {len(g['z'])}")
        return np.inf, False
    if np.abs(f["z"][:n] - g["z"][:n]).max() > 1e-9:
        print(f"  FAIL: record positions differ, max |dz| = "
              f"{np.abs(f['z'][:n]-g['z'][:n]).max():.3e} m")
        return np.inf, False

    worst = 0.0

    # Curves: elementwise relative for power (it spans decades and every point matters),
    # peak normalized for bunching (near zero at the quiet start, where elementwise
    # relative measures nothing but the quiet loading noise floor).
    rel_power = (np.abs(f["power"][:n] - g["power"][:n]) / np.abs(g["power"][:n])).max()
    rel_bunch = np.abs(f["bunching"][:n] - g["bunching"][:n]).max() / np.abs(g["bunching"][:n]).max()
    worst = max(worst, rel_power, rel_bunch)
    print(f"  power curve     ({n} records)   elementwise max rel = {rel_power:.3e}")
    print(f"  bunching curve                  peak normalized     = {rel_bunch:.3e}")

    # Final field dump, peak normalized.
    uf, ug = load_fld(fortran_fld), load_fld(genesis_fld)
    rel_fld = np.abs(uf - ug).max() / np.abs(ug).max()
    worst = max(worst, rel_fld)
    print(f"  final field                     peak normalized     = {rel_fld:.3e}")

    # Final particle dump, particle by particle. gamma relative to itself, transverse
    # peak normalized. These are gated.
    pf, pg = load_par(fortran_par), load_par(genesis_par)
    scales = {"gamma": np.abs(pg["gamma"]).max(),
              "x": np.abs(pg["x"]).max(), "y": np.abs(pg["y"]).max(),
              "px": np.abs(pg["px"]).max(), "py": np.abs(pg["py"]).max()}
    parts = []
    for k in ("gamma", "x", "y", "px", "py"):
        r = np.abs(pf[k] - pg[k]).max() / scales[k]
        worst = max(worst, r)
        parts.append(f"{k} {r:.1e}")
    print(f"  final particles                 per-particle        = " + ", ".join(parts))

    # theta is reported but NOT gated. It is the phase of a particle in its ponderomotive
    # bucket, and at saturation neighboring trajectories separate exponentially, so the
    # worst-particle theta difference measures the Lyapunov amplification of whatever
    # difference exists, not the size of that difference (design brief 9.1: chaotic growth
    # is not a usable comparison metric). The distribution tells the story: the median is
    # the typical particle, the max is the separatrix tail. theta's collective effect IS
    # gated, through the bunching curve and the final field above.
    dth = pf["theta"] - pg["theta"]
    print(f"  final theta (not gated; see comment)  max {np.abs(dth).max():.1e}, "
          f"rms {dth.std():.1e}, median {np.median(np.abs(dth)):.1e} rad")

    ok = worst <= tolerance
    print(f"  LARGEST RELATIVE DIFFERENCE: {worst:.6e}   (tolerance {tolerance:.1e})  "
          + ("ok" if ok else "FAIL"))
    print()
    return worst, ok


def compare_split(name, diag_a, diag_b, fld_a, fld_b, tolerance):
    """
    Fortran vs Fortran: split-weight run against the unsplit run.
    """
    print(f"--- {name} " + "-" * (74 - len(name)))
    a, b = load_fortran_diag(diag_a), load_fortran_diag(diag_b)
    worst = 0.0
    n = min(len(a["z"]), len(b["z"]))
    rel_power = (np.abs(a["power"][:n] - b["power"][:n]) / np.abs(b["power"][:n])).max()
    rel_bunch = np.abs(a["bunching"][:n] - b["bunching"][:n]).max() / np.abs(b["bunching"][:n]).max()
    ua, ub = load_fld(fld_a), load_fld(fld_b)
    rel_fld = np.abs(ua - ub).max() / np.abs(ub).max()
    worst = max(rel_power, rel_bunch, rel_fld)
    print(f"  power curve                     elementwise max rel = {rel_power:.3e}")
    print(f"  bunching curve                  peak normalized     = {rel_bunch:.3e}")
    print(f"  final field                     peak normalized     = {rel_fld:.3e}")
    ok = worst <= tolerance
    print(f"  LARGEST RELATIVE DIFFERENCE: {worst:.6e}   (tolerance {tolerance:.1e})  "
          + ("ok" if ok else "FAIL"))
    print()
    return worst, ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("workdir", help="Directory holding all run outputs")
    p.add_argument("--tol-tier1", type=float, default=1.0e-10)
    p.add_argument("--tol-tier2-genesis", type=float, default=1.0e-6)
    p.add_argument("--tol-tier2-bmad", type=float, default=1.0e-1)
    p.add_argument("--tol-split", type=float, default=1.0e-10)
    args = p.parse_args()
    w = args.workdir

    print("=" * 78)
    print("FEL steady-state benchmark: Bmad tracker against Genesis 1.3 Version 4")
    print("=" * 78)
    print()

    results = []
    for name, diag, out, ffld, gfld, fpar, gpar, tol in (
        ("tier1: FEL core, one undulator segment",
         f"{w}/tier1.diag.txt", f"{w}/Aramis1seg.out.h5",
         f"{w}/tier1-final.fld.h5", f"{w}/Aramis1seg-final.fld.h5",
         f"{w}/tier1-final.par.h5", f"{w}/Aramis1seg-final.par.h5",
         args.tol_tier1),
        ("tier2_genesis: full line, transcribed interludes",
         f"{w}/tier2g.diag.txt", f"{w}/Aramis.out.h5",
         f"{w}/tier2g-final.fld.h5", f"{w}/Aramis-final.fld.h5",
         f"{w}/tier2g-final.par.h5", f"{w}/Aramis-final.par.h5",
         args.tol_tier2_genesis),
        ("tier2_bmad: full line, Bmad seam interludes",
         f"{w}/tier2.diag.txt", f"{w}/Aramis.out.h5",
         f"{w}/tier2-final.fld.h5", f"{w}/Aramis-final.fld.h5",
         f"{w}/tier2-final.par.h5", f"{w}/Aramis-final.par.h5",
         args.tol_tier2_bmad),
    ):
        worst, ok = compare_tier(name, diag, out, ffld, gfld, fpar, gpar, tol)
        results.append((name, worst, ok))

    worst, ok = compare_split(
        "weight_split: nonuniform weights must be invisible",
        f"{w}/tier1s.diag.txt", f"{w}/tier1.diag.txt",
        f"{w}/tier1s-final.fld.h5", f"{w}/tier1-final.fld.h5",
        args.tol_split)
    results.append(("weight_split: nonuniform weights must be invisible", worst, ok))

    print("=" * 78)
    print("Summary")
    print("=" * 78)
    all_ok = True
    for name, worst, ok in results:
        print(f"  {'pass' if ok else 'FAIL'}  {worst:.6e}  {name}")
        all_ok = all_ok and ok
    print()
    print("The tier2_bmad number is a measured transport model difference, not an error:")
    print("Genesis samples the quad path-length term at mid-element, Bmad integrates it")
    print("exactly. tier2_genesis proves it: same code, Genesis's interlude model, and")
    print("the difference collapses by six orders of magnitude. See the README.")
    print()
    print("PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
