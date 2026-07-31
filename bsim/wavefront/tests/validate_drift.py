#!/usr/bin/env python3
"""
Python half of the Fortran against Python free-space drift comparison.

Run through run_validation.sh, which builds the Fortran side, sets PYTHONPATH to an
openPMD-beamphysics checkout and calls this script twice:

    validate_drift.py write   <input.h5>  --nx ... --nz ... --wavelength ... --dx ... --dz ...
    validate_drift.py compare <input.h5> <fortran_out.h5> --z ...

The first call builds a test wavefront and writes it in Genesis4 format. The second reads
the same input, drifts it with openPMD-beamphysics' drift_wavefront, reads the file the
Fortran program produced, and reports the largest relative difference between the two.

The input is deliberately not a centered symmetric Gaussian. A symmetric field hides a
whole class of mistakes: a transposed transverse grid, x and y swapped somewhere in the
HDF5 index ordering, or a sign error in the kernel's negative wavenumbers all leave a
symmetric field unchanged. The test field is therefore offset from the axis, elliptical,
tilted, and different in every slice.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile

import h5py
import numpy as np

from beamphysics.wavefront import Wavefront
from beamphysics.wavefront.propagators import drift_wavefront


def build_wavefront(nx: int, ny: int, nz: int, dx: float, dy: float, dz: float,
                    wavelength: float) -> Wavefront:
    """
    Build the test wavefront: an off-axis tilted elliptical Gaussian whose centroid,
    size and phase all vary from slice to slice.

    Asymmetry is the point. Every deliberate feature below breaks a symmetry that would
    otherwise let a bug through:

      - x0 != y0 and sigma_x != sigma_y  break the x/y interchange symmetry, so a
        transposed slice or a swapped dx/dy shows up.
      - the xy cross term tilts the ellipse, so a transpose is visible even if the two
        sizes were equal.
      - the linear phase gives the field a nonzero mean transverse wavenumber, so the
        centroid moves during the drift and the sign of the kernel matters.
      - the quadratic phase makes the field converge in x and diverge in y, so the two
        transverse axes evolve differently.
      - varying everything with slice index means a mistake in the slice axis, or a
        slice written into the wrong group, cannot cancel.
    """
    x = (np.arange(nx) - (nx - 1) / 2) * dx
    y = (np.arange(ny) - (ny - 1) / 2) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")

    k0 = 2 * np.pi / wavelength
    w = 8 * dx  # Reference transverse scale, comfortably resolved and well inside the grid.

    Ex = np.empty((nx, ny, nz), dtype=np.complex128)

    for iz in range(nz):
        f = iz / max(nz - 1, 1)  # 0 .. 1 across the pulse.

        x0 = (0.9 + 0.4 * f) * w  # Off axis, and moving with slice.
        y0 = -1.3 * w
        sx = (1.0 + 0.5 * f) * w
        sy = 1.7 * w  # Elliptical.

        u = (X - x0) / sx
        v = (Y - y0) / sy

        # Amplitude: tilted ellipse (the u*v term), amplitude varying along the pulse.
        amp = (0.4 + f) * np.exp(-(u**2 + v**2 + 0.6 * u * v) / 2)

        # Phase: a linear part that steers the beam, and a quadratic part with opposite
        # signs in x and y so the two axes focus and defocus respectively.
        tilt = k0 * (2.5e-5 * (X - x0) - 1.1e-5 * (Y - y0))
        curv = k0 * ((X - x0) ** 2 / (2 * 4.0) - (Y - y0) ** 2 / (2 * 7.0))

        Ex[:, :, iz] = amp * np.exp(1j * (tilt + curv + 0.37 * iz))

    # Scale to a physically plausible field amplitude in V/m.
    Ex *= 1e9 / np.abs(Ex).max()

    return Wavefront(Ex=Ex, dx=dx, dy=dy, dz=dz, wavelength=wavelength)


def file_layout(path: str) -> dict[str, tuple]:
    """
    Walk an HDF5 file and return {name: (kind, shape)} for every object in it, where kind is
    'group' or the NumPy dtype kind character of a dataset.

    The dtype is recorded in full so that width differences can be reported, but see
    compare_structure for which differences are treated as failures.
    """
    layout: dict[str, tuple] = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Group):
            layout[name] = ("group", None)
        else:
            layout[name] = (str(obj.dtype), obj.shape)

    with h5py.File(path, "r") as h5:
        h5.visititems(visit)

    return layout


# Scalar metadata whose numeric type may differ between the Fortran writer and
# Wavefront.write_genesis4 without that being an incompatibility. In each case the
# difference comes from NumPy inferring a type from a Python int, every reader involved
# converts on read, and the Fortran choice is the one Genesis itself writes.
#
#   gridpoints, slicecount  Genesis writes these through writeSingleNodeInt, which takes a
#                           vector<int>, so 32 bit. np.asarray([nx]) gives 64 bit. Both are
#                           read by H5LTread_dataset into a native int and by h5py directly.
#
#   refposition             Genesis writes this as a double with a unit attribute of "m",
#                           and reads it with H5Dread into H5T_NATIVE_DOUBLE. The Python
#                           writer's signature is refposition: float = 0, and np.asarray([0])
#                           on that default Python int gives int64 -- pass an actual float
#                           and it writes float64. So the Python default writes a position in
#                           meters as an integer. Harmless, since HDF5 converts on read, but
#                           it is the Python side that departs from the format here, not the
#                           Fortran side, and Fortran writes float64 deliberately.
_TYPE_DIFF_ALLOWED = {
    "gridpoints": "Genesis writes int32 here; NumPy infers int64",
    "slicecount": "Genesis writes int32 here; NumPy infers int64",
    "refposition": "Genesis writes float64 here; the Python writer's default arg is a "
                   "Python int, so NumPy infers int64",
}


def compare_structure(fortran_file: str, reference_file: str) -> bool:
    """
    Require the Fortran output file to have the same object names and shapes as a file
    written by Wavefront.write_genesis4, and the same dtypes except where
    _TYPE_DIFF_ALLOWED says otherwise.

    This is not redundant with the field comparison. Wavefront.from_genesis4 collects the
    slice groups by globbing 'slice*' and sorting the names, so it reads a file whose groups
    are numbered from any starting index and gets the right answer. Genesis does not: it
    builds the name with sprintf(name, "slice%6.6d/field-real", islice) for islice in
    1..slicecount, so a file numbered from zero, or from two, is silently empty to it. An
    off-by-one in the slice group names is therefore invisible to the field comparison and
    fatal in use, which is exactly the combination worth an explicit check.
    """
    got = file_layout(fortran_file)
    want = file_layout(reference_file)

    ok = True

    missing = sorted(set(want) - set(got))
    extra = sorted(set(got) - set(want))

    if missing:
        print(f"  FAIL: objects the reference has and the Fortran output does not: "
              f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}")
        ok = False
    if extra:
        print(f"  FAIL: objects the Fortran output has and the reference does not: "
              f"{', '.join(extra[:8])}{' ...' if len(extra) > 8 else ''}")
        ok = False

    noted = []
    for name in sorted(set(got) & set(want)):
        g_type, g_shape = got[name]
        w_type, w_shape = want[name]

        if g_shape != w_shape:
            print(f"  FAIL: {name}: Fortran shape {g_shape}, reference shape {w_shape}")
            ok = False
            continue

        if g_type == w_type:
            continue

        if name in _TYPE_DIFF_ALLOWED:
            noted.append(f"{name}: {g_type} against {w_type} ({_TYPE_DIFF_ALLOWED[name]})")
        else:
            print(f"  FAIL: {name}: Fortran dtype {g_type}, reference dtype {w_type}")
            ok = False

    # The values of the root scalars, not just their presence.
    #
    # Wavefront.from_genesis4 reads gridpoints and uses it to reshape, but it ignores
    # slicecount entirely -- it collects the slice groups by globbing and sorting. Genesis
    # does not ignore it: readFieldHDF5 reads slicecount and then asks for slices 1 through
    # that count, so a file understating it silently drops the tail of the pulse and a file
    # overstating it asks for a group that is not there. Nothing in the field comparison can
    # see either, so the values are checked here.

    scalars = ("gridpoints", "gridsize", "refposition", "wavelength", "slicecount",
               "slicespacing")
    with h5py.File(fortran_file, "r") as h5_f, h5py.File(reference_file, "r") as h5_r:
        for name in scalars:
            if name not in h5_f or name not in h5_r:
                continue   # Already reported as a missing object above.
            v_f = float(h5_f[name][0])
            v_r = float(h5_r[name][0])
            denom = abs(v_r) if v_r != 0 else 1.0
            if abs(v_f - v_r) / denom > 1.0e-14:
                print(f"  FAIL: {name}: Fortran wrote {v_f!r}, reference has {v_r!r}")
                ok = False

    if ok:
        n_slice = sum(1 for k in want if k.startswith("slice") and want[k][0] == "group")
        print(f"  {len(want)} objects match by name and shape, including {n_slice} slice "
              f"groups named slice000001 upward")
        print(f"  root scalars match by value: {', '.join(scalars)}")
    for note in noted:
        print(f"  note: {note}")

    return ok


def parse_check_lines(path: str) -> dict[str, float]:
    """
    Pull the 'CHECK <name> <value>' lines out of the Fortran program's stdout.
    """
    values: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) == 3 and fields[0] == "CHECK":
                values[fields[1]] = float(fields[2])
    return values


def report(name: str, a: np.ndarray, b: np.ndarray) -> float:
    """
    Compare two complex arrays and print the difference measures. Returns the peak
    normalized difference, that is, max|a-b| divided by max|b|.

    Peak normalization rather than elementwise relative error is the meaningful measure
    here. The field spans many orders of magnitude across the grid, and the round-off
    floor of an FFT based propagator is set by the largest value in the transform, not by
    the local one; an elementwise ratio in the far tail measures nothing but that floor
    divided by a number near zero.
    """
    scale = np.abs(b).max()
    if scale == 0:
        print(f"  {name}: reference is identically zero, nothing to compare")
        return np.inf

    abs_diff = np.abs(a - b)
    peak_norm = abs_diff.max() / scale

    # Elementwise relative difference, restricted to points carrying real amplitude, as a
    # cross-check that the agreement is not confined to the peak.
    mask = np.abs(b) > 1e-6 * scale
    elem_rel = (abs_diff[mask] / np.abs(b)[mask]).max() if mask.any() else np.nan

    print(f"  {name}:")
    print(f"    max |fortran - python|            = {abs_diff.max():.6e}")
    print(f"    max |python|                      = {scale:.6e}")
    print(f"    max |diff| / max |python|         = {peak_norm:.6e}   <-- largest relative difference")
    print(f"    max elementwise rel diff (>1e-6)  = {elem_rel:.6e}   ({mask.sum()} of {b.size} points)")

    return peak_norm


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    pw = sub.add_parser("write", help="Build the test wavefront and write it in Genesis4 format")
    pw.add_argument("input_file")
    pw.add_argument("--nx", type=int, default=64)
    pw.add_argument("--nz", type=int, default=5)
    pw.add_argument("--dx", type=float, default=2.0e-6)
    pw.add_argument("--dz", type=float, default=1.0e-9)
    pw.add_argument("--wavelength", type=float, default=1.0e-9)

    pc = sub.add_parser("compare", help="Drift in Python and compare against the Fortran output")
    pc.add_argument("input_file")
    pc.add_argument("fortran_file")
    pc.add_argument("--z", type=float, required=True)
    pc.add_argument("--tolerance", type=float, default=1.0e-12,
                    help="Maximum acceptable peak normalized difference")
    pc.add_argument("--fortran-log", default=None,
                    help="File holding the Fortran program's stdout. Its 'CHECK <name> <value>' "
                         "lines are recomputed here and verified.")
    pc.add_argument("--scalar-tolerance", type=float, default=1.0e-11,
                    help="Maximum acceptable relative difference for the CHECK scalars")

    args = p.parse_args()

    if args.command == "write":
        # nx == ny and dx == dy: the Genesis4 format requires a square grid.
        w = build_wavefront(nx=args.nx, ny=args.nx, nz=args.nz,
                            dx=args.dx, dy=args.dx, dz=args.dz,
                            wavelength=args.wavelength)
        w.write_genesis4(args.input_file)
        print(f"Wrote test wavefront to {args.input_file}")
        print(f"  shape {w.shape}, dx {w.dx:.6e} m, dz {w.dz:.6e} m, "
              f"wavelength {w.wavelength:.6e} m")
        print(f"  energy {w.energy:.6e} J, sigma_x {w.sigma_x:.6e} m, sigma_y {w.sigma_y:.6e} m")
        return 0

    # Compare.

    print("Python drift_wavefront reference")

    # Read the same file the Fortran program read, so that the two start from bitwise
    # identical input rather than from two separate constructions of the same field.
    w_in = Wavefront.from_genesis4(args.input_file)
    w_py = drift_wavefront(w_in, args.z)

    w_f = Wavefront.from_genesis4(args.fortran_file)

    print(f"  input:   {args.input_file}")
    print(f"  fortran: {args.fortran_file}")
    print(f"  drift:   {args.z:.6e} m")
    print(f"  shape:   python {w_py.shape}, fortran {w_f.shape}")

    if w_py.shape != w_f.shape:
        print(f"FAIL: shape mismatch, python {w_py.shape} vs fortran {w_f.shape}")
        return 1

    for name, val_py, val_f in (("dx", w_py.dx, w_f.dx),
                                ("dy", w_py.dy, w_f.dy),
                                ("dz", w_py.dz, w_f.dz),
                                ("wavelength", w_py.wavelength, w_f.wavelength)):
        if val_py != val_f:
            print(f"FAIL: {name} mismatch, python {val_py!r} vs fortran {val_f!r}")
            return 1

    print()
    print("File structure, against a file written by Wavefront.write_genesis4")
    with tempfile.TemporaryDirectory() as tmp:
        ref_file = os.path.join(tmp, "python_reference.fld.h5")
        w_py.write_genesis4(ref_file)
        structure_ok = compare_structure(args.fortran_file, ref_file)

    print()
    print("Field comparison, drifted")
    peak_norm = report("Ex", w_f.Ex, w_py.Ex)

    # An undrifted round trip, to separate an I/O error from a propagator error. If this
    # one is large but the drifted one is not, or vice versa, that says which half is wrong.
    print()
    print("Field comparison, input file read back (isolates HDF5 round trip from the drift)")
    report("Ex", Wavefront.from_genesis4(args.input_file).Ex, w_in.Ex)

    print()
    print(f"Energy: python {w_py.energy:.12e} J, fortran output file {w_f.energy:.12e} J, "
          f"relative difference {abs(w_f.energy - w_py.energy) / w_py.energy:.6e}")

    # The Fortran program's own values for the mirrored derived quantities, recomputed here
    # from the file it wrote. Without this the Fortran implementations of wavefront_energy,
    # wavefront_transverse_moments and the coordinate vectors underneath them would be
    # printed but never checked against anything.

    scalars_ok = True
    if args.fortran_log:
        print()
        print("Derived quantities computed by Fortran, recomputed here from its output file")

        # Two of these cannot reach the field comparison's round-off floor, for reasons that
        # are not defects in either implementation. Each gets a tolerance sized to its cause,
        # rather than one loose tolerance that would hide the next real problem.
        #
        # photon_energy. Bmad and scipy do not carry the same Planck constant. Bmad's
        # h_planck is 4.135667696e-15 eV s, truncated from the exact post-2019 CODATA value
        # of h/e, 4.135667696923859e-15. That is a relative difference of 2.2339e-10, which
        # is the entire discrepancy: nothing about the wavefront enters it.
        #
        # energy. A sum over every grid point, which is 20480 of them at the default size.
        # NumPy sums pairwise, with error growing as log(N) * eps; the Fortran sum intrinsic
        # sums sequentially, with error growing as N * eps, and N * eps here is 4.5e-12. The
        # two also differ in association: Fortran forms the intensity in W/m^2 and divides
        # the sum by c, where Python forms the energy density directly.
        expected = {
            "photon_energy": (w_f.photon_energy, 1.0e-9,
                              "Bmad and scipy carry different Planck constants; see the code comment"),
            "energy":        (w_f.energy, 1.0e-11,
                              "sequential against pairwise summation over all grid points"),
            "mean_x":        (w_f.mean_x, args.scalar_tolerance, ""),
            "mean_y":        (w_f.mean_y, args.scalar_tolerance, ""),
            "sigma_x":       (w_f.sigma_x, args.scalar_tolerance, ""),
            "sigma_y":       (w_f.sigma_y, args.scalar_tolerance, ""),
        }
        reported = parse_check_lines(args.fortran_log)

        missing = sorted(set(expected) - set(reported))
        if missing:
            print(f"  FAIL: the Fortran log reported no value for: {', '.join(missing)}")
            scalars_ok = False

        for name, (want, tol, note) in expected.items():
            if name not in reported:
                continue
            got = reported[name]
            denom = abs(want) if want != 0 else 1.0
            rel = abs(got - want) / denom
            status = "ok" if rel <= tol else "FAIL"
            if status == "FAIL":
                scalars_ok = False
            print(f"  {name:<14} fortran {got:+.12e}   python {want:+.12e}   "
                  f"rel {rel:.3e}  tol {tol:.1e}  {status}"
                  + (f"   [{note}]" if note else ""))

    print()
    print(f"LARGEST RELATIVE DIFFERENCE: {peak_norm:.6e}   (tolerance {args.tolerance:.1e})")

    if not np.isfinite(peak_norm) or peak_norm > args.tolerance or not scalars_ok or not structure_ok:
        print("FAIL")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
