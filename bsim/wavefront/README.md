# Fortran wavefront: representation, free-space propagation, Genesis4 HDF5

A `wavefront_struct` holding the complex transverse electric field of a radiation pulse, a
free-space propagator built on FFTW, and Genesis 1.3 Version 4 HDF5 input and output. It
mirrors the `Wavefront` class of
[openPMD-beamphysics](https://github.com/ChristopherMayes/openPMD-beamphysics)
(`beamphysics/wavefront/wavefront.py`) closely enough that the two can be compared field by
field, which is how it is validated.

There is no FEL physics here: no particles, no deposition, no undulator. This is the object
that a later FEL tracker will carry, and the transport step it will use outside undulators.

## Files

| Path | Contents |
|---|---|
| `bsim/modules/wavefront_mod.f90` | `wavefront_struct`, derived quantities, `wavefront_drift`, the cached FFTW plan, and `wavefront_drift_reference` |
| `bsim/modules/wavefront_hdf5_mod.f90` | `wavefront_read_genesis4`, `wavefront_write_genesis4` |
| `bsim/wavefront/wavefront_drift_test.f90` | Test driver, both modes |
| `bsim/wavefront/tests/run_validation.sh` | The whole validation, one command |
| `bsim/wavefront/tests/validate_drift.py` | Python side: builds the test input, drifts it, compares |
| `bsim/wavefront/tests/environment.yml` | Conda environment for the Python side |

## Building and running

```
cd <bmad-ecosystem>
conda env create -n bmad-build -f .github/bmad-build-env.yaml   # once
BUILD_PRODUCTION=N ./util/conda_compile

conda env create -f bsim/wavefront/tests/environment.yml        # once
./bsim/wavefront/tests/run_validation.sh
```

The script finds the binary, the Python interpreter and the openPMD-beamphysics checkout on
its own if they are in the usual places, and takes `--exe`, `--python` and `--beamphysics`
if they are not. It exits nonzero unless every comparison passes.

## What the validation checks

Two independent comparisons, because neither alone is sufficient.

**Against openPMD-beamphysics.** A test wavefront is built in Python and written in Genesis4
format; Fortran reads it, drifts it and writes it back; the same input is drifted in Python
with `drift_wavefront`; the two results are compared. This exercises the whole path against
an independent implementation. Five cases run: an even grid, an odd grid, a non-power-of-two
grid, a long drift and a zero drift. The structure and the root scalar values of the Fortran
output file are also compared against a file written by `Wavefront.write_genesis4`, because
`from_genesis4` globs and sorts the slice groups and so cannot see an error in their names or
in `slicecount`, while Genesis, which addresses them by index, would be broken by one.

**Against an in-Fortran reference propagator.** `wavefront_drift_reference` builds its own
kernel and performs a direct DFT rather than using FFTW. It exists because the Genesis4
format requires `nx = ny` and `dx = dy`, and on a square grid with equal spacings the
propagation kernel is symmetric under interchanging the two transverse axes. A transposed
transform, or an interchanged `dx` and `dy`, is therefore invisible to the first comparison
no matter how asymmetric the test field is — both were confirmed by mutation to pass it at
round-off. The reference runs on rectangular grids with unequal spacings and catches them.

Current result: agreement at **8.7e-16**, about four ulp of double precision.

The test field is deliberately asymmetric — off axis, elliptical, tilted, with a steering
phase and opposite curvature in the two planes, varying slice to slice. A centered symmetric
Gaussian is invariant under most of the mistakes worth catching.

### Mutation results

Every mutation below was built and run against the harness. All are caught.

| Mutation | Caught by |
|---|---|
| kernel sign flipped | both |
| kernel `4*pi` written as `2*pi` | both |
| `dx` and `dy` interchanged in the kernel | reference |
| negative half of the FFT wavenumber ordering dropped | both |
| two arguments of `fftw_plan_dft_2d` interchanged | reference only |
| `1/(nx*ny)` normalization dropped | both |
| slice transposed on HDF5 write, or on read | Python |
| `dx` factor dropped from the dfl scaling, either direction | Python |
| slice group names off by one, or zero based | file structure |
| `slicecount` off by one | file structure |
| imaginary part negated on write | Python |

One mutation is not caught and cannot be: interchanging `wf_fft_forward$` and
`wf_fft_backward$` leaves the result bit-identical, because the kernel is even in `k`. See
the comment at their declaration.

Separately, the reader checks each slice dataset's extent against `gridpoints` before reading
it, because `hdf5_read_dataset_real` goes through `H5LTread_dataset`, which takes no buffer
size and writes however much the file claims. A file whose `gridpoints` disagrees with its
slice datasets is refused by name instead of overrunning the buffer, which on a forgiving
allocator would read as working.

### Two differences from Python that are not defects

Both are reported by the harness with their cause, and each has its own tolerance rather than
being absorbed into one loose one.

`photon_energy` differs by 2.234e-10. Bmad's `h_planck` is `4.135667696e-15` eV·s, truncated
from the exact post-2019 CODATA value of `h/e`, `4.135667696923859e-15`. That accounts for
the whole difference.

`energy` differs by about 1e-12. It is a sum over every grid point, 20480 of them at the
default size. NumPy sums pairwise, with error growing as `log(N)·eps`; the Fortran `sum`
intrinsic sums sequentially, with error growing as `N·eps`, which is 4.5e-12 here.

### A finding about the Python writer

`Wavefront.write_genesis4` writes `refposition` as **int64** whenever the default argument is
used, because the default is the Python int `0` and `np.asarray([0])` infers an integer type.
Pass a float and it writes float64. Genesis writes this field as a double with a unit
attribute of `"m"`, and reads it with `H5Dread` into `H5T_NATIVE_DOUBLE`, so the integer file
is read correctly and nothing breaks — but a position in meters is being stored as an
integer. This Fortran writer writes float64 deliberately. The harness reports the difference
rather than failing on it. Worth fixing on the Python side.

## Design notes

**Field precision is one parameter.** `wf_rp` in `wavefront_mod` is the kind of every field
array and of every argument that carries field values. It is `rp` (double). Changing it to
single will not silently work: the plan cache's work buffer is declared `complex(wf_rp)` and
the double precision FFTW entry points are bound, so a single precision field fails to
compile at `fftw_execute_dft` with a type mismatch. That is the intended behaviour — the
`fftwf_*` entry points have to be bound deliberately, not fallen into.

**Nonzero curvature is refused, not ignored.** `wavefront_drift` accepts a `curvature`
argument and errors on a nonzero value. Curvature-corrected propagation with grid rescaling,
as in `drift_wavefront_advanced`, is a different propagator.

**The transverse boundary is periodic and unapodised**, because that is what a discrete
Fourier transform imposes. Field reaching the edge of the grid wraps around and re-enters.
Over a long drift the grid has to be large enough to contain the field. This is the same
limitation Genesis has.

**The transform is single threaded on purpose.** `fftw_plan_with_nthreads` is not called: the
parallelisation axis for a wavefront is the slice index, and threading a small 2D transform
underneath a loop over slices would only oversubscribe. `wavefront_fft2` is not thread safe,
since its plan cache and work buffer are module state; parallelising over slices wants one
cache per thread, which is a change confined to six variables and one routine.

**Slices are contiguous.** `Ex(:,:,iz)` is one slice and is contiguous in Fortran's column
major order, which is what both a per-slice 2D FFTW plan and a parallel loop over slices
want. The Python class uses the same `(nx, ny, nz)` shape, but in C order, so there its
slices are strided. The same property is why the Genesis4 per-slice datasets, which store x
fastest, need only a `reshape` here where the Python implementation needs a transpose.

## Not included

`pad`, `crop`, `auto_crop`, the k-space representation, `estimate_curvature`, the
openPMD-standard file format, and curvature-corrected propagation all exist on the Python
side and are not mirrored here.
