# Steady-state FEL tracker, validated against Genesis 1.3 Version 4

Deliverable 3 of the FEL port: a single-slice steady-state FEL tracker whose physics is
transcribed from Genesis 1.3 Version 4 (GPL permits transcription), embedded in Bmad's
lattice machinery by the seam of the design brief's section 4.1, and validated against
Genesis over its `benchmark/Benchmark1-SASE` configuration from bitwise-identical starting
states.

No time dependence, no slippage, no space charge, no wakes, no harmonics, no OpenMP.
Those are later deliverables. Per-particle weights, by contrast, are carried from day one
(brief section 5): the packed arrays store one, every reduction uses it, and the
split-weight check below tests the nonuniform case that no Genesis comparison can reach.

## Files

| Path | Contents |
|---|---|
| `bsim/modules/fel_beam_mod.f90` | Packed particle slices in Bmad coordinates plus per-particle weight, Genesis `.par.h5` dump read/write (converting), copy-only `coord_struct` conversion, weighted beam diagnostics with `N_eff` |
| `bsim/modules/fel_track_mod.f90` | The transcribed FEL step: transverse push with natural focusing, RK4 ponderomotive advance, source deposition, FFT field solve; plus the transcribed Genesis interlude model |
| `bsim/fel/fel_ss_test.f90` | The tracker: walks a Bmad lattice, FEL steps in undulator segments, seam everywhere else |
| `bsim/fel/tests/run_fel_benchmark.sh` | The whole validation, one command |
| `bsim/fel/tests/compare_fel.py` | Comparison: three tiers against Genesis plus the split-weight invariance check |
| `bsim/fel/tests/Aramis-ss.in`, `Aramis.lat` | Genesis deck: Benchmark1-SASE, modified as documented in the deck header |
| `bsim/fel/tests/Aramis-1seg.in`, `Aramis-1seg.lat` | Genesis deck: one undulator segment, importing the same dumps |
| `bsim/fel/tests/aramis.bmad`, `aramis_1seg.bmad` | The Bmad lattices |

## Running

```
cd <bmad-ecosystem>
BUILD_PRODUCTION=N ./util/conda_compile                      # builds fel_ss_test
./bsim/fel/tests/run_fel_benchmark.sh [--genesis <path to genesis4>]
```

The harness runs Genesis twice (full line and single segment), the Bmad tracker four
times, and prints the largest relative difference of each check. It fails loudly if the
genesis4 binary is missing; there is no comparison without it, so there is nothing to skip
to. Genesis must be built with FFTW, since the benchmark runs with `fft_fieldsolver=true`
(the Bmad tracker transcribes the FFT solver; Genesis's default ADI solver is out of
scope).

## Architecture

Inside elements named `UND*`, `fel_track_und_step` advances the coupled system in steps of
`delz`, in Genesis's exact order: transverse half step, RK4 advance of (theta, gamma) with
the field gathered once per step, transverse half step, then source deposition and the
`exp(K2 dz)` field solve. Bmad tracking is never used inside (the brief's rule:
`symp_lie_bmad` resolves the wiggle motion the period-averaged map assumes away).

Everywhere else -- the seam -- the packed slice converts to `coord_struct`s by plain
copies, `track1_bunch` tracks them, and the field goes through `wavefront_drift`.

**Coordinates.** The packed arrays store Bmad's `(x, px/p0, y, py/p0, z, pz)` plus a
per-particle weight (macroparticle charge in Coulombs, mapping to `coord_struct%charge`).
The ponderomotive phase is derived, not stored: Genesis's per-particle theta splits into
a common reference advance -- one scalar per beam, `phi0`, advanced once per step -- and
the particle-specific lag carried by Bmad's z:

```
theta_j = phi0 - ks*tau_j,    tau_j = -z_j/beta_j = c*(t_j - t_ref)
```

This is the reference-offset formulation the design brief's section 8 identifies as the
FP32-safe one, and it removes the brief's 6.4 hazard outright: z does not wrap, so slice
migration has no theta-wrap-plus-index update to get wrong. The longitudinal RK4 still
runs in Genesis's (theta, gamma) chart as per-step working variables -- theta <-> z is an
affine map, under which RK4 is exactly invariant, so the verbatim transcription survives;
the per-step gamma <-> pz round trip costs ~1 ulp of gamma per step, which is what moved
tier1 from 2.5e-12 to 2.8e-11 when this representation replaced the Genesis-coordinate
one. The particles live in the packed arrays at all times except inside `track1_bunch`;
`coord_struct` never enters the FEL step loop.

**Weights.** Every reduction is weighted: the source deposition scales per particle as
`c*w_j/slicelength` (Genesis's `current/N` for uniform weights), bunching is
`|sum w e^(i theta)|/sum w`, and `N_eff = (sum w)^2/sum w^2` is a per-slice diagnostic.
The Genesis dump format carries no weights, so imports are uniform and a Genesis
comparison can only exercise the uniform case; the split-weight check below covers the
rest.

The field inside the tracker is kept in Genesis's internal units, converting once at each
program boundary, so every coefficient of the dynamics is computed from the same numbers
Genesis computes it from. Genesis's own constants (`vacimp = 376.73`, note the truncation,
and `eev = 510998.95069`) are transcribed alongside the formulas and deliberately not
"corrected" to Bmad's values.

Undulator segments are marked by name (`UND*`) and their FEL parameters come from the
namelist, not from the lattice file. A real FEL element type with its own tracking method
is a later deliverable; this is the smallest scheme that exercises the seam.

## Validation: four checks, from one command

Both codes start from the same Genesis `&write` dumps, so the initial state is bitwise
identical and no loader is reproduced. Genesis records diagnostics once at the start and
once per integration step, with each interlude element being a single step; the tracker
records at the same z positions. Measured, on the numbers this tree was developed against:

| Tier | What runs | Largest relative difference |
|---|---|---|
| `tier1` | One undulator segment: the FEL core alone | **2.8e-11** (power curve 1.9e-12; per-particle final gamma 3.2e-15) |
| `tier2_genesis` | Full 6-FODO line, interludes via the transcribed Genesis model | **5.9e-8** (power curve 7.3e-9) |
| `tier2_bmad` | Full line, interludes via the Bmad seam | **5.0e-2** (power curve 1.3e-2) -- a measured model difference, see below |
| `weight_split` | tier1 rerun with every particle split into coincident w/3 + 2w/3 copies, against the unsplit run | **6.2e-13** |

Particle ordering is preserved by both codes in steady state, so the final dumps compare
particle by particle, not just statistically.

The `weight_split` check is Fortran against itself and exists because Genesis cannot test
the weighted paths: its dump format has no weights, so every cross-code comparison sees
the uniform case, where a bug like using one particle's weight for all is invisible.
Collective observables are linear in the weights, so the split run must reproduce the
unsplit one to round-off; the weight-dropping mutation fails this check at 2.4e-1.

### The tier2_bmad difference is a transport model difference, located and priced

The divergence localizes to the quadrupoles: after the first quad the bunching phase steps
by 2.1e-3 rad while the transverse beam sizes still agree to 2.5e-11. The mechanism:
Genesis advances theta through an interlude element as a single step whose path-length
term `(px^2+py^2)/(2 gamma^2)` is sampled once, at mid element (transverse half step, then
the theta step, then the second half; `TrackBeam::track` + the `BeamSolver` drift case).
Bmad's z advance integrates the same term exactly through the quad map. In a quad px
changes substantially across the element, so midpoint sampling differs from the integral
at the ~10 percent level of the ~0.01 rad px^2 contribution -- the observed 2e-3 rad per
quad. Through exponential gain and twelve quads this compounds to the 1.3e-2 power
difference at saturation.

The proof is `tier2_genesis`: the identical build with only the interlude model swapped to
the transcribed Genesis step collapses the difference by six orders of magnitude. Nothing
else changes between those runs, so nothing else contributes at that level. Where the two
transport models differ, Bmad's is the better one (the exact integral); the benchmark's
job is to price the difference against Genesis, not to prefer Genesis's answer.

Residual budget for `tier2_genesis`'s own 1e-8: rounding differences in the interlude
theta advance. Genesis evaluates `ks*(1 - 1/beta_z)` per particle -- a ~4e-9 cancellation
carrying ~1e-16 absolute rounding -- inside its RK4 bookkeeping; the transcription
evaluates the same expression but sums the step differently (the RK4 collapses exactly
when the slope is theta independent, and is collapsed). Multiplied by `ks*L` this is
microradians of per-particle phase noise per interlude, amplified through gain. The
per-particle theta medians tell the same story: 6.4e-14 rad for the typical particle, with
a chaotic separatrix tail (see below).

### theta is reported, not gated

The final per-particle theta difference is printed with max, rms and median but does not
gate the comparison. theta is a bucket phase: at saturation neighboring trajectories
separate exponentially, so the worst-particle difference measures Lyapunov amplification,
not implementation quality (the brief's 9.1 warning about chaotic growth, met in
practice: `tier2_genesis` has median 6.4e-14 rad against max 2.0e-5). Its collective
effect is gated, through the bunching curve and the final field.

### The harness bites

Verified by mutation, the FINDINGS.md 4.1 discipline: dropping the factor 2 on the source
term fails tier1 at 3.5e-1; dropping the conjugation in the field gather fails at 2.8e-2;
replacing `sqrt(faw2)` by `faw` in the deposition -- nearly degenerate for this beam, the
two differ at second order in the transverse offsets -- still fails, at 1.5e-10 against
the 1e-10 gate; and using one particle's weight for every deposition, invisible to every
Genesis-based tier, fails the split-weight check at 2.4e-1.

## Facts about Genesis this work pinned down

- Outside undulators Genesis does not subdivide into `delz` steps: each interlude element
  is one integration step of the element's full length (`Lattice::unrollLattice` pushes
  one entry per non-undulator layout segment). The step count over the benchmark is
  12*89 undulator steps + 36 interlude steps = 1104.
- `slippage` and `phaseshift` are no-ops for this benchmark: `Control::applySlippage`
  returns immediately when not time dependent, and the phaseshift array is all zero
  without phase shifter elements.
- A helical undulator defaults to `kx = ky = 0.5` (LatticeParser.cpp:329), scaled by
  `ku^2` in the unroll.
- The `&importbeam` / `&importfield` namelists take full filenames and make the
  shared-start methodology possible; `&write` before `&track` produces the dumps.
- Genesis's internal field unit: `dfl [sqrt(W)] = u * dgrid * eev / (ks * sqrt(vacimp))`
  (writeFieldHDF5.cpp:70), with `vacimp = 376.73` truncated and `eev = 510998.95069`.
- The quad transport is chromatic through per-particle `gammaz` with `foc^2 =
  k1*gamma0/gammaz`; Bmad's equivalent scaling is `k1*p0/p`. The two differ by ~4e-9
  relative at this energy (1 - beta0), far below the path-length-term difference.
