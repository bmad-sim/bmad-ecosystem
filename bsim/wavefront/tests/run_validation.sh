#!/bin/bash
#
# Validate the Fortran free-space wavefront propagator.
#
# One command runs everything and prints the largest relative difference found:
#
#   ./run_validation.sh
#
# Two independent comparisons are made.
#
# 1. Against openPMD-beamphysics. A test wavefront is built in Python and written in
#    Genesis4 HDF5 format; the Fortran program reads it, drifts it and writes it back; the
#    same input is drifted in Python with drift_wavefront and the two results compared. This
#    exercises the whole path -- HDF5 read, propagation, HDF5 write -- against an
#    independent implementation, and it is the comparison the deliverable is defined by.
#    Several grid sizes and drift distances are run, including an odd number of transverse
#    points and a negative drift.
#
# 2. Against an in-Fortran reference propagator, wavefront_drift_reference, which builds its
#    own kernel and performs a direct DFT instead of using FFTW. This exists because the
#    Genesis4 format requires nx = ny and dx = dy, and on a square grid with equal spacings
#    the propagation kernel is symmetric under interchanging the two transverse axes, so
#    comparison 1 structurally cannot see a transposed transform or an interchanged dx and
#    dy. Both were confirmed by mutation to pass comparison 1 at round-off. This one runs on
#    rectangular grids with unequal spacings and catches them.
#
# Prerequisites:
#
#   1. Bmad built, so that debug/bin/wavefront_drift_test (or production/bin) exists:
#        BUILD_PRODUCTION=N ./util/conda_compile       # from the bmad-ecosystem root
#
#   2. The Python environment defined alongside this script:
#        conda env create -f environment.yml
#
#   3. An openPMD-beamphysics checkout. Found automatically if it sits beside the
#      bmad-ecosystem checkout; otherwise pass --beamphysics <path>.
#
# Options:
#   --beamphysics <path>  openPMD-beamphysics checkout. Default: sibling of bmad-ecosystem.
#   --python <path>       Python interpreter. Default: the bmad-fel-validate conda env.
#   --exe <path>          wavefront_drift_test binary. Default: debug then production.
#   --work-dir <path>     Where to put the HDF5 files. Default: a temporary directory.
#   --tolerance <x>       Maximum acceptable relative difference. Default: 1e-12
#   --quick               Run only the first Python comparison case.
#
# Exit status is zero only if every comparison passes.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

BEAMPHYSICS=""
PYTHON=""
EXE=""
WORK_DIR=""
TOLERANCE="1e-12"
QUICK=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --beamphysics) BEAMPHYSICS="$2"; shift 2 ;;
    --python)      PYTHON="$2";      shift 2 ;;
    --exe)         EXE="$2";         shift 2 ;;
    --work-dir)    WORK_DIR="$2";    shift 2 ;;
    --tolerance)   TOLERANCE="$2";   shift 2 ;;
    --quick)       QUICK=1;          shift 1 ;;
    -h|--help)     sed -n '2,50p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

# Locate the Fortran binary. Debug first, since that is the build to develop against.

if [[ -z "$EXE" ]]; then
  for candidate in "$BMAD_ROOT/debug/bin/wavefront_drift_test" \
                   "$BMAD_ROOT/production/bin/wavefront_drift_test"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi

if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: wavefront_drift_test not found." >&2
  echo "Build it with:  cd $BMAD_ROOT && BUILD_PRODUCTION=N ./util/conda_compile" >&2
  echo "Or point at it with --exe <path>." >&2
  exit 1
fi

# Locate openPMD-beamphysics.

if [[ -z "$BEAMPHYSICS" ]]; then
  for candidate in "$BMAD_ROOT/../openPMD-beamphysics" "$HOME/Code/GitHub/openPMD-beamphysics"; do
    if [[ -d "$candidate/beamphysics/wavefront" ]]; then
      BEAMPHYSICS="$(cd "$candidate" && pwd)"
      break
    fi
  done
fi

if [[ -z "$BEAMPHYSICS" || ! -d "$BEAMPHYSICS/beamphysics/wavefront" ]]; then
  echo "Error: openPMD-beamphysics checkout not found (looked for beamphysics/wavefront)." >&2
  echo "Pass it with --beamphysics <path>." >&2
  exit 1
fi

# Locate Python.

if [[ -z "$PYTHON" ]]; then
  for candidate in "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$HOME/Code/miniforge3/envs/bmad-fel-validate/bin/python3"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi

if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "Error: Python interpreter not found." >&2
  echo "Create the environment with:  conda env create -f $SCRIPT_DIR/environment.yml" >&2
  echo "Or point at an interpreter with --python <path>." >&2
  exit 1
fi

# Work directory.

CLEAN_WORK_DIR=0
if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/wavefront_validation.XXXXXX")"
  CLEAN_WORK_DIR=1
fi
mkdir -p "$WORK_DIR"

cleanup() {
  if [[ "$CLEAN_WORK_DIR" == "1" ]]; then rm -rf "$WORK_DIR"; fi
}
trap cleanup EXIT

export PYTHONPATH="$BEAMPHYSICS${PYTHONPATH:+:$PYTHONPATH}"

echo "=============================================================================="
echo "Wavefront free-space drift validation"
echo "=============================================================================="
echo "  Fortran binary:      $EXE"
echo "  Python:              $PYTHON"
echo "  openPMD-beamphysics: $BEAMPHYSICS"
echo "  Work directory:      $WORK_DIR"
echo "  Tolerance:           $TOLERANCE"
echo

# ---------------------------------------------------------------------------------
# Comparison 1: against openPMD-beamphysics, through the Genesis4 file format.
#
# Case columns: label, nx (= ny), nz, drift [m].
#
#   even_grid       power of two grid, the usual case
#   odd_grid        odd nx, which indexes the negative wavenumber half of the FFT
#                     ordering differently from even nx
#   non_power_of_2  a size FFTW must handle with mixed radix rather than radix 2
#   negative_drift  backward propagation, which must be the inverse of forward
#   long_drift      far enough that the field diffracts substantially over the grid
#   zero_drift      a no-op, which must return the input unchanged
# ---------------------------------------------------------------------------------

CASES=(
  "even_grid       64  5   0.05"
  "odd_grid        63  4  -0.03"
  "non_power_of_2  48  3   0.02"
  "long_drift      64  2   0.40"
  "zero_drift      32  3   0.0"
)

if [[ "$QUICK" == "1" ]]; then
  CASES=("${CASES[0]}")
fi

OVERALL_STATUS=0
WORST="0"
SUMMARY=()

for case_spec in "${CASES[@]}"; do
  read -r LABEL NX NZ Z <<<"$case_spec"

  IN_FILE="$WORK_DIR/${LABEL}_in.fld.h5"
  OUT_FILE="$WORK_DIR/${LABEL}_fortran.fld.h5"
  LOG_FILE="$WORK_DIR/${LABEL}_fortran.log"

  echo "=============================================================================="
  echo "Case $LABEL:  nx = ny = $NX,  nz = $NZ,  drift = $Z m"
  echo "=============================================================================="

  echo "--- build the test wavefront in Python --------------------------------------"
  if ! "$PYTHON" "$SCRIPT_DIR/validate_drift.py" write "$IN_FILE" --nx "$NX" --nz "$NZ"; then
    echo "Case $LABEL FAILED: could not write the input file"
    SUMMARY+=("$LABEL FAIL (input)")
    OVERALL_STATUS=1
    continue
  fi
  echo

  echo "--- read, drift and write in Fortran ----------------------------------------"
  if ! "$EXE" "$IN_FILE" "$OUT_FILE" "$Z" 2>&1 | tee "$LOG_FILE"; then
    echo "Case $LABEL FAILED: the Fortran program returned nonzero"
    SUMMARY+=("$LABEL FAIL (fortran)")
    OVERALL_STATUS=1
    continue
  fi
  echo

  echo "--- drift the same input in Python and compare ------------------------------"
  CMP_OUT="$("$PYTHON" "$SCRIPT_DIR/validate_drift.py" compare "$IN_FILE" "$OUT_FILE" \
                --z "$Z" --tolerance "$TOLERANCE" --fortran-log "$LOG_FILE" 2>&1)"
  CMP_STATUS=$?
  echo "$CMP_OUT"
  echo

  CASE_DIFF="$(echo "$CMP_OUT" | sed -n 's/^LARGEST RELATIVE DIFFERENCE: *\([^ ]*\).*/\1/p' | head -1)"
  [[ -z "$CASE_DIFF" ]] && CASE_DIFF="nan"

  if [[ $CMP_STATUS -eq 0 ]]; then
    SUMMARY+=("$LABEL pass  $CASE_DIFF")
  else
    SUMMARY+=("$LABEL FAIL  $CASE_DIFF")
    OVERALL_STATUS=1
  fi

  WORST="$(WORST_A="$WORST" WORST_B="$CASE_DIFF" "$PYTHON" -c '
import math, os
a = float(os.environ["WORST_A"])
b = float(os.environ["WORST_B"])
b = b if math.isfinite(b) else float("inf")
print(repr(max(a, b)))')"
done

# ---------------------------------------------------------------------------------
# Comparison 2: against the in-Fortran direct DFT reference, on rectangular grids.
# ---------------------------------------------------------------------------------

echo "=============================================================================="
echo "Rectangular grid check: FFT propagator against the direct DFT reference"
echo "=============================================================================="
SELF_OUT="$("$EXE" -self_check 2>&1)"
SELF_STATUS=$?
echo "$SELF_OUT"
echo

SELF_DIFF="$(echo "$SELF_OUT" | sed -n 's/.*LARGEST RELATIVE DIFFERENCE.*: *\([^ ]*\).*/\1/p' | head -1)"
[[ -z "$SELF_DIFF" ]] && SELF_DIFF="nan"

if [[ $SELF_STATUS -eq 0 ]]; then
  SUMMARY+=("rectangular_grid pass  $SELF_DIFF")
else
  SUMMARY+=("rectangular_grid FAIL  $SELF_DIFF")
  OVERALL_STATUS=1
fi

WORST="$(WORST_A="$WORST" WORST_B="$SELF_DIFF" "$PYTHON" -c '
import math, os
a = float(os.environ["WORST_A"])
b = float(os.environ["WORST_B"])
b = b if math.isfinite(b) else float("inf")
print(repr(max(a, b)))')"

# ---------------------------------------------------------------------------------

echo "=============================================================================="
echo "Summary"
echo "=============================================================================="
for line in "${SUMMARY[@]}"; do
  printf '  %s\n' "$line"
done
echo
printf '  LARGEST RELATIVE DIFFERENCE OVER ALL CASES: %s   (tolerance %s)\n' "$WORST" "$TOLERANCE"
echo

if [[ $OVERALL_STATUS -eq 0 ]]; then
  echo "PASS"
else
  echo "FAIL"
fi
echo "=============================================================================="

exit $OVERALL_STATUS
