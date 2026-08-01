#!/bin/bash
#
# Validate the Bmad steady-state FEL tracker against Genesis 1.3 Version 4.
#
# One command runs everything: Genesis over the Benchmark1-SASE lattice (writing the
# initial dumps both codes start from), Genesis over a single undulator segment importing
# the same dumps, the Bmad tracker over both lattices (the full line twice, once with the
# Bmad seam and once with the transcribed Genesis interlude model), and a three-tier
# comparison printing the largest relative difference of each tier.
#
#   ./run_fel_benchmark.sh
#
# Prerequisites:
#
#   1. Bmad built, so debug/bin/fel_ss_test (or production/bin) exists:
#        BUILD_PRODUCTION=N ./util/conda_compile        # from the bmad-ecosystem root
#
#   2. A genesis4 binary (CPU is fine). Default location below; override with --genesis.
#      Built per the Genesis repository's instructions; the benchmark needs FFTW support
#      compiled in, since it runs with fft_fieldsolver = true.
#
#   3. Python with numpy and h5py: the bmad-fel-validate environment
#      (conda env create -f ../../wavefront/tests/environment.yml).
#
# Options:
#   --genesis <path>    genesis4 binary. Default: ~/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4
#   --exe <path>        fel_ss_test binary. Default: debug then production.
#   --python <path>     Python interpreter. Default: the bmad-fel-validate conda env.
#   --work-dir <path>   Where to run. Default: a temporary directory (kept on failure).
#
# Exit status is zero only if every tier passes its tolerance.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

GENESIS="$HOME/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4"
EXE=""
PYTHON=""
WORK_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)  GENESIS="$2";  shift 2 ;;
    --exe)      EXE="$2";      shift 2 ;;
    --python)   PYTHON="$2";   shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    -h|--help)  sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

# The comparison is against Genesis; without the binary there is nothing to compare and
# the only correct behavior is to fail loudly, not to skip.

if [[ ! -x "$GENESIS" ]]; then
  echo "Error: genesis4 binary not found at: $GENESIS" >&2
  echo "Point at one with --genesis <path>." >&2
  exit 1
fi

if [[ -z "$EXE" ]]; then
  for candidate in "$BMAD_ROOT/debug/bin/fel_ss_test" "$BMAD_ROOT/production/bin/fel_ss_test"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: fel_ss_test not found. Build with: cd $BMAD_ROOT && BUILD_PRODUCTION=N ./util/conda_compile" >&2
  exit 1
fi

if [[ -z "$PYTHON" ]]; then
  for candidate in "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$HOME/Code/miniforge3/envs/bmad-fel-validate/bin/python3"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "Error: Python not found. Create the environment:" >&2
  echo "  conda env create -f $BMAD_ROOT/bsim/wavefront/tests/environment.yml" >&2
  exit 1
fi

KEEP_WORK_DIR=1
if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/fel_benchmark.XXXXXX")"
  KEEP_WORK_DIR=0
fi
mkdir -p "$WORK_DIR"

echo "=============================================================================="
echo "FEL steady-state benchmark"
echo "=============================================================================="
echo "  fel_ss_test: $EXE"
echo "  genesis4:    $GENESIS"
echo "  python:      $PYTHON"
echo "  workdir:     $WORK_DIR"
echo

cp "$SCRIPT_DIR/Aramis-ss.in" "$SCRIPT_DIR/Aramis.lat" \
   "$SCRIPT_DIR/Aramis-1seg.in" "$SCRIPT_DIR/Aramis-1seg.lat" \
   "$SCRIPT_DIR/aramis.bmad" "$SCRIPT_DIR/aramis_1seg.bmad" "$WORK_DIR/"

cd "$WORK_DIR" || exit 1

# Without FI_PROVIDER=tcp the MPI runtime's provider search adds tens of seconds
# (FINDINGS.md / design brief 11.1).
export FI_PROVIDER=tcp

echo "--- Genesis: full line (writes the shared initial dumps) ---------------------"
if ! "$GENESIS" Aramis-ss.in > genesis-full.log 2>&1; then
  echo "Genesis full-line run FAILED; log tail:" >&2
  tail -20 genesis-full.log >&2
  exit 1
fi
tail -3 genesis-full.log
echo

echo "--- Genesis: single segment (imports the same dumps) -------------------------"
if ! "$GENESIS" Aramis-1seg.in > genesis-1seg.log 2>&1; then
  echo "Genesis single-segment run FAILED; log tail:" >&2
  tail -20 genesis-1seg.log >&2
  exit 1
fi
tail -3 genesis-1seg.log
echo

make_nml () {
  cat > "$1" <<NML
&fel_ss_params
  lat_file = "$2"
  beam_file = "Aramis-initial.par.h5"
  field_file = "Aramis-initial.fld.h5"
  out_root = "$3"
  gamma0 = 11357.82
  delz = 0.045
  und_aw = 0.84853
  und_lambdau = 0.015
  und_kx = 0.5, und_ky = 0.5
  und_helical = T
  interlude_model = "$4"
&end
NML
}

make_nml tier1.nml  aramis_1seg.bmad tier1  bmad
make_nml tier2.nml  aramis.bmad      tier2  bmad
make_nml tier2g.nml aramis.bmad      tier2g genesis

for tier in tier1 tier2 tier2g; do
  echo "--- fel_ss_test: $tier -------------------------------------------------------"
  if ! "$EXE" $tier.nml > fel-$tier.log 2>&1; then
    echo "fel_ss_test $tier FAILED; log tail:" >&2
    tail -20 fel-$tier.log >&2
    exit 1
  fi
  tail -4 fel-$tier.log
  echo
done

"$PYTHON" "$SCRIPT_DIR/compare_fel.py" "$WORK_DIR"
STATUS=$?

if [[ $STATUS -eq 0 && $KEEP_WORK_DIR -eq 0 ]]; then
  rm -rf "$WORK_DIR"
elif [[ $STATUS -ne 0 ]]; then
  echo "Outputs kept in: $WORK_DIR"
fi

exit $STATUS
