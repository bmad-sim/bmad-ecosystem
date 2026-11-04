#!/bin/bash



if ! command -v conda &>/dev/null; then
  echo "Error: conda binary not found"
  exit 1
fi

CONDA_ROOT=$(conda info --base)
source "$CONDA_ROOT/etc/profile.d/conda.sh"

conda activate bmad-build

# Activating from another env swaps its bin dir in place, which can leave it behind
# /opt/local/bin (MacPorts) or Homebrew, so force the env's tools first.
export PATH="$CONDA_PREFIX/bin:$PATH"

# conda-forge's ld64/tapi can't parse the macOS 27 SDK's .tbd stubs (new arm64e.x1 arch),
# so build against the 26.5 SDK until conda-forge catches up.
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk

# compilers>=2 no longer exports CC/CXX on activation, so CMake would fall back to Apple's
# /usr/bin/cc (no -fopenmp). The prefixed clang reads its .cfg for conda include/lib/rpath.
export CC=arm64-apple-darwin20.0.0-clang
export CXX=arm64-apple-darwin20.0.0-clang++
export FC=gfortran




export LDFLAGS="$LDFLAGS -Wl,-no_compact_unwind"

echo "** CONDA PREFIX IS: $CONDA_PREFIX"

BUILD_PRODUCTION=${BUILD_PRODUCTION=Y}
ARTIFACT_FOLDER="production"

git checkout origin/main -- util/dist_prefs

cat <<EOF >>util/dist_prefs
export DIST_F90_REQUEST="gfortran"
export ACC_PLOT_PACKAGE="pgplot"
export ACC_PLOT_DISPLAY_TYPE="X"
export ACC_ENABLE_OPENMP="Y"
export ACC_ENABLE_MPI="N"
export ACC_FORCE_BUILTIN_MPI="N"
export ACC_ENABLE_GFORTRAN_OPTIMIZATION="$BUILD_PRODUCTION"
export ACC_ENABLE_SHARED="Y"
export ACC_ENABLE_SHARED_ONLY="N"
export ACC_ENABLE_FPIC="Y"
export ACC_ENABLE_PROFILING="N"
export ACC_SET_GMAKE_JOBS="12"
export ACC_CONDA_BUILD="Y"
export ACC_CONDA_BUILD_TESTS="Y"
export ACC_CONDA_PATH="$CONDA_PREFIX"
export ACC_USE_MACPORTS="N"
EOF

echo "**** Invoking dist_source_me"
source util/dist_source_me

cp $CONDA_PREFIX/share/gnuconfig/config.* .

# build production if BUILD_PRODUCTION is set to Y
if [[ "$BUILD_PRODUCTION" == "Y" ]]; then
  echo "**** Invoking dist_build_production"
  util/dist_build_production
else
  echo "**** Invoking dist_build_debug"
  util/dist_build_debug
fi

# --- Optional: build code_examples (not included in main build) ---
# export LIBRARY_PATH=$CONDA_PREFIX/lib
# cd code_examples/production
# cmake ..
# make plot_continuous_color-exe
