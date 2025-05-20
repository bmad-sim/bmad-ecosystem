#!/bin/bash

# echo "******* DEBUG ********"
# dpkg -l pkg-config
# ls -l /usr/share/aclocal/pkg.m4
# aclocal --print-ac-dir
# export ACLOCAL_PATH=/usr/share/aclocal
# echo "******* DEBUG ********"

echo "**** Setup Preferences"

echo "Number of processors: $(nproc)"

cat <<EOF >>./util/dist_prefs
export DIST_F90_REQUEST="gfortran"
export ACC_PLOT_PACKAGE="pgplot"
export ACC_PLOT_DISPLAY_TYPE="X"
export ACC_ENABLE_OPENMP="$USE_MPI"
export ACC_ENABLE_MPI="$USE_MPI"
export ACC_FORCE_BUILTIN_MPI="Y"
export ACC_ENABLE_GFORTRAN_OPTIMIZATION="Y"
export ACC_ENABLE_SHARED="N"
export ACC_ENABLE_SHARED_ONLY="$SHARED"
export ACC_ENABLE_FPIC="Y"
export ACC_ENABLE_PROFILING="N"
export ACC_SET_GMAKE_JOBS="$(nproc)"
EOF

if [ "$USE_CONDA" == "1" ]; then
  [ -z "$CONDA_PREFIX" ] && {
    echo "CONDA_PREFIX unset?"
    exit 1
  }

  echo "* Using conda to build: $CONDA_PREFIX"

  cat <<EOF >>./util/dist_prefs
export ACC_CONDA_BUILD="Y"
export ACC_CONDA_BUILD_TESTS="Y"
export ACC_CONDA_PATH="$CONDA_PREFIX"
export ACC_USE_MACPORTS="N"
EOF

fi

echo "**** Invoking dist_source_me"
source ./util/dist_source_me

echo "**** Invoking dist_build_production"
./util/dist_build_production
