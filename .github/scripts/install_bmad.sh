#!/bin/bash
set -e

# echo "******* DEBUG ********"
# dpkg -l pkg-config
# ls -l /usr/share/aclocal/pkg.m4
# aclocal --print-ac-dir
# export ACLOCAL_PATH=/usr/share/aclocal
# echo "******* DEBUG ********"

echo "**** Setup Preferences"

cat <<EOF >> ./util/dist_prefs
export DIST_F90_REQUEST="gfortran"
export ACC_PLOT_PACKAGE="pgplot"
export ACC_PLOT_DISPLAY_TYPE="X"
export ACC_ENABLE_OPENMP="$USE_MPI"
export ACC_ENABLE_MPI="$USE_MPI"
export ACC_FORCE_BUILTIN_MPI="N"
export ACC_ENABLE_GFORTRAN_OPTIMIZATION="Y"
export ACC_ENABLE_SHARED="$SHARED"
export ACC_ENABLE_SHARED_ONLY="$SHARED"
export ACC_ENABLE_FPIC="Y"
export ACC_ENABLE_PROFILING="N"
export ACC_SET_GMAKE_JOBS="2"
EOF

echo "**** Invoking dist_source_me"
source ./util/dist_source_me

echo "**** Invoking dist_build_production"
./util/dist_build_production
