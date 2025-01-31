#!/bin/bash

# echo "******* DEBUG ********"
# dpkg -l pkg-config
# ls -l /usr/share/aclocal/pkg.m4
# aclocal --print-ac-dir
# export ACLOCAL_PATH=/usr/share/aclocal
# echo "******* DEBUG ********"

echo "# Installing bmad"

echo "## Setup Preferences"

echo "Number of processors: $(nproc)"

cat <<EOF >>./util/dist_prefs
export DIST_F90_REQUEST="gfortran"
export ACC_PLOT_PACKAGE="pgplot"
export ACC_PLOT_DISPLAY_TYPE="X"
export ACC_ENABLE_OPENMP="$USE_MPI"
export ACC_ENABLE_MPI="$USE_MPI"
export ACC_FORCE_BUILTIN_MPI="Y"
export ACC_ENABLE_GFORTRAN_OPTIMIZATION="Y"
export ACC_ENABLE_SHARED="$SHARED"
export ACC_ENABLE_SHARED_ONLY="$SHARED"
export ACC_ENABLE_FPIC="Y"
export ACC_ENABLE_PROFILING="N"
export ACC_SET_GMAKE_JOBS="$(nproc)"
EOF

echo "## Invoking dist_source_me"

source ./util/dist_source_me

echo -e "\n## Invoking dist_build_production\n"

echo '<details>'
echo '<summary>Build output</summary>'
echo ''
echo '```'

./util/dist_build_production

echo '```'
echo '</details>'
