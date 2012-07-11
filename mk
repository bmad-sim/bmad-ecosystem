#!/bin/bash

#/nfs/acc/libs/Linux_x86_64_intel/extra/bin/cmake -DCMAKE_BUILD_TYPE=None .; make -j2

THIS_SCRIPT=`basename $0`

# Determine if the current working directory is a buildable project
if ( [ ! -e CMakeLists.txt ] ) then
  echo "This working directory is not a project supported by the ACC build system."
  exit 1
fi

if [[ ${THIS_SCRIPT} == "mk" ]]
then
  BUILD_TYPE="production"
fi

if [[ ${THIS_SCRIPT} == "mkd" ]]
then
  BUILD_TYPE="debug"
fi


# If the local out-of-source build directory does not yet exist, create it.
if ( [ ! -d ${BUILD_TYPE} ] ) then
  mkdir ${BUILD_TYPE}
fi

cd ${BUILD_TYPE}
#cmake -DCMAKE_BUILD_TYPE=None ..
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ..
make -j2 $@

