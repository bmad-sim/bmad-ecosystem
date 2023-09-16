#!/bin/bash

# Downloads the latest open_spacecharge_mod.f90 open_spacecharge_core_mod.f90
# from https://github.com/RobertRyne/OpenSpaceCharge
#
# Only open_spacecharge_mod.f90 is actually needed for Bmad
#

VERSION=1.0.1
PACKAGE=OpenSpaceCharge

F=v$VERSION.tar.gz
DIR=$PACKAGE-$VERSION
RELEASE=https://github.com/RobertRyne/$PACKAGE/archive/$F

echo "Getting $RELEASE"

wget $RELEASE
tar -xf $F
rm $F

for file in open_spacecharge_mod.f90 open_spacecharge_core_mod.f90
do
 echo "cp $DIR/code/$file $ACC_ROOT_DIR/bmad/space_charge/$file"
 cp $DIR/code/$file $ACC_ROOT_DIR/bmad/space_charge/$file
done

rm -r $DIR
