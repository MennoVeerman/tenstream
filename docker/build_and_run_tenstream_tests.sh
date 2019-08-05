#!/bin/bash

. $WORKDIR/.profile

mkdir -p $WORKDIR/TenstreamLUT_cpy
cd $WORKDIR/TenstreamLUT_cpy
cp -rns $WORKDIR/TenstreamLUT/* .
echo "-lut_basename $(pwd)/LUT" >> $HOME/.petscrc

cd $WORKDIR/tenstream || exit
rm -rf $WORKDIR/tenstream/build
mkdir build && cd build && cmake .. -DSYST=$SYST -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DNETCDF_DIR=$NETCDF_DIR || exit
make -j all || exit
#bin/createLUT 3_10 || exit
#bin/createLUT 8_10 || exit
#bin/createLUT rectilinear_wedge_5_8 || exit
make check || exit
