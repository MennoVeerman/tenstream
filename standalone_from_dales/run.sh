#!/bin/bash

TENSTREAM_SRC=$HOME/tenstream
EXE=${TENSTREAM_SRC}/build/bin/run_cld_fld
PTH=/scratch-shared/<...>/<...>/<...>

RUN="srun -p short -n 240 -t 00:15:00 --constraint=haswell"
CLDFILE=${PTH}/cloud_field.nc
ATMFILE=${PTH}/afglus_100m.dat
OUTDIR=${PTH}

BASE_OPT="-thermal no"

OUT=$OUTDIR/sw_tenstream.nc
OPT=$BASE_OPT

[ ! -e $OUT ] && $RUN $EXE -cld $CLDFILE -atm_filename $ATMFILE -out $OUT $OPT



