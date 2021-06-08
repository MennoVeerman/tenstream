#!/bin/bash
BIN=/home/mveerman/tenstream/build/bin/run_cld_fld
#make -j -C $(dirname $BIN)/.. $(basename $BIN) 

RUN="srun -p short -n 240 -t 00:15:00 --constraint=haswell"
TENSTREAM_SRC=$HOME/tenstream
CLDFILE=/scratch-shared/mveerman/dales_20160815/rrtmgp_100m/cloud_field.nc
ATMFILE=${TENSTREAM_SRC}/examples/mytest/afglus_100m.dat
OUTDIR=/scratch-shared/mveerman/tenstream_offline

BASE_OPT="-thermal no"
OUT=$OUTDIR/output_sw.nc
OPT=$BASE_OPT
echo $RUN $BIN -cld $CLDFILE -atm_filename $ATMFILE -out $OUT $OPT
[ ! -e $OUT ] && $RUN $BIN -cld $CLDFILE -atm_filename $ATMFILE -out $OUT $OPT

echo $OUT
echo $OPT



