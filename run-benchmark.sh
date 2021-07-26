#!/bin/sh
S=2000000
N=256
LOG="`git rev-parse --short HEAD`.log"
rm -f $LOG
while [ $N -lt $S ] ; do
  printf "%7d " $N | tee -a $LOG
  ./benchmark $S $N | tee -a $LOG
  N=`expr 2 \* $N`
done
