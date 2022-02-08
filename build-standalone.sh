#!/bin/sh
set -o errexit
DIR=`mktemp -d`
trap "rm -rf $DIR" 0 2

VER=empi-`git describe --tags --exact-match`-lin64
mkdir -p $DIR/$VER
cp README.md LICENCE $DIR/$VER
cmake -DSTANDALONE=TRUE -DWITH_CUDA=FALSE -S. -B$DIR
( cd $DIR \
  && make -j8 empi \
  && mv empi $VER/empi-lin64 \
  && zip -9r $VER.zip $VER/
)
mv $DIR/$VER.zip .
