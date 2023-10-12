#!/bin/sh
set -o errexit
DIR=`mktemp -d`
trap "rm -rf $DIR" 0 2

OS=`uname -s | tr '[A-Z]' '[a-z]'`
if [ "$OS" = darwin ] ; then
  OS=macos
  STANDALONE=FALSE
else
  STANDALONE=TRUE
fi

VER=empi-`git describe --tags --exact-match`-$OS
mkdir -p $DIR/$VER
cp README.md LICENCE $DIR/$VER
cmake -DSTANDALONE=$STANDALONE -DWITH_CUDA=FALSE -S. -B$DIR
( cd $DIR \
  && make -j8 empi \
  && mv empi $VER/empi-$OS \
  && zip -9r $VER.zip $VER/
)
mv $DIR/$VER.zip .
