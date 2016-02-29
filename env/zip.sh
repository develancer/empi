#!/bin/sh
# script (and the whole directory) for release packing

if [ ! -f zip.sh ] ; then
  echo "script must be run from its own directory"
  exit 1
fi

if [ $# -lt 1 ] ; then
  echo "USAGE: $0 version"
  exit 1
fi

VER="$1"
for SRC in * ; do
  if [ -d $SRC ] ; then
    DIR="empi-$VER-$SRC"
    ZIP="$DIR.zip"
    ln -s $SRC $DIR
    rm -f $ZIP
    zip -9r $ZIP $DIR
    rm $DIR
  fi
done
