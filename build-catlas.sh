#! /bin/bash
set -x
set -e

mkdir tr-cross
mv tr-cross.{gxt,mxt} tr-cross/
./spacegraphcats/build-catlas.py tr-cross 5

mkdir acido-short
mv acido-short.{gxt,mxt} acido-short/
./spacegraphcats/build-catlas.py acido-short 5
