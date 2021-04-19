#!/usr/bin/bash

## setup the packing generation
## Requires PackingGeneration.exe in path
## $1 is the title of the directory
## $2 may be the diameters.txt file 


set -eo pipefail

putConfigFile(){
cat > generation.conf <<EOFMARKER
Particles count: $1
Packing size: $2 $3 $4
Generation start: $5
Seed: $6
Steps to write: 1000
Boundaries mode: 1
Contraction rate: 1e-005
Generation mode: 1

1. boundaries mode: 1 - bulk; 2 - ellipse (inscribed in XYZ box, Z is length of an ellipse); 3 - rectangle
2. generationMode = 1 (Poisson, R) or 2 (Poisson in cells, S)
EOFMARKER
}

title="$1"
diameters="$2"

[ -z $1 ] && echo "usage: <title> <diameters.txt>" && exit

read -p "#particles: " -e -i "$nbeads" nbeads
read -p "box size [x]: " -e -i "$boxx" boxx
read -p "box size [y]: " -e -i "$boxy" boxy
read -p "box size [z]: " -e -i "$boxz" boxz

ROOT=${PWD}
mkdir -p $title/{fba,ls,lsgd,mesh}

if [[ -n $diameters ]]; then
    cp $diameters $title/fba/diameters.txt
    cp $diameters $title/ls/diameters.txt
    cp $diameters $title/lsgd/diameters.txt
fi

cd $title/fba
putConfigFile $nbeads $boxx $boxy $boxz 1 $RANDOM
PackingGeneration.exe -fba |& tee out.log
cd $ROOT

cd $title/ls
cp ../fba/packing.xyzd .
putConfigFile $nbeads $boxx $boxy $boxz 0 $RANDOM
PackingGeneration.exe -ls |& tee out.log
cd $ROOT

cd $title/lsgd
cp ../ls/packing.xyzd .
putConfigFile $nbeads $boxx $boxy $boxz 0 $RANDOM
PackingGeneration.exe -lsgd |& tee out.log
cd $ROOT
