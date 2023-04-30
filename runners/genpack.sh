#!/usr/bin/bash

## setup the packing generation
## Requires PackingGeneration.exe in path
## $1 is the title of the directory
## $2 may be the diameters.txt file 

## This is a preliminary script to just use the 3 algorithms one by one.
## To get a target porosity, generate a first packing with fba
## Try to keep the Calc porosity within 0.4-0.6. (Using contraction rate)
## The average bead radius is controlled by the 'particles count'. 
## (diameters.txt only serves to provide a ratio of bead diamters)
## Then run ls with random seeds at a required contraction rate to get the porosity you desire.
## see optimizer.sh

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

putGenmeshFile(){
cat > genmesh.in <<EOFMARKER
# dryRun 1

## Packing
packing         	packing_fixed.xyzd
packingPrecision	8
translateOffsets	0 0 0
preScalingFactor	1

# Periodicity
periodic xyz
periodicOffsets auto
periodicInlet 2
periodicOutlet 2

## Bead selection range
nBeads $nbeads

## Bead modifications
reduced	0.9997

## Container Geometry
autoContainment	0
containerShape	1 
box 0 0 0 $boxx $boxy $boxz

## Mesh size
lc_beads 	0.2
lc_out   	0.2
lc_bridge	0.2

## GMSH settings
General.NumThreads 0
Mesh.MaxNumThreads 0
Mesh.ScalingFactor 0.0001
Mesh.Generate      2
Mesh.Algorithm     5
Mesh.Algorithm3D   10
EOFMARKER
}

title="$1"
diameters="$2"

contractionrate=1e-5

[ -z $1 ] && echo "usage: <title> <diameters.txt>" && exit

read -p "#particles: " -e -i "$nbeads" nbeads
read -p "box size [x]: " -e -i "$boxx" boxx
read -p "box size [y]: " -e -i "$boxy" boxy
read -p "box size [z]: " -e -i "$boxz" boxz
read -p "contraction rate: " -e -i "$contractionrate" contractionrate

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
dumpy --dpacking packing.xyzd --nfo packing.nfo -w packing_fixed.xyzd
cd $ROOT

cd $title/mesh
cp ../lsgd/packing_fixed.xyzd .
cp ../lsgd/packing.nfo .
putGenmeshFile
cd $ROOT
