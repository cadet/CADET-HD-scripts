#!/usr/bin/bash

################
# Script to run a periodic inlet/outlet batch
################

set -euo pipefail

# RUN="PARALLEL"              ## SERIAL/PARALLEL. UNSET to not run xns
# HOST="REMOTE"                ## if LOCAL, mpiexec
# NODES=1
# CORES=32
# XNS="xns-mpid-per"
# MODULES="intel-para Boost flex"

RUN=""              ## SERIAL/PARALLEL. UNSET to not run xns
HOST="LOCAL"                ## if LOCAL, mpiexec
NODES=1
CORES=32
XNS="xns-mpid"
MODULES=""

get_ndouble(){
    echo $(( $(stat --printf="%s" $1) / ( ${2:-1} * 8 ) ))
}

## SETUP
MESHFILEPREFIX="$1"
ROOT="$PWD"

mkdir -p inlet outlet column
cp "$MESHFILEPREFIX"_column.msh2 column
cp "$MESHFILEPREFIX"_inlet.msh2 inlet
cp "$MESHFILEPREFIX"_outlet.msh2 outlet

chroma -c 
sed -i '/periodic/s/no/xy/' chroma.json
sed -i '/decompose/s/no/yes/' chroma.json
sed -i "/job-cpu/c\ \ \ \ \"job-cpus\": \"$CORES\"," chroma.json
cp chroma.json inlet
cp chroma.json outlet
cp chroma.json column

cp ~/templates/xns.{flow,mass}.in inlet
cp ~/templates/xns.stitch.flow.in column/xns.flow.in
cp ~/templates/xns.stitch.mass.in column/xns.mass.in
cp ~/templates/xns.stitch.flow.in outlet/xns.flow.in
cp ~/templates/xns.stitch.mass.in outlet/xns.mass.in

cd $ROOT/inlet
chroma "$MESHFILEPREFIX"_inlet.msh2

cd $ROOT/column
chroma "$MESHFILEPREFIX"_column.msh2

cd $ROOT/outlet
chroma "$MESHFILEPREFIX"_outlet.msh2

runSim(){
    [ -n "$MODULES" ] && module load $MODULES
    [[ $RUN == "SERIAL" ]] && $XNS
    if [[ $HOST == "LOCAL" ]]; then
        [[ $RUN == "PARALLEL" ]] && mpiexec -np $CORES $XNS
    else
        [[ $RUN == "PARALLEL" ]] && srun -N $NODES -n $CORES -A jibg12 -p dc-cpu-devel $XNS
    fi
}

## SIM INLET
cd $ROOT/inlet/FLOW/sim
runSim

cd $ROOT/inlet/MASS/mesh
mapflow -tet ../../FLOW/mesh ../../FLOW/sim
cp flowfield $ROOT/inlet/MASS/sim

cd $ROOT/inlet/MASS/sim
runSim

## STITCH
cd $ROOT
dumpy --stitch inlet column

## SIM COLUMN
cd $ROOT/column/FLOW/sim
nn_rng=$(get_ndouble rng.xyz.in 3)
sed -i "/nn_rng/cnn_rng $nn_rng" xns.in
runSim

cd $ROOT/column/MASS/mesh
mapflow -tet ../../FLOW/mesh ../../FLOW/sim
cp flowfield $ROOT/column/MASS/sim

cd $ROOT/column/MASS/sim
nn_rng=$(get_ndouble rng.xyz.in 3)
sed -i "/nn_rng/cnn_rng $nn_rng" xns.in
runSim

cd $ROOT
dumpy --stitch column outlet

## SIM OUTLET
cd $ROOT/outlet/FLOW/sim
nn_rng=$(get_ndouble rng.xyz.in 3)
sed -i "/nn_rng/cnn_rng $nn_rng" xns.in
runSim

cd $ROOT/outlet/MASS/mesh
mapflow -tet ../../FLOW/mesh ../../FLOW/sim
cp flowfield $ROOT/outlet/MASS/sim

cd $ROOT/outlet/MASS/sim
nn_rng=$(get_ndouble rng.xyz.in 3)
sed -i "/nn_rng/cnn_rng $nn_rng" xns.in
runSim
