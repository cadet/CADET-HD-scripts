#!/usr/bin/env bash

## A "script" to be sourced by the jobscript on JURECA to make life easier
## Assumes that commands are run from the sim directory of the current column.
## Assumes directory structure:
##          - ROOT
##              - Column1
##                  - FLOW
##                      - mesh
##                      - sim
##                  - MASS
##                      - mesh
##                      - sim
##              - Column2
##                  ...


NEN=4
NSD=3

## mixd double records. Find nn from mxyz, using: mdrec <mxyz> <nsd>
get_mixd_double_nrec(){
        echo $(( $(stat --printf="%s" $1) / ( ${2:-1} * 8 ) ))
}


## mixd integer records. Find ne from mien, using: mdrec <mien> <nen>
get_mixd_int_nrec(){
    echo $(( $(stat --printf="%s" $1) / ( ${2:-1} * 4 ) ))
}

get_nts(){
    NN=$(get_mixd_double_nrec ../mesh/mxyz $NSD)
    NTS=$(get_mixd_double_nrec $1 $(( 2*NN )) )
    echo $NTS
}

mapflow_wrapper(){
    # usage: mapflow sim
    CWD=${PWD}
    cd ../mesh
    mapflow -tet ../../FLOW/mesh ../../FLOW/$1
    cd "$CWD"
    cp ../mesh/flowfield .
}

stitch_flow(){
    ## Usage: stitch_flow <column> <sim-dirname>
    ## Usage: stitch_flow inlet sim
    ROOT=${PWD}
    cd ../../../$1/FLOW/$2
    stitchperiodic data.out -r 2 -t 2 -n 4
    cd $ROOT
    cp ../../../$1/FLOW/$2/rng.data rng.data.in
    cp ../../../$1/FLOW/$2/rng.xyz rng.xyz.in
    NN_RNG=$(get_mixd_double_nrec rng.xyz.in $NSD)
    sed -i "/nn_rng/cnn_rng $NN_RNG" xns.in
}

stitch_mass(){
    ## Usage: stitch_mass <column> <sim-dirname>
    ## Usage: stitch_mass inlet sim
    ROOT=${PWD}
    cd ../../../$1/MASS/$2
    NTS=$(get_nts data.all)
    stitchperiodic data.all -r 2 -t $NTS -n 2 -s
    cd $ROOT
    cp ../../../$1/MASS/$2/rng.data rng.data.in
    cp ../../../$1/MASS/$2/rng.xyz rng.xyz.in
    NN_RNG=$(get_mixd_double_nrec rng.xyz.in $NSD)
    sed -i "/nn_rng/cnn_rng $NN_RNG" xns.in
}
