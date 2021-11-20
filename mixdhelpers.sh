#!/usr/bin/env bash

## Meant to be a rewrite/replacement to jobtools.sh
## Help get nrec from mixd mesh and data files easily

# USAGE: Source this file, then
# EXAMPLE: get_nrec -m ../mesh data.all -n 2

DOUBLESIZE=8
INTSIZE=4
NSD=3

function die(){
    echo "ERROR: $@"
    return
}

function get_filesize(){
    [[ -f "$1" ]] && stat --printf="%s" $1 || die "$1: Not a file!"
}

function get_nrec() {
    MESHDIR=$PWD
    NDF=1

    ALL_ARGS="$@"

    POSITIONAL=()
    while [[ $# -gt 0 ]]
    do
        key="$1"

        case $key in
            -m|--meshdir)
                MESHDIR="$2"
                shift # past value
                shift # past value
                ;;
            -n|--ndf)
                NDF="$2"
                shift # past value
                shift # past value
                ;;
            -tet|--tet)
                NEN=4
                shift # past value
                ;;
            -tri|--tri)
                NEN=3
                shift # past value
                ;;
            --nsd)
                NSD="$2"
                shift # past value
                shift # past value
                ;;
            *)    # unknown option
                POSITIONAL+=("$1") # save it in an array for later
                shift # past argument
                ;;
        esac
    done
    set -- "${POSITIONAL[@]}" # restore positional parameters

    minf=$(find "$MESHDIR" -iname "*minf*" -type f | head -n 1)
    mien=$(find "$MESHDIR" -iname "*mien*" -type f | head -n 1)
    mxyz=$(find "$MESHDIR" -iname "*mxyz*" -type f | head -n 1)

    DATAFILE="$1"

    FILESIZE=$(get_filesize "$DATAFILE")

    NN=$(get_nn $ALL_ARGS)

    NREC=$(( $FILESIZE / ( $NN * $NDF * $DOUBLESIZE) ))
    echo "$NREC"

}

function get_nn(){

    MESHDIR=$PWD
    ALL_ARGS="$@"

    POSITIONAL=()
    while [[ $# -gt 0 ]]
    do
        key="$1"

        case $key in
            -m|--meshdir)
                MESHDIR="$2"
                shift # past value
                shift # past value
                ;;
            -tet|--tet)
                NEN=4
                shift # past value
                ;;
            --nsd)
                NSD="$2"
                shift # past value
                shift # past value
                ;;
            *)    # unknown option
                POSITIONAL+=("$1") # save it in an array for later
                shift # past argument
                ;;
        esac
    done
    set -- "${POSITIONAL[@]}" # restore positional parameters

    mxyz=$(find "$MESHDIR" -iname "*mxyz*" -type f | head -n 1)

    NN=$(( $(get_filesize "$mxyz") / ( $NSD * $DOUBLESIZE ) ))
    echo "$NN"

}

findup() {
    # USAGE: findup -iname "*minf*" -type f 
    # USAGE: findup -iname "mesh" -type d 

    # findpath="$1"
    # shift 1

    findpath="$PWD"

    while [[ $findpath != "/" ]];
    do
        out=$(find "$findpath" "$@" | head -n 1)
        [ -n "$out" ] && echo "$out" && return

        # Note: if you want to ignore symlinks, use "$(realpath -s "$findpath"/..)"
        findpath="$(readlink -f "$findpath"/..)"
    done

}

mapflow_v2(){
    # USAGE: mapflow_v2
    # DETAILS: 
    #   - Run from flow-simname/mass-simname/ dir.
    #   - Finds closest "mesh" dir above
    #   - Assumes the following dir structure
    #   - Moves flowfield to velocity dir.
    # DIRS: 
    #   - mesh
    #       - mass
    #       - flow
    #   - flow-simname
    #       - <flow data>
    #       - mass
    #           - <mass data>

    CWD=${PWD}
    MESHDIR_MASS="$(findup -type d -iname mesh)/mass"
    MESHDIR_FLOW="$(findup -type d -iname mesh)/flow"
    cd "$MESHDIR_MASS"
    mapflow -tet "$MESHDIR_FLOW" $(readlink -f "$CWD/..")
    cd "$CWD"
    mv "$MESHDIR_MASS/flowfield" ..
}
