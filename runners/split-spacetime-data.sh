#!/usr/bin/env bash

# Script to split 2-domain space-time mesh+data 
# into separate bulk and bed mesh+data sets

function die(){
    echo -e "ERROR: $@" >&2
    exit -1
}

ROOT="$PWD"

MASS_MESH_DIR=
MASS_SIM_DIR=
OUTPUT_MESH_ROOT="splitmeshes"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -m|--mesh-dir)
            MASS_MESH_DIR="$2"
            shift; shift; exit
            ;;
        -s|--sim-dir)
            MASS_SIM_DIR="$2"
            shift; shift; exit
            ;;
        -o|--output-mesh-dir)
            OUTPUT_MESH_ROOT="$2"
            shift; shift; exit
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

[[ -d "$MASS_SIM_DIR" ]] || die "MASS_SIM_DIR either not set, or doesn't exist."
[[ -d "$MASS_MESH_DIR" ]] || die "MASS_MESH_DIR either not set, or doesn't exist."

echo cd "$MASS_MESH_DIR"
echo rmmat -tet -st "${OUTPUT_MESH_ROOT}/mesh-bulk" 2 # Remove packed bed region
echo cp nmap nmap.bulk
echo rmmat -tet -st mesh-bed 1 # Remove interstitial region
echo cp nmap nmap.bed
echo cd "$ROOT"

echo cd "$MASS_SIM_DIR"
echo mixdsplit -m "${MASS_MESH_DIR}/minf" -N "$MASS_MESH_DIR/nmap.bulk" -o interstitial_c.all
echo mixdsplit -m "${MASS_MESH_DIR}/minf" -N "$MASS_MESH_DIR/nmap.bed" -o bed_c.all -i 0
echo mixdsplit -m "${MASS_MESH_DIR}/minf" -N "$MASS_MESH_DIR/nmap.bed" -o bed_q.all -i 1
echo cd "$ROOT"
