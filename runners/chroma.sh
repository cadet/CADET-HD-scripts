#!/usr/bin/env bash

# The chroma script to handle simulation workflow was unnecessarily complicated
# Since my experience in bash is now better than a few years ago, here's a simpler 
# script that is better suited to the job. It should also make more sense to readers.

# Usage: ./chroma.sh -n 1280 -ddp "srun -A jibg12 -N1 -n48 -p dc-cpu-devel" mesh_column.msh2

filter_integer() {
    if [[ $1 =~ ^[[:digit:]]+$ ]]; then
        echo "$1"
    fi
}

function proclaim()
{
    echo "################################################################################"
    center_and_fold "$@"
    echo "################################################################################"
}

function center_and_fold() {
    local string="$1"
    local length=${#string}
    local block_size=80

    if (( length >= block_size )); then
        # Fold the string to fit within 80 characters
        fold -s -w "$block_size" <<< "$string" | sed -e :a -e '/^\n*$/{$d;N;ba' -e '}'
    else
        # Calculate padding for centering
        local pad=$(( (block_size - length) / 2 ))
        printf "%*s%s%*s\n" $pad '' "$string" $pad ''
    fi
}

function die(){
    echo -e "ERROR: $@" >&2
    exit -1
}

function setup_dir()
{
    [ -d "$1" ] && rm -rf "$1"
    mkdir -p "$1"
}

function spacetimeify()
{
    if [[ -f "$1" ]]; then
        echo "spacetimeifying $1"
        cp "$1" "${1}.space"
        cat "${1}.space" >> "$1"
    fi
}

function handle_mass_periodicity()
{
    ## Some code in XNS: genmesh.F / genmtbl()/genmprd() etc allows us
    ## to just duplicate the semidiscrete data to generate the spacetime data.
    ## So we just generate mtbl and mprd for semidiscrete meshes and then double it.
    ## To be "proper" we'd have to add an offset (nnspace) to the second half of the mtbl/mprd data
    echo "Handling mass periodicity"
    if [[ "$PERIODICITY" == "XY" ]]; then
        genmprd 0 0 x -readmtbl
        genmprd 0 0 y -readmprd -readmtbl
    elif [[ "$PERIODICITY" == "XYZ" ]]; then
        genmprd 0 0 x -readmtbl
        genmprd 0 0 y -readmprd -readmtbl
        genmprd 0 0 z -readmprd -readmtbl
    elif [[ "$PERIODICITY" == "Z" ]]; then
        genmprd 0 0 z -readmtbl
    fi
}

function handle_flow_periodicity()
{
    if [[ "$PERIODICITY" == "XY" ]]; then
        genmprd 0 0 x
        genmprd 0 0 y -readmprd
    elif [[ "$PERIODICITY" == "XYZ" ]]; then
        genmprd 0 0 x
        genmprd 0 0 y -readmprd
        genmprd 0 0 z -readmprd
    elif [[ "$PERIODICITY" == "Z" ]]; then
        genmprd 0 0 z
    fi
}

# Globals
PERIODICITY=
NMESHPARTS=
ETYPE=tet && NEN=4 && MESH_ORDER=1
PARTICLES_SURFACE_GROUP_IDX=4
PARTICLES_VOLUME_GROUP_IDX=2
DECOMPOSE_DISPATCH_PREFIX=
ROOT_DIR=$PWD
MASS_MESH_DIR=MASS/mesh
FLOW_MESH_DIR=FLOW/mesh
MFILES=('mien' 'mxyz' 'mrng' 'minf' 'mtbl' 'mmat' 'mprd' 'mtbl.space' 'mxyz.space' 'mtbl.dual' 'minf.space' 'mprd.space')
LEGACY_GMSH2MIXD=0
DECOMPOSE_COMMAND="decompose.metis"

## Alternative options. See below.
# ETYPE=tetP2 && NEN=10 && MESH_ORDER=2 
# DECOMPOSE_DISPATCH_PREFIX="srun -N1 -n48 -A jibg12 -p dc-cpu-devel"

## Commandline args processing
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -tetP2|-o2)
            ETYPE=tetP2 && NEN=10 && MESH_ORDER=2
            shift # past value
            ;;
        -ddd|--decompose-default-dispatch)
            DECOMPOSE_DISPATCH_PREFIX="srun -N1 -n48 -A jibg12 -p dc-cpu-devel"
            DECOMPOSE_COMMAND="$DECOMPOSE_DISPATCH_PREFIX $DECOMPOSE_COMMAND"
            shift # past value
            ;;
        -ddp|--decompose-dispatch-prefix)
            DECOMPOSE_DISPATCH_PREFIX="$2"
            [ -n "$2" ] && DECOMPOSE_COMMAND="$DECOMPOSE_DISPATCH_PREFIX $DECOMPOSE_COMMAND"
            shift # past value
            shift # past value
            ;;
        -n|--nmeshparts)
            NMESHPARTS="$2"
            shift # past value
            shift # past value
            ;;
        -l|--legacy-gmsh2mixd)
            LEGACY_GMSH2MIXD=1
            shift
            ;;
        -p|--periodicity)
            PERIODICITY="$2"
            shift
            shift
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

[ -n "$1" ] && MESHFILE="$1" || die "No mesh file given"

if [[ $(hostname) =~ "jureca" ]]; then 
    export MODULEPATH=/p/project/cjibg12/modulefiles:/p/software/jurecadc/supercomputer_modules:/p/software/jurecadc/productionstages
    module load Stages/2022
    module load GCC/11.2.0 ParaStationMPI/5.5.0-1 ParMETIS/4.0.3-double Boost/1.78.0 VTK/9.1.0 flex/2.6.4  
fi

## Cleanup
for mfile in "${MFILES[@]}"; do
    rm -f "$mfile"
done

setup_dir "$FLOW_MESH_DIR"
setup_dir "$MASS_MESH_DIR"

proclaim "Starting mesh conversion"

if [[ "$LEGACY_GMSH2MIXD" == 1 ]]; then
    gmsh2mixd_3D -d "$PARTICLES_SURFACE_GROUP_IDX" -m sd "$MESHFILE"
else
    gmsh2mixdv2 -d "$PARTICLES_SURFACE_GROUP_IDX" -o "$MESH_ORDER" "$MESHFILE"
fi

handle_mass_periodicity

proclaim "Creating spacetime meshes"
spacetimeify mxyz
spacetimeify mtbl
spacetimeify mprd

# update minf
cp minf minf.space
NN=$(awk '/^nn/{print $2}' minf)
NNST=$(( NN * 2 ))
sed -i 's/nn/# nn' minf
echo "nn     $NNST" >> minf

proclaim "Moving mesh files"

## move to MASS/mesh
for mfile in "${MFILES[@]}"; do
    echo "moving $mfile"
    mv "$mfile" "$MASS_MESH_DIR"
done

proclaim "Removing packed-bed to create interstitial mesh"

cd "$MASS_MESH_DIR"
rmmat -"$ETYPE" -st ../../FLOW/mesh "$PARTICLES_VOLUME_GROUP_IDX"
cd $ROOT_DIR

cd "$FLOW_MESH_DIR"
handle_flow_periodicity
cd $ROOT_DIR

cd "$FLOW_MESH_DIR"
if [[ -n $(filter_integer "$NMESHPARTS") ]]; then
    proclaim "Partitioning FLOW mesh"
    gendual -e "$ETYPE" && genneim --nen "$NEN"
    $DECOMPOSE_COMMAND -e "$ETYPE" -w -p "$NMESHPARTS"
fi
cd "$ROOT_DIR"

cd "$MASS_MESH_DIR"
if [[ -n $(filter_integer "$NMESHPARTS") ]]; then
    proclaim "Partitioning MASS mesh"
    gendual -e "$ETYPE" && genneim --nen "$NEN" -s
    $DECOMPOSE_COMMAND -e "$ETYPE" -w -p "$NMESHPARTS" -s
fi
cd "$ROOT_DIR"
