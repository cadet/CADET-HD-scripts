#!/usr/bin/env bash

# Script to extract (separately) the concentration data in the particle and bulk domains
# Output data can be used with paravision to create mass curves in different domains.
# Mass curves can be used as reference data in chromoo fits.

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

function ensure_run(){
    [ -x $(command -v "$1") ] || die "no such command: $1"
    "$@"
    [[ $? == 0 ]] || die "command exited with error: $@"
}

function ensure_commands()
{
    for ARG in "$@"; do
        [[ -x $(command -v "$ARG") ]] || die "no such command: $1"
    done
}

function ensure_files()
{
    for ARG in "$@"; do
        [[ -f "$ARG" ]] || die "File not found: $ARG"
    done
}

function ensure_match()
{
    CMP_STR="$1"
    shift
    for ARG in "$@"; do
        [[ "$ARG" =~ $CMP_STR ]] || die "$ARG doesn't match $CMP_STR"
    done
}

function check_files()
{
    for ARG in "$@"; do
        [[ -f "$ARG" ]] || return -1
    done
    return 0
}

function ensure_dirs()
{
    for ARG in "$@"; do
        [[ -d "$ARG" ]] || die "Dir not found: $ARG"
    done
}

function split_mesh()
{
    cd "$MASS_MESH_DIR"
    if ! check_files "${OUTPUT_MESH_ROOT}/mesh-bulk/"{minf,mxyz,mien,mrng}; then
        ensure_run rmmat -tet -st "${OUTPUT_MESH_ROOT}/mesh-bulk" 2 # Remove packed bed region
        ensure_run cp nmap nmap.bulk
    fi
    if ! check_files "${OUTPUT_MESH_ROOT}/mesh-bed/"{minf,mxyz,mien,mrng}; then
        ensure_run rmmat -tet -st "${OUTPUT_MESH_ROOT}/mesh-bed" 1 # Remove interstitial region
        ensure_run cp nmap nmap.bed
    fi
    cd "$ROOT"
}

function split_data()
{
    DATA_FILE="$1"
    DATA_DIR="$(dirname $DATA_FILE)"

    if [[ -n "$OUTPUT_DATA_ROOT" ]]; then
        OUTPUT_PREFIX="$(realpath $OUTPUT_DATA_ROOT)"
    else
        OUTPUT_PREFIX="$(realpath $DATA_DIR)"
    fi

    mkdir -p "$OUTPUT_PREFIX"
    ensure_dirs "$OUTPUT_PREFIX"

    cd "$DATA_DIR" || die "Bad cd: $DATA_DIR" 
    if [[ "$SPLIT_BULKC" == 1 ]] && [[ ! -f "${OUTPUT_PREFIX}/bulk_c.all" ]]; then 
        ensure_run mixdsplit -m "${MASS_MESH_DIR}/minf" -N "$MASS_MESH_DIR/nmap.bulk" -o "${OUTPUT_PREFIX}"/bulk_c.all & 
    fi
    if [[ "$SPLIT_BEDQ" == 1 ]] && [[ ! -f "${OUTPUT_PREFIX}/bed_c.all" ]]; then
        ensure_run mixdsplit -m "${MASS_MESH_DIR}/minf" -N "$MASS_MESH_DIR/nmap.bed" -o "${OUTPUT_PREFIX}"/bed_c.all -i 0 & 
    fi
    if [[ "$SPLIT_BEDC" == 1 ]] && [[ ! -f "${OUTPUT_PREFIX}/bed_q.all" ]]; then
        ensure_run mixdsplit -m "${MASS_MESH_DIR}/minf" -N "$MASS_MESH_DIR/nmap.bed" -o "${OUTPUT_PREFIX}"/bed_q.all -i 1 &
    fi
    cd "$ROOT"

    wait
}

ROOT="$PWD"

MASS_MESH_DIR=
MASS_SIM_DIR=
OUTPUT_MESH_ROOT="meshes_split"
OUTPUT_DATA_ROOT=

DATA_FILES=$(find . -type f -iname data.all | sort)

SPLIT_BULKC=
SPLIT_BEDC=
SPLIT_BEDQ=

SCRIPT_PATH=$(readlink -f -- "$BASH_SOURCE")

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -m|--mesh-dir)
            MASS_MESH_DIR="$2"
            shift; shift;
            ;;
        -s|--sim-dir)
            MASS_SIM_DIR="$2"
            DATA_FILES=$(find "$MASS_SIM_DIR" -type f -iname data.all | sort)
            shift; shift;
            ;;
        -mo|--output-mesh-dir)
            OUTPUT_MESH_ROOT="$2"
            shift; shift;
            ;;
        -do|--output-data-dir)
            OUTPUT_DATA_ROOT="$2"
            shift; shift;
            ;;
        -d|--data)
            shift;
            DATA_FILES=()
            N_ARGS=0
            for ARG in "$@"; do 
                [[ "$ARG" =~ ^- ]] && break
                DATA_FILES+=("$1")
                shift;
            done
            ;;
        --bulkc)
            SPLIT_BULKC=1
            shift;
            ;;
        --bedc)
            SPLIT_BEDC=1
            shift;
            ;;
        --bedq)
            SPLIT_BEDQ=1
            shift;
            ;;
        -h|--help)
            echo "The following args are processed:" 
            grep -Po '\s*-\w+\|[^)]*' "$SCRIPT_PATH"
            exit
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

[[ -d "$MASS_MESH_DIR" ]] || die "MASS_MESH_DIR either not set, or doesn't exist."
[[ "${#DATA_FILES[@]}" != 0 ]] || die "No data files found. Use -s or -d args."

ensure_commands rmmat mixdsplit 

mkdir -p "${OUTPUT_MESH_ROOT}/mesh-bulk" "${OUTPUT_MESH_ROOT}/mesh-bed"

proclaim "Splitting mesh"
split_mesh

proclaim "Splitting data"
for DATA_FILE in "${DATA_FILES[@]}"; do 
    split_data "$DATA_FILE" &
done

read -r -d '' CONFIG_BULKC << EOFMARKER
title bulk_c
outpath bulk_c

minf ${OUTPUT_MESH_ROOT}/mesh-bulk/minf
mxyz ${OUTPUT_MESH_ROOT}/mesh-bulk/mxyz
mien ${OUTPUT_MESH_ROOT}/mesh-bulk/mien

elemtype tet
nrec 99999
data $(echo ${DATA_FILES[@]} | sed 's|data.all|bulk_c.all|')
ndf 1
# dt 50
# dtFile ..
EOFMARKER

echo -e "$CONFIG_BULKC" > mixd2pvtu.bulk_c.in
echo -e "$CONFIG_BULKC" | sed 's|bulk|bed|' > mixd2pvtu.bed_c.in
echo -e "$CONFIG_BULKC" | sed 's|_c|_q|;s|bulk|bed|' > mixd2pvtu.bed_q.in

wait
