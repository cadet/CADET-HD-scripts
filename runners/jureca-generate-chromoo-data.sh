#!/usr/bin/env bash

set -e

function die(){
    echo -e "ERROR: $@" >&2
    exit -1
}

function ensure_dirs()
{
    for ARG in "$@"; do
        [[ -d "$ARG" ]] || die "Dir not found: $ARG"
    done
}

function get_clip_positions_from_packinfo()
{
    ## Get bed start and end positions to clip geometry
    cd ../..
    if [[ ! -f "input.yaml_packing_processed.json" ]]; then
        module load gmsh/4.11.0-2ac03e-copymesh.lua
        source ~/cjibg12/miniconda3/bin/activate dev
        pack-info input.yaml
        source ~/cjibg12/miniconda3/bin/deactivate
        module unload gmsh/4.11.0-2ac03e-copymesh.lua
    fi
    bed_zmin=$(jq '.post_scale_data.bed_zmin' input.yaml_packing_processed.json | sed 's/-/ -/') ## -ve values should be escaped for the shell
    bed_zmax=$(jq '.post_scale_data.bed_zmax' input.yaml_packing_processed.json)
    cd "$BASE"
}

## NOTE: Run this in a mass sim directory. We assume the corresponding flow directory is named similarly in ../../FLOW/$SIM_DIR
# - Split mesh + data
# - generate pvtus
# - Run paravision at different NRAD and SHELLTYPE

# Execution strategy
# 1. Run this script directly, and let the ${JRUN_COMMAND} calls create jobs, OR
# 2. echo commands (including cds), dump into text file, and use ${JRUN_COMMAND}w to execute as required

# Unfortunately (currently), running these as separate commands would require NRADS * SHELLTYPES * 8 jrun commands
# A more efficient approach would be to run commands in batch on a smaller number of nodes.

BASE="$PWD"
NRADS=(1 2 5 10)
# SHELLTYPES=("EQUIDISTANT" "EQUIVOLUME")
SHELLTYPES=("EQUIDISTANT")
SIM_DIR=$(basename $(realpath .))

FLOW_SIM_DIR="../../FLOW/${SIM_DIR}"
FLOW_OUTPUT_DIR="output"
BULK_OUTPUT_DIR="bulk_c"
BEDC_OUTPUT_DIR="bed_c"
BEDQ_OUTPUT_DIR="bed_q"

# Remember to add a space before -ve values so that they aren't mistaken for the start of an arg flag
bed_zmin=
bed_zmax=

RUN_FULL=false

## Default jrun command, no chunking of files for parallel runs
JRUN_COMMAND=jrun
OUTPUT_SUFFIX=

# Example case to allow chunking files for parallel jrun execs. 
# Use this with high NRAD values, else we hit our 24 hour time limit
# NOTE: Remove echos in actual run.
# chunk_files.sh -e '*.pvtu' -n 10 -k 3 -c echo jrun -v -ne -n -C pvrun -np 48 radial_shell_integrate --nrad 5 --shelltype EQUIDISTANT -o OUTPUT_{ICHUNK}_FULL_U.DV 
JRUN_COMMAND="chunk_files.sh -e '*.pvtu' -n 10 -k 3 -C echo jrun"
OUTPUT_SUFFIX="_{ICHUNK}"

# ## Split mesh and data
# # WARNING: make sure there are not other sims with the same mass running this in parallel. If so, ensure that the mesh is split already.
# split-mass-mesh-data.sh -m ../mesh -s . -mo ../meshes_split --bulkc --bedc --bedq -fp "${SIM_DIR}_"
# ## Run mixd2pvtu generated from above command
# fd mixd2pvtu.b -x jrun -ne -n -v -c 'srun -n 48 mixd2pvtu {}'

## TODO: This is bad. FIXME
# sleep 3
# while [[ -n $(jrun --running) ]]; do
#     sleep 10
# done

if ensure_dirs "${BEDC_OUTPUT_DIR}" ; then
    cd "${BEDC_OUTPUT_DIR}" || exit
    TIME="24:00:00"
    for NRAD in "${NRADS[@]}"; do
        for SHELLTYPE in ${SHELLTYPES[@]}; do
            if [[ "$RUN_FULL" == "true" ]] ; then
                ${JRUN_COMMAND} -T ${TIME} --title bedc_${NRAD} -v -n -ne -C pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} --scale 0.75 -o ${SHELLTYPE}_N${NRAD}_FULL_U.DV_SCALED_0.75${OUTPUT_SUFFIX}
            fi
        done
    done
    cd "${BASE}"
fi

if ensure_dirs "${BEDQ_OUTPUT_DIR}" ; then
    cd "${BEDQ_OUTPUT_DIR}" || exit
    TIME="24:00:00"
    for NRAD in "${NRADS[@]}"; do
        for SHELLTYPE in ${SHELLTYPES[@]}; do
            if [[ "$RUN_FULL" == "true" ]] ; then
                ${JRUN_COMMAND} -T ${TIME} --title bedq_${NRAD} -v -n -ne -C pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} --scale 0.25 -o ${SHELLTYPE}_N${NRAD}_FULL_U.DV_SCALED_0.25${OUTPUT_SUFFIX}
            fi
        done
    done
    cd "${BASE}"
fi


if ensure_dirs "${BULK_OUTPUT_DIR}" ; then
    cd "${BULK_OUTPUT_DIR}" || exit
    TIME="24:00:00"
    for NRAD in "${NRADS[@]}"; do
        for SHELLTYPE in ${SHELLTYPES[@]}; do
            if [ -n "${bed_zmin}" ] && [ -n "${bed_zmax}" ]; then
                ${JRUN_COMMAND} -T ${TIME} --title bulk_clip_${NRAD} -v -n -ne -C pvrun radial_shell_integrate --project clip twice "'${bed_zmin}/${bed_zmax}'" 'z#' --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_CLIPPED_U.DV${OUTPUT_SUFFIX}
            fi

            if [[ "$RUN_FULL" == "true" ]] ; then
                ${JRUN_COMMAND} -T ${TIME} --title bulk_full_${NRAD} -v -n -ne -C pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_FULL_U.DV${OUTPUT_SUFFIX}
            fi
        done
    done
    cd "$BASE"
fi

if ensure_dirs "${FLOW_SIM_DIR}" ; then
    # Generate flow results if not already done
    cd "${FLOW_SIM_DIR}" || exit
    if [[ -f "mixd2pvtu.in" ]]; then
        outpath=$(grep "outpath" mixd2pvtu.in | awk '{print $2}')
        outpath=${outpath:-output}
        if [ -d "$outpath" ]; then
            cd "$outpath"
        else
            jrun -v -n -ne -C "srun -n 48 mixd2pvtu mixd2pvtu.in"
            ## TODO: Wait for jrun to finish
            cd "$outpath"
        fi
    else
        echo "No mixd2pvtu.in in FLOW dir."
        exit -1
    fi
    cd ${BASE}
fi

if ensure_dirs "${FLOW_OUTPUT_DIR}" ; then
    cd "${FLOW_OUTPUT_DIR}" || exit
    TIME="24:00:00"
    for NRAD in "${NRADS[@]}"; do
        for SHELLTYPE in ${SHELLTYPES[@]}; do
            if [ -n "${bed_zmin}" ] && [ -n "${bed_zmax}" ]; then
                ${JRUN_COMMAND} -T ${TIME} --title flow_clip_${NRAD} -v -n -ne -C pvrun radial_shell_integrate --project clip twice "'${bed_zmin}/${bed_zmax}'" 'z#' --nrad ${NRAD} --shelltype ${SHELLTYPE} --divide-by-length -o ${SHELLTYPE}_N${NRAD}_CLIPPED_FLOWRATES${OUTPUT_SUFFIX}
                ${JRUN_COMMAND} -T ${TIME} --title por_clip_${NRAD} -v -n -ne -C pvrun radial_porosity --project clip twice "'${bed_zmin}/${bed_zmax}'" 'z#' --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_CLIPPED${OUTPUT_SUFFIX}
            fi
            if [[ "$RUN_FULL" == "true" ]] ; then
                ${JRUN_COMMAND} -T ${TIME} --title flow_full_${NRAD} -v -n -ne -C pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} --divide-by-length -o ${SHELLTYPE}_N${NRAD}_FULL_FLOWRATES${OUTPUT_SUFFIX}
                ${JRUN_COMMAND} -T ${TIME} --title por_full_${NRAD} -v -n -ne -C pvrun radial_porosity --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_FULL${OUTPUT_SUFFIX}
            fi
        done
    done
    cd "${BASE}"
fi
