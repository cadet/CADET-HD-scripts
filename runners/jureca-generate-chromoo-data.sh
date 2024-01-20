#!/usr/bin/env bash

## NOTE: Run this in a mass sim directory. We assume the corresponding flow directory is named similarly in ../../FLOW/$SIM_DIR
# - Split mesh + data
# - generate pvtus
# - Run paravision at different NRAD and SHELLTYPE

# Execution strategy
# 1. Run this script directly, and let the ${JRUN_COMMAND} calls create jobs
# 2. echo commands (including cds), dump into text file, and use ${JRUN_COMMAND}w to execute as required

# Unfortunately (currently), running these as separate commands would require NRADS * SHELLTYPES * 8 jrun commands
# A more efficient approach would be to run commands in batch on a smaller number of nodes.

BASE="$PWD"
NRADS=(1 2 5 10)
SHELLTYPES=("EQUIDISTANT" "EQUIVOLUME")
JRUN_COMMAND="jrun"
SIM_DIR=$(basename $(realpath .))

## Split mesh and data
# WARNING: make sure there are not other sims with the same mass running this in parallel. If so, ensure that the mesh is split already.
split-mass-mesh-data.sh -m ../mesh -s . -mo ../meshes_split --bulkc --bedc --bedq -fp "${SIM_DIR}_"
## Run mixd2pvtu generated from above command
fd mixd2pvtu.b -x ${JRUN_COMMAND} -ne -n -v -c 'srun -n 48 mixd2pvtu {}'

sleep 3
while [[ -n $(jrun --running) ]]; do
    sleep 10
done

cd "output_bed_c"
TIME="05:00:00"
for NRAD in "${NRADS[@]}"; do 
    for SHELLTYPE in ${SHELLTYPES[@]}; do 
        ${JRUN_COMMAND} -T ${TIME} --title bedc_${NRAD} -v -c "pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} --scale 0.75 -o ${SHELLTYPE}_N${NRAD}_FULL_U.DV_SCALED_0.75" -n -ne
    done 
done
cd "${BASE}"

cd "output_bed_q"
TIME="05:00:00"
for NRAD in "${NRADS[@]}"; do 
    for SHELLTYPE in ${SHELLTYPES[@]}; do 
        ${JRUN_COMMAND} -T ${TIME} --title bedq_${NRAD} -v -c "pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} --scale 0.25 -o ${SHELLTYPE}_N${NRAD}_FULL_U.DV_SCALED_0.25" -n -ne
    done 
done
cd "${BASE}"

# Compartmentalization doesn't make too much sense when we consider MASS (and not flux/chromatograms) as reference data.
# bedc and bedq are literally the same. Only bulkc and bulkq will vary, but even that, only minorly. We need timesteps to be fine as well for this.
# Only that the reduced order model will now be more (geometrically) accurate when we account for the correct packed bed length.
# Nonetheless, we generate results here.

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

cd "output_bulk_c"
TIME="05:00:00"
for NRAD in "${NRADS[@]}"; do 
    for SHELLTYPE in ${SHELLTYPES[@]}; do 
        ${JRUN_COMMAND} -T ${TIME} --title bulk_clip_${NRAD} -v -c "pvrun radial_shell_integrate --project clip twice '${bed_zmin}/${bed_zmax}' 'z#' --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_CLIPPED_U.DV" -n -ne
        ${JRUN_COMMAND} -T ${TIME} --title bulk_full_${NRAD} -v -c "pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_FULL_U.DV" -n -ne
    done 
done
cd "$BASE"

cd "../../FLOW/${SIM_DIR}"

if [[ -f "mixd2pvtu.in" ]]; then
    outpath=$(grep "outpath" mixd2pvtu.in | awk '{print $2}')
    outpath=${outpath:-output}
    if [ -d "$outpath" ]; then
        cd "$outpath"
    else
        jrun -v -c "srun -n 48 mixd2pvtu mixd2pvtu.in" -ne -n
        cd "$outpath"
    fi
else
    echo "No mixd2pvtu.in in FLOW dir."
    exit -1
fi

TIME="00:10:00"
for NRAD in "${NRADS[@]}"; do 
    for SHELLTYPE in ${SHELLTYPES[@]}; do 
        ${JRUN_COMMAND} -T "${TIME}" --title flow_clip_${NRAD} -v -c "pvrun radial_shell_integrate --project clip twice '${bed_zmin}/${bed_zmax}' 'z#' --nrad ${NRAD} --shelltype ${SHELLTYPE} --divide-by-length -o ${SHELLTYPE}_N${NRAD}_CLIPPED_FLOWRATES" -n -ne
        ${JRUN_COMMAND} -T "${TIME}" --title flow_full_${NRAD} -v -c "pvrun radial_shell_integrate --nrad ${NRAD} --shelltype ${SHELLTYPE} --divide-by-length -o ${SHELLTYPE}_N${NRAD}_FULL_FLOWRATES" -n -ne
        ${JRUN_COMMAND} -T "${TIME}" --title por_clip_${NRAD} -v -c "pvrun radial_porosity --project clip twice '${bed_zmin}/${bed_zmax}' 'z#' --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_CLIPPED" -n -ne
        ${JRUN_COMMAND} -T "${TIME}" --title por_full_${NRAD} -v -c "pvrun radial_porosity --nrad ${NRAD} --shelltype ${SHELLTYPE} -o ${SHELLTYPE}_N${NRAD}_FULL" -n -ne
    done
done
cd "${BASE}"
