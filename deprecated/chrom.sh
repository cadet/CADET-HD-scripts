#!/bin/bash

set -e

## Author: Jayghosh Rao
## Description: Chromatography mesh -> job submission
## Built for JURECA/Slurm

## ====================================================
# Be aware that this script may not be entirely safe.
# It wasn't entirely tested or sanitized. 
## ====================================================

# chrom.sh (initially meshprep.sh) 
# Takes msh2 as input and broadly performs the following actions:
# 1. Convert gmsh to mixd
# 2. Create two meshes: FLOW and MASS
# 3. Partition meshes
# 4. Submit jobs (with MASS depending on FLOW)

# TODO: [ ] Make checkstage work with set -e
# TODO: [ ] auto mixd2xdmf after simulations?
# TODO: [ ] auto chrompost after simulations?
# TODO: [x] Check and show squeue on mem512 before confirming decompose on backend
# TODO: [x] Ask for backend (devel/mem512/batch)
# TODO: [x] Velocity Control 
# TODO: [x] test templates
# TODO: [x] load modules
# TODO: [x] check for sub.sh
# TODO: [x] allow different startpoints for code exec
# TODO: [x] modularize script
# TODO: [ ] interactive/auto modes
# TODO: [x] Vipe the input files
# TODO: [x] Cleanup mode
# TODO: [ ] Log output


if [ -z "$1" ]; then
    echo "No mesh provided!"
    echo "Usage: chrom.sh <mesh.msh2>"
    exit
fi

TOOLS=(gmsh2mixd shiftmixd scalemixd rmmat gennmat decompose.metis sub.sh)

exitFlag=0;

for TOOL in ${TOOLS[@]}; do 
    if [ -x "$(command -v $TOOL )" ]; then
        echo "Found $TOOL!"
    else 
        echo "ERR: $TOOL doesn't exist."
        exitFlag=1;
    fi
done

[[ -d ~/templates ]] || { echo "No templates found."; exitFlag=1; }

if ((exitFlag == 1)); then
    echo "Please install prerequisites!"
    exit
fi

# nparts=240
shiftvals="-5 -5 -707.5"
scalevals="1e-4 1e-4 1e-4"
decompBackend="no"
backend="mem512"
MSendtime="4000"
MSdtexp="2.0 + heaviside(t-100.0)*(t/100.0 - 1.0)"
FSrngdexp="1  3  2.0*2.09e-4*( 1.0 - (x*x + y*y)/( 5.01e-4 * 5.01e-4 ) )"

job_name=${PWD##*/}
job_cpus=240
job_acct="jibg12"
job_part="batch"
job_flow_time="02:00:00"
job_mass_time="12:00:00"
xns_path="~/bin/xns-mpi"

read -p "Enter Translate Offsets (dx dy dz): " -e -i "$shiftvals" shiftvals
read -p "Enter Scale Factors (dx dy dz): " -e -i "$scalevals" scalevals
# read -p "Enter number of partitions: " -e -i "$nparts" nparts
read -p "Decompose in backend?: " -e -i "$decompBackend" decompBackend


if [[ "$decompBackend" == "yes" ]]; then 
    squeue -p mem512
    read -p "Backend: " -e -i "$backend" backend
fi

if [ -z "$backend" ]; then
    decompBackend = "no"
fi

read -p "job_name: " -e -i "$job_name" job_name
read -p "job_cpus: " -e -i "$job_cpus" job_cpus
read -p "job_acct: " -e -i "$job_acct" job_acct
read -p "job_part: " -e -i "$job_part" job_part
read -p "job_flow_time: " -e -i "$job_flow_time" job_flow_time
read -p "job_mass_time: " -e -i "$job_mass_time" job_mass_time
read -p "xns_path: " -e -i "$xns_path" xns_path

echo "Setting up directories..."
echo "========================="
mkdir -p FLOW/mesh MASS/mesh
cp ~/templates/decompose.mass.in MASS/mesh
cp ~/templates/decompose.flow.in FLOW/mesh
cp ~/templates/xns.mass.in MASS/
cp ~/templates/xns.flow.in FLOW/
cp ~/templates/mixd2xdmf.mass.in MASS/
cp ~/templates/mixd2xdmf.flow.in FLOW/
cp ~/templates/job.mass.sh MASS
cp ~/templates/job.flow.sh FLOW
echo 


$EDITOR FLOW/xns.flow.in MASS/xns.mass.in

function check_stage(){
    module list 2>&1 | grep "Stages/2018b" &>/dev/null && return 0 || return 1
}
check_stage
correct_stage=$?

if [ correct_stage ]; then 
    echo "Setting environment (2018a)"
    echo "==========================="
    module --force purge
    module use /usr/local/software/jureca/OtherStages
    module load Stages/2018a
    module load intel-para Boost flex
    echo 
fi

# read -p "FLOW/sim rngdexp: " -e -i "$FSrngdexp" FSrngdexp
# read -p "MASS/sim endtime: " -e -i "$MSendtime" MSendtime
# read -p "MASS/sim dtexp: " -e -i "$MSdtexp" MSdtexp

function convert_to_mixd(){
    echo "Converting msh to mixd..."
    echo "========================="
    gmsh2mixd -m sd -d 4 "$1"

    transform_mixd
}

function transform_mixd(){
    if [ ! -z "$shiftvals" ]; then
        echo "Shifting mesh..."
        shiftmixd -tet $shiftvals
    fi
    if [ ! -z "$scalevals" ]; then
        echo "Scaling mesh..."
        scalemixd -tet $scalevals
    fi

    gen_spacetime_mesh
}

function gen_spacetime_mesh(){

    echo "Creating space-time mesh..."
    echo "==========================="
    cp mxyz mxyz.space
    cp mtbl mtbl.space
    cat mxyz.space >> mxyz
    cat mtbl.space >> mtbl
    NN=$(tail -n 1 minf | awk '{ print $2 }')
    let NN_ST="$NN * 2"
    sed -i 's/nn/#nn/' minf 
    echo "nn     $NN_ST" >> minf
    echo 

    mv mien mxyz mrng minf mtbl mmat MASS/mesh
    # dir_setup
    remove_beads
}


function remove_beads(){

    echo "Removing bead material..."
    echo "========================="
    cd MASS/mesh
    rmmat -tet -st ../../FLOW/mesh 2
    cd - 
    echo 

    partition_flow

}

function partition_flow(){

    echo "Decomposing FLOW mesh..."
    echo "========================"
    cd FLOW/mesh

    if [[ "$decompBackend" == "no" ]]; then  
        gendual -e tet && genneim --nen 4 && decompose.metis -e tet -w -p $job_cpus
    elif [[ "$decompBackend" == "yes" ]]; then
        gendual -e tet && genneim --nen 4 && srun -N1 -n24 -A jibg12 -p $backend decompose.metis -e tet -w -p $job_cpus
    else 
        echo "== Not Partitioning Mesh =="
    fi

    cd - 
    echo 

    partition_mass
}

function partition_mass(){
    echo "Decomposing MASS mesh..."
    echo "========================"
    cd MASS/mesh

    if [[ "$decompBackend" == "no" ]]; then  
        gendual -e tet && genneim --nen 4 -s &&  decompose.metis -e tet -w -p $job_cpus -s
    elif [[ "$decompBackend" == "yes" ]]; then
        gendual -e tet && genneim --nen 4 -s && srun -N1 -n24 -A jibg12 -p $backend decompose.metis -e tet -w -p $job_cpus -s
    else 
        echo "== Not Partitioning Mesh =="
    fi

    echo "Generating nmat file..."
    echo "======================="
    gennmat 4 st
    cd - 
    echo 

    flow_setup

}

function flow_setup() {

    echo "Setting up FLOW/sim directory..." 
    echo "================================" 

    mkdir FLOW/sim
    cp FLOW/xns.flow.in FLOW/sim/xns.in
    cp FLOW/job.flow.sh FLOW/sim/job.sh
    cd FLOW/sim

    sed -i "s/<job_name>/$job_name/g" job.sh
    sed -i "s/<job_cpus>/$job_cpus/g" job.sh
    sed -i "s/<job_acct>/$job_acct/g" job.sh
    sed -i "s/<job_part>/$job_part/g" job.sh
    sed -i "s/<job_flow_time>/$job_flow_time/g" job.sh
    sed -i "s:<xns_path>:$xns_path:g" job.sh


    # sed -i "s/^rngdexp 0/rngdexp $FSrngdexp/" xns.in

    # FJOUT=$(sub.sh -r)
    # FJID=$(echo $FJOUT | awk '{print $NF}')
    # echo "$FJOUT"
    cd -

    mass_setup
}

function mass_setup() {
    
    echo "Setting up MASS/sim directory..." 
    echo "================================" 

    mkdir MASS/sim
    cp MASS/xns.mass.in MASS/sim/xns.in
    cp MASS/job.mass.sh MASS/sim/job.sh
    cd MASS/sim

    sed -i "s/<job_name>/$job_name/g" job.sh
    sed -i "s/<job_cpus>/$job_cpus/g" job.sh
    sed -i "s/<job_acct>/$job_acct/g" job.sh
    sed -i "s/<job_part>/$job_part/g" job.sh
    sed -i "s/<job_mass_time>/$job_mass_time/g" job.sh
    sed -i "s:<xns_path>:$xns_path:g" job.sh
    sed -i "/^#SBATCH -J/i #SBATCH -d afterok:$FJID" job.sh

    # sed -i "s/^endtime 0/endtime $MSendtime/" xns.in
    # sed -i "s:^dtexp 0:dtexp $MSdtexp:" xns.in
    # sub.sh -r
    cd - 

    echo "DONE!"
    echo "====="

}

function full_clean() {
    rm -rf FLOW MASS mxyz* mtbl*
}


POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -pf)
            partition_flow
            shift
            ;;
        -pm)
            partition_mass
            shift
            ;;
        -fs)
            flow_setup
            shift
            ;;
        -ms)
            mass_setup
            shift
            ;;
        --fclean)
            full_clean
            shift
            ;;
        -f)
            convert_to_mixd $2
            shift
            shift
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            convert_to_mixd $1
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters
