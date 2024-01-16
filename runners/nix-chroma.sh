#!/usr/bin/env bash

## Use the following shebang to run the script under a nix develop shell
#!/usr/bin/env -S nix develop --impure /home/jayghoshter/dev/tools/pymesh/flake.nix --command zsh

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

function run_simulation_stage()
{
    if [[ "$DISPATCH" == "JURECA" ]]; then
        ensure_match "(jureca|jrc.*)" $(hostname)
    fi

    local SIM_STAGE="$1"
    local SIM_DIR="$2"

    proclaim "Running simulation $SIM_STAGE/$SIM_DIR"

    [ -n "$SIM_STAGE" ] && [ -n "$SIM_DIR" ] || die "Bad sim stage or dir: $SIM_STAGE $SIM_DIR"
    [ -n "$NMESHPARTS" ] && [ -n "$REMOTE" ] || die "Bad globals NMESHPARTS: $NMESHPARTS or REMOTE: $REMOTE"

    local SIM_STAGE_UPPER=$(echo "$SIM_STAGE" | tr '[:lower:]' '[:upper:]')
    local SIM_STAGE_LOWER=$(echo "$SIM_STAGE" | tr '[:upper:]' '[:lower:]')

    ensure_files "xns.${SIM_STAGE_LOWER}.in"
    ensure_dirs "$SIM_STAGE_UPPER"/{,mesh}

    mkdir -p "$SIM_STAGE_UPPER/$SIM_DIR"
    if ! check_files "$SIM_STAGE_UPPER/$SIM_DIR/xns.in"; then
        check_files "xns.${SIM_STAGE_LOWER}.in" && cp "xns.${SIM_STAGE_LOWER}.in" "$SIM_STAGE_UPPER/$SIM_DIR/xns.in"
    fi

    ensure_files "$SIM_STAGE_UPPER/$SIM_DIR/xns.in"
    [ -f "mixd2pvtu.${SIM_STAGE_LOWER}.in" ] && cp "mixd2pvtu.${SIM_STAGE_LOWER}.in" "$SIM_STAGE_UPPER/$SIM_DIR/mixd2pvtu.in"

    [[ "$SIM_STAGE_UPPER" == "MASS" ]] && mapflow_wrapper "$SIM_DIR"

    if [[ "$DISPATCH" == "REMOTE" ]] ; then 
        echo "Pushing files to $REMOTE"
        ensure_run mirror -m -y push $REMOTE -f "$SIM_STAGE_UPPER/mesh" "$SIM_STAGE_UPPER/$SIM_DIR"
        echo "Submitting job $SIM_STAGE_UPPER/$SIM_DIR on $REMOTE"
        cd "$SIM_STAGE_UPPER/$SIM_DIR"
        local JRUN_OUT=$(ensure_run mirror cmd 'jrun -x -nt '$NTPN' -np '$NMESHPARTS' -n -ne' --target $REMOTE)
        JOB_ID=$(echo "$JRUN_OUT" | tail -n 1 | grep Submitted | awk '{print $2}')
        echo "$SIM_STAGE_UPPER simulation dispatched with JobID: $JOB_ID"
        cd "$BASE"
    elif [[ "$DISPATCH" == "JURECA" ]] ; then
        echo "Submitting job $SIM_STAGE_UPPER/$SIM_DIR"
        cd "$SIM_STAGE_UPPER/$SIM_DIR"
        local JRUN_OUT=$(jrun -x -nt $NTPN -np $NMESHPARTS -n -ne)
        JOB_ID=$(echo "$JRUN_OUT" | tail -n 1 | grep Submitted | awk '{print $2}')
        echo "$SIM_STAGE_UPPER simulation dispatched with JobID: $JOB_ID"
        cd "$BASE"
    elif [[ "$DISPATCH" == "LOCAL" ]]; then
        echo "Running job $SIM_STAGE_UPPER/$SIM_DIR"
        cd "$SIM_STAGE_UPPER/$SIM_DIR"
        ensure_run $XNS_LOCAL_COMMAND < xns.in
        cd "$BASE"
    fi
}

function wait_for_results()
{
    if [[ "$DISPATCH" == "JURECA" ]]; then
        ensure_match "(jureca|jrc.*)" $(hostname)
    elif [[ "$DISPATCH" == "LOCAL" ]]; then
        return
    fi

    local SIM_STAGE="$1"
    local SIM_DIR="$2"

    [ -n "$SIM_STAGE" ] && [ -n "$SIM_DIR" ] || die "Bad sim stage or dir: $SIM_STAGE $SIM_DIR"
    local SIM_STAGE_UPPER=$(echo "$SIM_STAGE" | tr '[:lower:]' '[:upper:]')
    local SIM_STAGE_LOWER=$(echo "$SIM_STAGE" | tr '[:upper:]' '[:lower:]')

    [[ -n "$JOB_ID" ]] || die "Please provide a JOB_ID to wait for"
    proclaim "Waiting for simulation"

    cd "$SIM_STAGE_UPPER/$SIM_DIR"
    
    if [[ "$DISPATCH" == "REMOTE" ]] ; then 
        mirror pull $REMOTE -y
        while mirror cmd "jrun --is-queued $JOB_ID" --target "$REMOTE"
        do
            echo "Job $JOB_ID is still queued. Waiting $SIM_COMPLETION_CHECK_TIME seconds before checking again."
            sleep $SIM_COMPLETION_CHECK_TIME
            mirror pull $REMOTE -y
        done
    elif [[ "$DISPATCH" == "JURECA" ]] ; then 
        while jrun --is-queued "$JOB_ID"
        do
            echo "Job $JOB_ID is still queued. Waiting $SIM_COMPLETION_CHECK_TIME seconds before checking again."
            sleep $SIM_COMPLETION_CHECK_TIME
        done
    fi

    cd "$BASE"
}

function convert_results()
{
    local SIM_STAGE="$1"
    local SIM_DIR="$2"

    local SIM_STAGE_UPPER=$(echo "$SIM_STAGE" | tr '[:lower:]' '[:upper:]')
    local SIM_STAGE_LOWER=$(echo "$SIM_STAGE" | tr '[:upper:]' '[:lower:]')

    proclaim "Converting results"
    [ -f "$SIM_STAGE_UPPER/$SIM_DIR/mixd2pvtu.in" ] || die "No mixd2pvtu file found in dir: $SIM_STAGE_UPPER/$SIM_DIR"

    cd "$SIM_STAGE_UPPER/$SIM_DIR"
    ensure_run $MIXD2PVTU_COMMAND mixd2pvtu.in
    cd "$BASE"
}

function mapflow_wrapper()
{
    local SIM_DIR="${1:-sim}"

    [[ "$PWD" == "$BASE" ]] || die "mapflow_wrapper must be run from base: $BASE"

    ensure_files "FLOW/mesh/"{mxyz,mien,mrng,minf}
    ensure_files "FLOW/$SIM_DIR/data.out"
    ensure_files "MASS/mesh/"{mxyz,mien,mrng,minf}
    ensure_dirs "MASS/mesh" "MASS/$SIM_DIR"

    cd "MASS/mesh"
    mapflow -tet ../../FLOW/mesh ../../FLOW/$SIM_DIR
    cp flowfield "../$SIM_DIR"
    cd "$BASE"
}

function generate_mesh()
{
    if [[ ! -f mesh_column.msh2 ]]; then
        proclaim "Generating mesh"
        ensure_run mesh input.yaml
    fi


}

# TODO: Consider moving chroma.sh code into this function fully. That would
# reduce the overhead of maintaining 2 scripts and passing arguments and
# parameters between them. Error code handling also becomes easier In order to
# just prepare the mesh locally, I could then just run this script with the
# appropriate flags.
function prepare_mesh()
{
    [[ -f mesh_column.msh2 ]] || die "No such file: mesh_column.msh2"
    if ! check_files FLOW/mesh/{mxyz,mien,mrng,minf} MASS/mesh/{mxyz,mien,mrng,minf} ; then
        # TODO: Ensure chroma.sh is error-code compliant. 
        # TODO: Ensure that all component tools are error-code compliant
        if [[ -n "$DISPATCH_PREFIX" ]]; then
            ensure_run chroma.sh mesh_column.msh2 -n $NMESHPARTS -l --dispatch-prefix "$DISPATCH_PREFIX"
        else
            ensure_run chroma.sh mesh_column.msh2 -n $NMESHPARTS -l
        fi
    fi
}

function setup_dir()
{
    [[ -n "$1" ]] || die "Empty arg: $1"
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

function spacetimeify_mesh()
{
    echo "Creating spacetime meshes"
    spacetimeify mxyz
    spacetimeify mtbl
    spacetimeify mprd

    # update minf
    cp minf minf.space
    NN=$(awk '/^nn/{print $2}' minf)
    NNST=$(( NN * 2 ))
    sed -i 's/nn/# nn/' minf
    echo "nn     $NNST" >> minf
}

function handle_mass_periodicity()
{
    ## Some code in XNS: genmesh.F / genmtbl()/genmprd() etc allows us
    ## to just duplicate the semidiscrete data to generate the spacetime data.
    ## So we just generate mtbl and mprd for semidiscrete meshes and then double it.
    ## To be "proper" we'd have to add an offset (nnspace) to the second half of the mtbl/mprd data
    ensure_commands genmprd
    if [[ "$PERIODICITY" == "XY" ]]; then
        echo "Handling MASS periodicity: $PERIODICITY"
        genmprd 0 0 x -readmtbl
        genmprd 0 0 y -readmprd -readmtbl
    elif [[ "$PERIODICITY" == "XYZ" ]]; then
        echo "Handling MASS periodicity: $PERIODICITY"
        genmprd 0 0 x -readmtbl
        genmprd 0 0 y -readmprd -readmtbl
        genmprd 0 0 z -readmprd -readmtbl
    elif [[ "$PERIODICITY" == "Z" ]]; then
        echo "Handling MASS periodicity: $PERIODICITY"
        genmprd 0 0 z -readmtbl
    fi
}

function handle_flow_periodicity()
{
    ensure_commands genmprd
    if [[ "$PERIODICITY" == "XY" ]]; then
        echo "Handling FLOW periodicity: $PERIODICITY"
        genmprd 0 0 x
        genmprd 0 0 y -readmprd
    elif [[ "$PERIODICITY" == "XYZ" ]]; then
        echo "Handling FLOW periodicity: $PERIODICITY"
        genmprd 0 0 x
        genmprd 0 0 y -readmprd
        genmprd 0 0 z -readmprd
    elif [[ "$PERIODICITY" == "Z" ]]; then
        echo "Handling FLOW periodicity: $PERIODICITY"
        genmprd 0 0 z
    fi
}

function prepare_mesh_myself()
{
    local MFILES=('mien' 'mxyz' 'mrng' 'minf' 'mtbl' 'mmat' 'mprd' 'mtbl.space' 'mxyz.space' 'mtbl.dual' 'minf.space' 'mprd.space')
    local FLOW_MESH_DIR="FLOW/mesh"
    local MASS_MESH_DIR="MASS/mesh"
    local PARTICLES_SURFACE_GROUP_IDX=4
    local PARTICLES_VOLUME_GROUP_IDX=2
    local MESHFILE="mesh_column.msh2"

    [[ -f "$MESHFILE" ]] || die "No such file: $MESHFILE"
    check_files FLOW/mesh/{mxyz,mien,mrng,minf} MASS/mesh/{mxyz,mien,mrng,minf} && echo "Mesh files exist. Skipping preparation." && return

    proclaim "Preparing mesh"
    echo "Cleaning up files"
    for mfile in "${MFILES[@]}"; do
        rm -f "$mfile"
    done

    echo "Setting up directories"
    setup_dir "$FLOW_MESH_DIR"
    setup_dir "$MASS_MESH_DIR"

    proclaim "Starting mesh conversion"
    if [[ "$LEGACY_GMSH2MIXD" == 1 ]]; then
        ensure_run gmsh2mixd_3D -d "$PARTICLES_SURFACE_GROUP_IDX" -m sd "$MESHFILE"
    else
        ensure_run gmsh2mixdv2 -d "$PARTICLES_SURFACE_GROUP_IDX" -o "$MESH_ORDER" "$MESHFILE"
    fi
    proclaim "Mesh conversion complete"

    handle_mass_periodicity
    spacetimeify_mesh

    echo "Moving mesh files"
    for mfile in "${MFILES[@]}"; do
        echo "moving $mfile"
        mv "$mfile" "$MASS_MESH_DIR"
    done

    echo "Removing packed-bed to create interstitial mesh"
    cd "$MASS_MESH_DIR"
    rmmat -"$ETYPE" -st ../../FLOW/mesh "$PARTICLES_VOLUME_GROUP_IDX"
    gennmat "$NEN" st
    cd "$BASE"

    cd "$FLOW_MESH_DIR"
    handle_flow_periodicity
    cd "$BASE"
}

function decompose_mesh()
{
    local SIM_STAGE="$1"
    local SIM_DIR="$2"

    local SIM_STAGE_UPPER=$(echo "$SIM_STAGE" | tr '[:lower:]' '[:upper:]')
    local SIM_STAGE_LOWER=$(echo "$SIM_STAGE" | tr '[:upper:]' '[:lower:]')

    NMP_04=$(printf "%04d\n" "$NMESHPARTS") 
    local SPACETIME_ARG=
    [[ "$SIM_STAGE_UPPER" == "MASS" ]] && SPACETIME_ARG="-s"

    cd "${SIM_STAGE_UPPER}/mesh"

    if ! check_files neim ; then
        ensure_run gendual -e "$ETYPE" 
        ensure_run genneim --nen "$NEN" "$SPACETIME_ARG"
    fi

    if ! check_files {mprm,nprm}.${NMP_04} ; then
        ensure_run $DECOMPOSE_COMMAND -e "$ETYPE" -p "$NMESHPARTS" -w $SPACETIME_ARG
    fi

    cd "$BASE"
}

function driver()
{
    ## This function uses globals only

    if [[ "$MODE" == "MESH" ]]; then
        generate_mesh
    elif [[ "$MODE" == "PREPARE" ]]; then
        generate_mesh
        prepare_mesh_myself
    elif [[ "$MODE" == "RUN" ]]; then
        generate_mesh
        prepare_mesh_myself
        local SIM_STAGE="${SIM_STAGES[0]}"
        echo "WARNING: Due to FLOW -> MESH dependency, RUN mode can only dispatch one simulation stage at a time."
        echo "WARNING: Running $SIM_STAGE simulation."
        decompose_mesh "$SIM_STAGE"
        run_simulation_stage "$SIM_STAGE" "$SIM_NAME"
    elif [[ "$MODE" == "RUNWAIT" ]]; then
        generate_mesh
        prepare_mesh_myself
        for SIM_STAGE in ${SIM_STAGES[@]}; do 
            decompose_mesh "$SIM_STAGE"
            run_simulation_stage "$SIM_STAGE" "$SIM_NAME"
            ## Necessary to wait for FLOW stage results
            wait_for_results "$SIM_STAGE" "$SIM_NAME"
            convert_results "$SIM_STAGE" "$SIM_NAME"
        done
    elif [[ "$MODE" == "WAIT" ]]; then
        for SIM_STAGE in ${SIM_STAGES[@]}; do 
            wait_for_results "$SIM_STAGE" "$SIM_NAME"
            convert_results "$SIM_STAGE" "$SIM_NAME"
        done
    elif [[ "$MODE" == "CONVERT" ]]; then
        for SIM_STAGE in ${SIM_STAGES[@]}; do 
            convert_results "$SIM_STAGE" "$SIM_NAME"
        done
    fi

    # cd "$SIM_STAGE_UPPER/$SIM_DIR/output"
    # cd "$BASE"
}

## Globals
BASE="$PWD"

## Number of tasks per node to use for XNS job on JURECA
NTPN=128

## Number of mesh parts, or number of processes to run
NMESHPARTS=48

## Remote name. Used primarily by `mirror` and `ssh` to sync data and send commands
REMOTE="jureca"

## Sleep time for checks while waiting for simulation
SIM_COMPLETION_CHECK_TIME=60

## Name of the simulation to be performed. Uses this as directory name.
SIM_NAME="sim"

## Simulation stages to perform
SIM_STAGES=(FLOW MASS)

## MODE can be one of [ MESH, PREPARE, RUN, WAIT, RUNWAIT, CONVERT ]
MODE="RUNWAIT"

# DISPATCH can be one of [ JURECA , REMOTE, LOCAL ]. 
# LOCAL => Run fully locally on a machine/node 
# REMOTE => Run script locally, sync data to/from remote. Use jrun to submit data on remote. Currently tightly coupled with JURECA.
# JURECA => Run script locally on JURECA. Submit job with jrun.
DISPATCH="JURECA" 
JOB_ID=

# Sometimes, tasks need to be dispatched in parallel on compute nodes on JURECA.
# For example, decompose and mixd2pvtu.
# This prefix can specify dispatch parameters, e.g.,
# DISPATCH_PREFIX='srun -N1 -n64 -A jibg12 -p dc-cpu-devel --time=02:00:00'
# We currently
DISPATCH_PREFIX=

## Commands with potential for MPI dispatched runs
DECOMPOSE_COMMAND="decompose.metis"
MIXD2PVTU_COMMAND="mixd2pvtu"
XNS_LOCAL_COMMAND=/home/jayghoshter/local/simulation/bin/xns-mpi

ETYPE=tet && NEN=4 && MESH_ORDER=1
PERIODICITY=
LEGACY_GMSH2MIXD=0

SOFTWARE_STAGE=2022


## Commandline args processing
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -n|--nmeshparts)
            NMESHPARTS="$2"
            shift # past value
            shift # past value
            ;;
        -nt|--ntpn)
            NTPN="$2"
            shift # past value
            shift # past value
            ;;
        -s|--sim)
            SIM_NAME="$2"
            shift # past value
            shift # past value
            ;;
        -m|--mesh)
            MODE="MESH"
            shift
            ;;
        -p|--prepare)
            MODE="PREPARE"
            shift
            ;;
        -c|--convert)
            MODE="CONVERT"
            shift
            ;;
        -d|--dispatch)
            DISPATCH=$(echo "$2" | tr '[:lower:]' '[:upper:]')
            ensure_match "^(JURECA|REMOTE|LOCAL)$" "$DISPATCH"
            shift
            shift
            ;;
        -r|--run|--no-wait)
            MODE="RUN"
            shift
            ;;
        -w|--wait)
            MODE="WAIT"
            JOB_ID=$(filter_integer "$2")
            [[ -n "$JOB_ID" ]] || die "Bad JOB_ID provided to --wait"
            shift
            shift
            ;;
        -ddp|--default-dispatch-prefix)
            DISPATCH_PREFIX="srun -A jibg12 -N1 -n48"
            shift # past value
            ;;
        -dp|--dispatch-prefix)
            DISPATCH_PREFIX="$2"
            shift # past value
            shift # past value
            ;;
        -l|--legacy-gmsh2mixd)
            LEGACY_GMSH2MIXD=1
            shift
            ;;
        --periodicity)
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

if [[ -n "$1" ]]; then
    ensure_match "^(flow|mass|FLOW|MASS)$" "$1"
    SIM_STAGES=("$1")
fi 

ensure_match "^(MESH|PREPARE|RUN|WAIT|RUNWAIT|CONVERT)$" "$MODE"

echo "Dispatch mode: $DISPATCH"
echo "Dispatch prefix: $DISPATCH_PREFIX"
echo "Remote: $REMOTE"

if [[ $(hostname) =~ (jureca|jrc.*) ]]; then 
    source /p/software/jurecadc/lmod/8.4.1/init/zsh
    export MODULEPATH=/p/project/cjibg12/modulefiles:/p/software/jurecadc/supercomputer_modules:/p/software/jurecadc/productionstages

    if [[ "$SOFTWARE_STAGE" == 2022 ]]; then
        ## Unfortunately, on JURECA, ParMETIS-double is only available on 2022 stage
        ## So we use this as the default
        module load Stages/2022 GCC/11.2.0 ParaStationMPI/5.5.0-1 ParMETIS/4.0.3-double Boost/1.78.0 VTK/9.1.0 flex/2.6.4  

    elif [[ "$SOFTWARE_STAGE" == 2024 ]]; then
        module load Stages/2024 GCC/12.3.0 ParaStationMPI/5.9.2-1 ParMETIS/4.0.3 Boost/1.82.0 VTK/9.3.0 flex/2.6.4
    fi

    # for pymesh
    module load gmsh/4.11.0-2ac03e-copymesh.lua
    source ~/cjibg12/miniconda3/bin/activate dev 
fi

if [[ -n "$DISPATCH_PREFIX" ]]; then
    DECOMPOSE_COMMAND="$DISPATCH_PREFIX $DECOMPOSE_COMMAND"
    MIXD2PVTU_COMMAND="$DISPATCH_PREFIX $MIXD2PVTU_COMMAND"
    XNS_LOCAL_COMMAND="$DISPATCH_PREFIX $XNS_LOCAL_COMMAND"
fi

driver
