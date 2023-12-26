#!/usr/bin/env zsh

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
    local SIM_STAGE="$1"
    local SIM_DIR="$2"

    proclaim "Running simulation $SIM_STAGE/$SIM_DIR"

    [ -n "$SIM_STAGE" ] && [ -n "$SIM_DIR" ] || die "Bad sim stage or dir: $SIM_STAGE $SIM_DIR"
    [ -n "$NMESHPARTS" ] && [ -n "$REMOTE" ] || die "Bad globals NMESHPARTS: $NMESHPARTS or REMOTE: $REMOTE"

    local SIM_STAGE_UPPER=$(echo "$SIM_STAGE" | tr '[:lower:]' '[:upper:]')
    local SIM_STAGE_LOWER=$(echo "$SIM_STAGE" | tr '[:upper:]' '[:lower:]')

    ensure_files "xns.${SIM_STAGE_LOWER}.in"
    ensure_dirs "$SIM_STAGE_UPPER"/{,mesh}

    mkdir "$SIM_STAGE_UPPER/$SIM_DIR"
    [ -f "xns.${SIM_STAGE_LOWER}.in" ] && cp "xns.${SIM_STAGE_LOWER}.in" $SIM_STAGE_UPPER/$SIM_DIR/xns.in
    [ -f "mixd2pvtu.${SIM_STAGE_LOWER}.in" ] && cp mixd2pvtu.${SIM_STAGE_LOWER}.in $SIM_STAGE_UPPER/$SIM_DIR/mixd2pvtu.in

    [[ "$SIM_STAGE_UPPER" == "MASS" ]] && mapflow_wrapper "$SIM_DIR"

    if [[ "$DISPATCH" == "REMOTE" ]] ; then 
        ensure_run mirror -m -y --delete push $REMOTE -f "$SIM_STAGE_UPPER"
        proclaim "Submitting job on $REMOTE"
        cd "$SIM_STAGE_UPPER/$SIM_DIR"
        local JRUN_OUT=$(ensure_run mirror cmd 'jrun -x -nt '$NTPN' -np '$NMESHPARTS' -n -ne' --target $REMOTE)
        JOB_ID=$(echo "$JRUN_OUT" | tail -n 1 | grep Submitted | awk '{print $2}')
        cd "$BASE"
    elif [[ "$DISPATCH" == "JURECA" ]] ; then
        cd "$SIM_STAGE_UPPER/$SIM_DIR"
        local JRUN_OUT=$(jrun -x -nt $NTPN -np $NMESHPARTS -n -ne)
        JOB_ID=$(echo "$JRUN_OUT" | tail -n 1 | grep Submitted | awk '{print $2}')
        cd "$BASE"
    fi
}

function wait_for_results()
{
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
    ensure_run mixd2pvtu mixd2pvtu.in
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

function prepare_mesh()
{
    [[ -f mesh_column.msh2 ]] || die "No such file: mesh_column.msh2"
    if ! check_files FLOW/mesh/{mxyz,mien,mrng,minf} MASS/mesh/{mxyz,mien,mrng,minf} ; then
        # TODO: Ensure chroma.sh is error-code compliant. 
        # TODO: Ensure that all component tools are error-code compliant
        ensure_run chroma.sh mesh_column.msh2 -n $NMESHPARTS -l 
    fi
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
        ensure_run gendual -e "$ETYPE" && genneim --nen "$NEN" "$SPACETIME_ARG"
    fi

    if ! check_files {mprm,nprm}.${NMP_04} ; then
        ensure_run "$DECOMPOSE_COMMAND" -e "$ETYPE" -p "$NMESHPARTS" -w $SPACETIME_ARG
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
        prepare_mesh
    elif [[ "$MODE" == "RUN" ]]; then
        generate_mesh
        prepare_mesh
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

## MODE can be one of [ MESH, PREPARE, RUN, WAIT ]
MODE="RUN"

# DISPATCH can be one of [ JURECA , REMOTE, LOCAL ]. 
# LOCAL => Run locally on a machine/node (NOT IMPLEMENTED YET)
# REMOTE => Run script locally, sync data to/from remote. Use jrun to submit data. Currently tightly coupled with JURECA.
# JURECA => Run script locally on JURECA. Submit job with jrun.
DISPATCH="JURECA" 
JOB_ID=

## Currently duplicated decomposition from chroma.sh
DECOMPOSE_COMMAND="decompose.metis"
ETYPE=tet && NEN=4 && MESH_ORDER=1

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
            DISPATCH="$2"
            ensure_match "^(JURECA|REMOTE)$" "$DISPATCH"
            shift
            ;;
        -w|--wait)
            MODE="WAIT"
            JOB_ID=$(filter_integer "$2")
            [[ -n "$JOB_ID" ]] || die "Bad JOB_ID provided to --wait"
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

ensure_match "^(MESH|PREPARE|RUN|WAIT|CONVERT)$" "$MODE"

if [[ "$DISPATCH" == "JURECA" ]]; then
    ensure_match "jureca" $(hostname)
fi

if [[ $(hostname) =~ "jureca" ]]; then 
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

driver
