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
        [[ "$ARG" =~ "$CMP_STR" ]] || die "$ARG doesn't match $CMP_STR"
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

function run_simulation_stage_on_remote()
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

    ensure_run mirror -m -y --delete push $REMOTE -f "$SIM_STAGE_UPPER"

    proclaim "Submitting job on $REMOTE"
    cd "$SIM_STAGE_UPPER/$SIM_DIR"
    ensure_run mirror cmd 'jrun -x -np '$NMESHPARTS' -n -ne' --target $REMOTE
    cd "$BASE"

}

function wait_for_results()
{
    local SIM_STAGE="$1"
    local SIM_DIR="$2"

    [ -n "$SIM_STAGE" ] && [ -n "$SIM_DIR" ] || die "Bad sim stage or dir: $SIM_STAGE $SIM_DIR"
    local SIM_STAGE_UPPER=$(echo "$SIM_STAGE" | tr '[:lower:]' '[:upper:]')
    local SIM_STAGE_LOWER=$(echo "$SIM_STAGE" | tr '[:upper:]' '[:lower:]')

    proclaim "Waiting for simulation on $REMOTE"

    cd "$SIM_STAGE_UPPER/$SIM_DIR"
    
    mirror pull $REMOTE -y
    while mirror cmd 'jrun --running' --target "$REMOTE"
    do
        echo "No simulation results found, waiting $SIM_COMPLETION_CHECK_TIME before trying again"
        sleep $SIM_COMPLETION_CHECK_TIME
        mirror pull $REMOTE -y
    done

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
        ensure run gendual -e "$ETYPE" && genneim --nen "$NEN" "$SPACETIME_ARG"
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
            run_simulation_stage_on_remote "$SIM_STAGE" "$SIM_NAME"
            wait_for_results "$SIM_STAGE" "$SIM_NAME"
            convert_results "$SIM_STAGE" "$SIM_NAME"
        done
    elif [[ "$MODE" == "WAIT" ]]; then
        for SIM_STAGE in ${SIM_STAGES[@]}; do 
            wait_for_results "$SIM_STAGE" "$SIM_NAME"
        done
    fi

    # cd "$SIM_STAGE_UPPER/$SIM_DIR/output"
    # cd "$BASE"
}

## Globals
NMESHPARTS=48
BASE="$PWD"
REMOTE="jureca"
SIM_COMPLETION_CHECK_TIME=60
SIM_NAME="sim"
SIM_STAGES=(FLOW MASS)
MODE="RUN"
DECOMPOSE_COMMAND="decompose.metis"
ETYPE=tet && NEN=4 && MESH_ORDER=1

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
        -w|--wait)
            MODE="WAIT"
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

driver
