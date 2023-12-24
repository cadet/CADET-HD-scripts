#!/usr/bin/env -S nix develop --impure /home/jayghoshter/templates/nix-flakes/pymesh/flake.nix --command zsh

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
    "$@"
    [[ $? == 0 ]] || die "command exited with error: $@"
}

NMESHPARTS=48
BASE="$PWD"
REMOTE="jureca"
SIM_STAGE="FLOW"
SIM_DIR="sim"
SIM_COMPLETION_CHECK_TIME=60

SIM_STAGE_UPPER=$(echo "$SIM_STAGE" | tr '[:lower:]' '[:upper:]')
SIM_STAGE_LOWER=$(echo "$SIM_STAGE" | tr '[:upper:]' '[:lower:]')

if [[ ! -f mesh_column.msh2 ]]; then
    proclaim "Generating mesh"
    ensure_run mesh input.yaml
fi

[[ -f mesh_column.msh2 ]] || die "No such file: mesh_column.msh2"

chroma.sh mesh_column.msh2 -n $NMESHPARTS -l 
mkdir "$SIM_STAGE_UPPER/$SIM_DIR"

[ -f "xns.${SIM_STAGE_LOWER}.in" ] && cp "xns.${SIM_STAGE_LOWER}.in" $SIM_STAGE_UPPER/$SIM_DIR/xns.in
[ -f "mixd2pvtu.${SIM_STAGE_LOWER}.in" ] && cp mixd2pvtu.${SIM_STAGE_LOWER}.in $SIM_STAGE_UPPER/$SIM_DIR
mirror -m -y --delete push $REMOTE -f "$SIM_STAGE_UPPER"

proclaim "Submitting job on $REMOTE"
cd "$SIM_STAGE_UPPER/$SIM_DIR"
mirror cmd 'jrun -x -np '$NMESHPARTS' -n -ne' --target $REMOTE
cd "$BASE"

proclaim "Waiting for simulation on $REMOTE"
mirror pull $REMOTE -y
while [[ ! -f "$SIM_STAGE_UPPER/$SIM_DIR/data.out" ]]
do
    echo "No simulation results found, waiting $SIM_COMPLETION_CHECK_TIME before trying again"
    sleep $SIM_COMPLETION_CHECK_TIME
    mirror pull $REMOTE -f "$SIM_STAGE_UPPER" -y
done

proclaim "Converting results"
cd "$SIM_STAGE_UPPER/$SIM_DIR"
ensure_run mixd2pvtu mixd2pvtu.${SIM_STAGE_LOWER}.in
cd "$BASE"

# cd "$SIM_STAGE_UPPER/$SIM_DIR/output"
# cd "$BASE"
