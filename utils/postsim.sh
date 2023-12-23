#!/usr/bin/env bash
## Postsim script
# Run this after a simulation run in order to make restarting easy
# Copies results into a subdirectory, modifies xns.in and .xns.time files to appropriate time
#
# usage: ./postsim.sh -n [ndf:int]

function die(){
    echo -e "ERROR: $@" >&2
    exit -1
}

findup()
{
    path="$1"
    shift 1
    while [[ $path != / ]];
    do
        find "$path" -maxdepth 1 -mindepth 1 "$@"
        # Note: if you want to ignore symlinks, use "$(realpath -s "$path"/..)"
        path="$(readlink -f "$path"/..)"
    done
}

filter_integer() {
    if [[ $1 =~ ^[[:digit:]]+$ ]]; then
        echo "$1"
    fi
}

NDF=2
XNS_CONF_FILE=xns.in

## Commandline args processing
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -n|--ndf)
            NDF=$2
            shift # past value
            shift # past value
            ;;
        -x|--xns-conf)
            XNS_CONF_FILE=$2
            shift # past value
            shift # past value
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


# Process and find last executed timestep
# Extract solution for this timestep into data.in
MESHDIR=$(findup . -type d -iname "mesh") && [[ -d "$MESHDIR" ]] || die "No mesh dir found!"
minffile="$MESHDIR/minf" && [[ -f "$minffile" ]] || die "No minf file found!"
echo "Using $minffile"
NN=$(awk '/^nn/{print $2}' "$minffile")

[[ -f data.all ]] || die "No data.all found!"
NREC=$(( $(stat --printf="%s" data.all) / ( NN * NDF * 8 ) ))
echo "Found $NREC records in data.all"

LASTREC=$(( $NREC - 1 ))
echo "Extracting last timestep: $LASTREC"
extractTS -m "$minffile" -f data.all -n $NDF -o data.in -t "$LASTREC"

STARTTIME=$(cat .xns.time)
DIR_PREFIX=$(date '+%Y_%m_%d_%H_%M')
STARTTIME_TRIMMED=$(echo "$starttime" | tr -d '[:blank:]')
DIR="${DIR_PREFIX}_${STARTTIME_TRIMMED%.*}_${NREC}"

mkdir "$DIR"
echo "created dir: $DIR"

## Move solution and log files
MV_FILES=(chromatogram data.all xns.log)
for mvfile in "${MV_FILES[@]}" ; do
    if [ -f "$mvfile" ]; then 
        mv "$mvfile" "$DIR"
    fi
done

## copy config files
CP_FILES=(xns.in xns.flow.in xns.mass.in .xns.time)
for cpfile in "${CP_FILES[@]}" ; do
    if [ -f "$cpfile" ]; then
        cp "$cpfile" "$DIR"
    fi
done

sed -i '/^restart/c\restart on' $XNS_CONF_FILE
sed -i "/^starttime/c\starttime $STARTTIME" $XNS_CONF_FILE
