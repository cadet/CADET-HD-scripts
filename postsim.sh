#!/bin/bash

set -euo pipefail

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


DIR="${1:-prev}"

## Old method to find Timestep
# TS=$(grep 'Convergence' xns.log | tail -2 | head -1 | awk '{print $6}')
# echo "Extracting timestep: $TS"
# exts -m ../mesh/minf -f data.all -n 2 -o data.in -t "$TS"

MESHDIR=$(findup . -type d -iname "mesh")
[ -n "$MESHDIR" ] && minffile="$MESHDIR/minf" || exit
[ -f "$minffile" ] && [ -n "$minffile" ] || exit
echo "Using $minffile"
NN=$(awk '/^nn/{print $2}' "$minffile")
NREC=$(( $(stat --printf="%s" data.all) / ( NN * 2 * 8 ) ))
echo "Found $NREC records in data.all"
LASTREC=$(( $NREC - 1 ))
echo "Extracting last timestep: $LASTREC"
exts -m "$minffile" -f data.all -n 2 -o data.in -t "$LASTREC"

mkdir "$DIR"
mv chromatogram data.all out-* xns.log "$DIR"
cp xns.in "$DIR"
cp .xns.time "$DIR"

starttime=$(cat .xns.time)

sed -i '/^restart/c\restart on' xns.in
sed -i "/^starttime/c\starttime $starttime" xns.in
