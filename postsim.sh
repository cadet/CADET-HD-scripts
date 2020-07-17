#!/bin/bash

DIR="${1:-prev}"

# TS=$(grep 'Convergence' $(ls -t *err | head -1) | tail -2 | head -1 | awk '{print $6}')
TS=$(grep 'Convergence' xns.log | tail -2 | head -1 | awk '{print $6}')
echo "Extracting timestep: $TS"
exts -m ../mesh/minf -f data.all -n 2 -o data.in -t "$TS"
mkdir "$DIR"
mv chromatogram data.all out-* xns.log "$DIR"
cp xns.in "$DIR"
cp .xns.time "$DIR"

starttime=$(cat .xns.time)

sed -i '/^restart/c\restart on' xns.in
sed -i "/^starttime/c\starttime $starttime" xns.in
# sed -i "s/ninner 1000/ninner 500/" xns.in
