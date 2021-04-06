#!/usr/bin/env bash
## Should work for msh2 files with 3D meshes
## I can probably rewrite the whole thing in awk. 
## Started as a quick and fun test to see if I can do this without too much effort. 
## Also a good lesson on the power of unix tools (mostly).

## Get nodes
sed '1,/\$Nodes/{N;d};/\$EndNodes/,$d' "$1" | awk '{for (i=2; i<=NF; i++) print $i}' > mxyz.asc

## elm-number elm-type number-of-tags < tag > ... node-number-list

## TODO: Awk with bash variables: -tet => 4 -tri =>3 and respective nens to print
## TODO: pull items of a specific tag

## if elm-type is 4 (lin tet), print the last 4 columns one by one
sed '1,/\$Elements/{N;d};/\$EndElements/,$d' "$1" | awk '$2 == 4 {for (i=NF-3; i<=NF; i++) print $i}' > mien.asc

# ## if elm-type is 2 (lin tri), print the last 3 columns one by one
# sed '1,/\$Elements/{N;d};/\$EndElements/,$d' "$1" | awk '$2 == 2 {for (i=NF-2; i<=NF; i++) print $i}' > mien.asc

## Convert input ascii file into binary files. 
## dumpy's default endianness is -e '>' to handle mixd files.
dumpy mxyz.asc -a -v d -w mxyz 
dumpy mien.asc -a -v i -w mien 


## minf file
NN=$(grep '\$Nodes' "$1" -A1 | tail -n1)
NE=$(grep '\$Elements' "$1" -A1 | tail -n1)
echo -e "nn $NN\nne $NE" >> minf 

## TODO: mrng
## TODO: mtbl & doubling nodes
