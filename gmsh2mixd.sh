#!/usr/bin/env bash
## Should work for msh2 files with 3D meshes
## I can probably rewrite the whole thing in awk. 
## Started as a quick and fun test to see if I can do this without too much effort. 
## Also a good lesson on the power of unix tools (mostly).

ELEM="TET"
NEN=4
ETYPE=4

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -tri)
            ELEM="TRI"
            NEN=3
            ETYPE=2
            shift # past value
            ;;
        -tet)
            ELEM="TET"
            NEN=4
            ETYPE=4
            shift # past value
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


## Get nodes
# sed '1,/\$Nodes/{N;d};/\$EndNodes/,$d' "$1" | awk '{for (i=2; i<=NF; i++) print $i}' > mxyz.asc
sed '1,/\$Nodes/d;/\$EndNodes/,$d' "$1" | tail -n +2 | awk '{for (i=2; i<=NF; i++) print $i}' > mxyz.asc


## elm-number elm-type number-of-tags < tag > ... node-number-list

sed '1,/\$Elements/d;/\$EndElements/,$d' "$1" | tail -n +2 | awk -v n="$NEN" -v e="$ETYPE" '$2 == e {for (i=NF-n+1; i<=NF; i++) print $i}' > mien.asc 
# sed '1,/\$Elements/d;/\$EndElements/,$d' "$1" | tail -n +2 | awk '$2 == 2 {for (i=NF-2; i<=NF; i++) print $i}' > mien.asc

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
## TODO: pull elements of a specific tag
