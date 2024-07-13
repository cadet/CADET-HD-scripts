#!/usr/bin/env bash

TOOLS=(genmesh pymesh mixdtools paravision chromoo scripts)
URL=git@github.com:modsim

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        --https)
            URL=https://github.com/modsim; shift ;;
        *)
            POSITIONAL+=("$1"); shift ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

for TOOL in "${TOOLS[@]}"; do
    git clone "$URL/ChromaHD-$TOOL" "$TOOL"
done
