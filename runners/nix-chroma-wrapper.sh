#!/usr/bin/env bash

# An example wrapper script to generate and run multiple simulations

dirs=(0.03 0.04)
stabs=(1.0 0.5 0.1 0.05 0.01)
ROOT="$PWD"

## Prepare mesh in all given subdirs
for dir in "${dirs[@]}"; do 
    cd "$dir"
        nix-chroma.sh -n 640 --prepare &
    cd "$ROOT"
done
wait

## Patch config and run simulation 
for dir in "${dirs[@]}"; do 
    cd "$dir"
    for stab in "${stabs[@]}"; do 
        mkdir -p FLOW/sim_stab_${stab} 
        cp xns.flow.in FLOW/sim_stab_${stab}/xns.in
        sed -i "/^tau_momentum_factor/c\tau_momentum_factor ${stab}" FLOW/sim_stab_${stab}/xns.in
        sed -i "/^tau_continuity_factor/c\tau_continuity_factor ${stab}" FLOW/sim_stab_${stab}/xns.in
        nix-chroma.sh flow -n 640 -nt 128 -s sim_stab_${stab} &
    done
    cd "$ROOT"
done
wait
