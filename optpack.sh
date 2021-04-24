#!/usr/bin/bash

## Optimize a packing generated with PackingGeneration.exe to match a particular porosity

exec |& tee out.opti.log

while :
do
    sed -i "/Seed/s/^.*/Seed: $RANDOM/" generation.conf
    rm packing.nfo
    cp -f packing.fba.xyzd packing.xyzd
    PackingGeneration.exe -ls |& tee out.ls.log
    val=$(grep 'Calc. porosity' out.ls.log | tail -n 1 | awk '{print $4}' | cut -b 1-5)
    if [[ "$val" == "0.407" ]]; then
        break;
    fi
done
