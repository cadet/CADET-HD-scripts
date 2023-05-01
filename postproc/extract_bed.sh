#!usr/bin/env bash

# A template of commands to run for extracting packed bed mesh with data from the full MASS simulations

exit

set -euo pipefail

cd $MASS_MESH_DIR
mkdir bed-mesh

# Removes interstitial elements from the full mass mesh. Stores results into ./bed-mesh
rmmat -tet -st bed-mesh 1

cd $MASS_SIM_DIR

# split the full mass data file into bed_c and bed_q respectively
mixdsplit -m $MASS_MESH_DIR -o bed_c.all -i 0
mixdsplit -m $MASS_MESH_DIR -o bed_q.all -i 1

# Then use the following mixd2pvtu config template to generate the pvtu files
# srun -N 1 -n 47 --time=02:00:00 mixd2pvtu mixd2pvtu.in
######### MIXD2PVTU.in config example
# outpath bed_c
# title bed_c
#
# mien /p/project/cjibg12/simulations/10k-beads/AC-poly-0.06-re/MASS/mesh/bed-mesh/mien
# minf /p/project/cjibg12/simulations/10k-beads/AC-poly-0.06-re/MASS/mesh/bed-mesh/minf
# mxyz /p/project/cjibg12/simulations/10k-beads/AC-poly-0.06-re/MASS/mesh/bed-mesh/mxyz
#
# elemtype tet
# nrec 99999
# data /p/project/cjibg12/simulations/10k-beads/AC-poly-0.06-re/MASS/sim/01/bed_c.all /p/project/cjibg12/simulations/10k-beads/AC-poly-0.06-re/MASS/sim/02/bed_c.all 
# ndf 1
# dtFile /p/home/jusers/rao2/jureca/cjibg12/simulations/10k-beads/AC-poly-0.06-re/MASS/sim/timesteps_01_02.txt
