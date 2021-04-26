#!/usr/bin/env bash

## JURECA jobscript wrapper
## usage: Create + Edit + Submit: $sub.sh 
## usage: Edit + Submit: $sub.sh -e
## usage: Pager: $sub.sh -p
## usage: Submit existing (resubmit): $sub.sh -r

die(){
    echo "$@"
    exit
}

title="HDChroma"
timevar=""
nproc=""
tasks_node="128"
partition="dc-cpu"
account="jibg12"
xns="~/bin/xns-"
mapflow="sim"
stitch=""
dep=""
depCookie=""
mail=""
mailCookie=""

if [ "$1" == "-p" ]; then       #pager
    less job.sh
elif [ "$1" == "-e" ]; then     #edit
	$EDITOR job.sh
	sleep 1 
    sbatch ./job.sh
elif [ "$1" != "-r" ]; then     #not resubmit

read -p "title: " -e -i "$title" title
read -p "time (H:M:S): " -e -i "$timevar" timevar
read -p "nproc: " -e -i "$nproc" nproc

[ -f "../mesh/mprm.$(printf '%04d' $nproc)" ] || die "No partition file found: ../mesh/mprm.$(printf '%04d' $nproc)"

read -p "tasks_node: " -e -i "$tasks_node" tasks_node
read -p "partition: " -e -i "$partition" partition
read -p "account: " -e -i "$account" account
read -p "xns: " -e -i "$xns" xns
read -p "mapflow: " -e -i "$mapflow" mapflow
read -p "dep: " -e -i "$dep" dep
read -p "mail: " -e -i "$mail" mail

## dependent simulation job id
[ -n "$dep" ] && depCookie="#SBATCH -d afterany:$dep"

[ -n "$mail" ] && mailCookie="#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.rao@fz-juelich.de"

[ -n "$mapflow" ] && mapflowCMD="mapflow_wrapper $mapflow"

cat > job.sh <<EOFMARKER
#!/usr/bin/env bash

$depCookie
$mailCookie
#SBATCH -J $title
#SBATCH -n $nproc
#SBATCH --ntasks-per-node=$tasks_node
#SBATCH -o out-%j.out 
#SBATCH -e out-%j.err 
#SBATCH --time=$timevar
#SBATCH -A $account
#SBATCH --partition=$partition

# source functions
source ~/cjibg12/dev/tools/scripts/jobtools.sh

# module loads
module load intel-para Boost flex

$mapflowCMD

## stitch flow
# stitch_flow inlet sim

## stitch mass
# stitch_mass inlet sim

# run MPI application below (with srun)
srun -n $nproc --ntasks-per-node=$tasks_node $xns < xns.in
EOFMARKER

    sleep 1
    $EDITOR job.sh
    sbatch ./job.sh

else 
    ## resubmit
    sleep 1
    sbatch ./job.sh
fi
