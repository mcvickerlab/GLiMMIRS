#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

### define simulation parameters
#outdir=/iblm/netapp/data1/jezhou/crisprQTL/sim_base_data/
outdir=data/simulated/base
d=2 # nguides per target
ngenes=13000
ntargets=1000
ncells=50000
lambda=15
guide_disp='1 10 100' # noisy guide efficiencies to simulate
mkdir -p $outdir

h5=$outdir/sim.h5

if [ -f "$h5" ] ; then
    rm "$h5"
fi

# run script 
Rscript src/simulations/data/sim_counts.R --out $outdir \
--genes $ngenes \
--cells $ncells \
--targets $ntargets \
--lambda $lambda \
--guide_disp $guide_disp \
--d $d

# # message the user on slack if possible
# exit_code="$?"
# if command -v 'slack' &>/dev/null; then
#     if [ "$exit_code" -eq 0 ]; then
# 		slack "sim_counts finished successfully" &>/dev/null
# 	else
# 		slack "sim_counts exited with error code $exit_code"
# 	fi
# fi
# exit "$exit_code"
