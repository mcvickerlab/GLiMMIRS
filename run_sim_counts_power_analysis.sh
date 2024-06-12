#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

#outdir=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_interactions_2guides_power_analysis_01-09-2023/
outdir=data/simulated/int_power_analysis
d=2 # nguides per target
ngenes=13000
ntargets=1000
ncells=50000
lambdas='15 25 50 75 100'
effects='0.5 1 3 5 7'

mkdir -p $outdir

h5=$outdir/sim.h5

if [ -f "$h5" ] ; then
    rm "$h5"
fi

Rscript src/data/sim_counts_interactions_power_analysis.R --out $outdir \
--genes $ngenes \
--targets $ntargets \
--cells $ncells \
--d $d \
--lambda $lambdas \
--effect $effects

## message the user on slack if possible
#exit_code="$?"
#if command -v 'slack' &>/dev/null; then
#    if [ "$exit_code" -eq 0 ]; then
#		slack "sim_counts_interactions for power analysis 09-Jan-2023 finished successfully" &>/dev/null
#	else
#		slack "sim_counts_interactions for power analysis 09-Jan-2023 exited with error code $exit_code"
#	fi
#fi
#exit "$exit_code"
