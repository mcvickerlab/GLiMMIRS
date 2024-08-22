#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

# outdir=$HOME/crisprqtl_sim/sim_data
outdir=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_interactions_2guides_power_analysis_28-11-2023/

mkdir -p $outdir

h5=$outdir/sim.h5

if [ -f "$h5" ] ; then
    rm "$h5"
fi

# # simulate data resembling pilot 
# Rscript $PWD/sim_counts_interactions_power_analysis.R --out $outdir \
# --guide_disp 1 10 100 --lambda 15 25 50 75 100 \
# --effect 0.5 1 3 5 7 \
# --png --tiff --pdf

# simulate data resembling at scale experiment
Rscript $PWD/sim_counts_interactions_power_analysis.R --out $outdir \
--lambda 15 20 30 50  \
--cells 50000 \
--genes 1000 \
--effect 0.5 1 2 3 4 5 6 \
--png --tiff --pdf

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "sim_counts_interactions for power analysis 28-Nov-2023 (positive interaction effects) finished successfully" &>/dev/null
	else
		slack "sim_counts_interactions for power analysis 28-Nov-2023 (positive interaction effects) exited with error code $exit_code"
	fi
fi
exit "$exit_code"
