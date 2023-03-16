#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

# outdir=$HOME/crisprqtl_sim/sim_data
outdir=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_4guides_discrete_and_continuous_10-10-2022_4guides/


mkdir -p $outdir

h5=$outdir/sim.h5

if [ -f "$h5" ] ; then
    rm "$h5"
fi

Rscript $PWD/sim_counts.R --out $outdir --guide_disp 1 10 100 --d 4


# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "sim_counts finished successfully" &>/dev/null
	else
		slack "sim_counts exited with error code $exit_code"
	fi
fi
exit "$exit_code"
