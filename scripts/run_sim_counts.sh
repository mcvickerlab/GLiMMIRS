#! /bin/bash
#$ -V
#$ -cwd

<<<<<<< HEAD
# outdir=$HOME/crisprqtl_sim/sim_data
outdir=$HOME/crisprQTL/simulated_data
=======
outdir=$HOME/crisprqtl_sim/sim_data
>>>>>>> 484e4f3229fb141ed4efa5485bbaa02b24273cbf

mkdir -p $outdir

Rscript $PWD/sim_counts.R --out $outdir

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