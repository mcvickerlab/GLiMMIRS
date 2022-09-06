#! /bin/bash
#$ -V
#$ -cwd

# outdir=$HOME/crisprqtl_sim/sim_data
outdir=$HOME/crisprQTL/sim_performance_D100_4guides_targeting

mkdir -p $outdir

h5=$HOME/crisprQTL/simulated_data_4guides_x1prob/sim.h5

if [ -f "$outdir/x1_with_D100_efficiencies.h5" ] ; then
    rm "$outdir/x1_with_D100_efficiencies.h5"
fi

Rscript $PWD/test_model_on_simulated_with_estimated_efficiency.R --h5 $h5 \
--targeting --out $outdir --guide_disp 100 

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "test_model_on_simulated_with_estimated_efficiency finished successfully" &>/dev/null
	else
		slack "test_model_on_simulated_with_estimated_efficiency exited with error code $exit_code"
	fi
fi
exit "$exit_code"
