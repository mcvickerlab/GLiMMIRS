#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err


### model specifications:
### using 4 gRNAs per target
### evaluating only targeted genes
### using continuous X1 values
### using noisy guide efficiencies (D=10)
### using counts simulated with continuous X1 values


# outdir=$HOME/crisprqtl_sim/sim_data
outdir=/iblm/netapp/data1/jezhou/crisprQTL/sim_performance_D10-noisy_4guides_cont-x1_cont-counts_targeting
mkdir -p $outdir

h5=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_4guides_discrete_and_continuous/sim.h5

# if [ -f "$outdir/x1_with_D100_efficiencies.h5" ] ; then
#     rm "$outdir/x1_with_D100_efficiencies.h5"
# fi

Rscript $PWD/test_model_on_simulated_with_noisy_efficiency.R --h5 $h5 \
--targeting --out $outdir --d 4 --x1 continuous --counts continuous --guide_disp 10

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "test_model_on_simulated_with_true_efficiency using 4 guides, continuous x1 with noisy guide efficiency (D=10), continuous counts finished successfully" &>/dev/null
	else
		slack "test_model_on_simulated_with_true_efficiency using 4 guides, continuous x1 with noisy guide efficiency (D=10), continuous counts exited with error code $exit_code"
	fi
fi
exit "$exit_code"
