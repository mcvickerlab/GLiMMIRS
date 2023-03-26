#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err


### model specifications:
### using 4 gRNAs per target
### evaluating only targeted genes
### using continuous X1 values CALCULATED WITH noisy guide efficiencies (D=100)
### using counts simulated with continuous X1 values


# outdir=$HOME/crisprqtl_sim/sim_data
outdir=/iblm/netapp/data1/jezhou/crisprQTL/sim_performance_D100-noisy_4guides_cont-x1_cont-counts_targeting
mkdir -p $outdir

h5=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_4guides_discrete_and_continuous_10-10-2022_4guides/sim.h5

# if true, delete existing directory with same ename 
override=true

if [ override ]; then
	rm -rf $outdir && mkdir -p $outdir
fi

Rscript $PWD/test_model_on_simulated_with_noisy_efficiency.R --h5 $h5 \
--targeting --out $outdir --d 4 --x1 continuous --counts continuous --guide_disp 100 \
--pseudocount 0.1

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "test_model_on_simulated_with_true_efficiency using 4 guides, continuous x1 with noisy guide efficiency (D=100), continuous counts finished successfully" &>/dev/null
	else
		slack "test_model_on_simulated_with_true_efficiency using 4 guides, continuous x1 with noisy guide efficiency (D=100), continuous counts exited with error code $exit_code"
	fi
fi
exit "$exit_code"