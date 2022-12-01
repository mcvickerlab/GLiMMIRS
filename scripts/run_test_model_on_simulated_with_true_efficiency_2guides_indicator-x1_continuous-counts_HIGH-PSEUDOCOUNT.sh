#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err


### model specifications:
### using 2 gRNAs per target
### evaluating only targeted genes
### using indicator X1 values
### using counts simulated with continuous X1 values
### using pseudocount of 0.5

# outdir=$HOME/crisprqtl_sim/sim_data
outdir=/iblm/netapp/data1/jezhou/crisprQTL/sim_performance_true_efficiency_2guides_indic-x1_cont-counts_targeting_pseudocount-HIGH
mkdir -p $outdir

h5=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_2guides_discrete_and_continuous_10-10-2022/sim.h5

# if true, delete existing directory with same ename 
override=true

if [ override ]; then
	rm -rf $outdir && mkdir -p $outdir
fi

Rscript $PWD/test_model_on_simulated_with_true_efficiency.R --h5 $h5 \
--targeting --out $outdir --d 2 --x1 indicator --counts continuous \
--pseudocount 0.01

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "test_model_on_simulated_with_true_efficiency using 2 guides, indicator x1, continuous counts, pseudocount=0.5 finished successfully" &>/dev/null
	else
		slack "test_model_on_simulated_with_true_efficiency using 2 guides, indicator x1, continuous counts, pseudocount=0.5 exited with error code $exit_code"
	fi
fi
exit "$exit_code"
