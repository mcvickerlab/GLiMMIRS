#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

outdir=/iblm/netapp/data1/jezhou/crisprQTL/GLiMMIRS-base_sim_performance
mkdir -p $outdir

h5=/iblm/netapp/data1/jezhou/crisprQTL/sim_base_data/sim.h5

# if true, delete existing directory with same name 
override=true
if [ override ]; then
	rm -rf $outdir && mkdir -p $outdir
fi

### model specifications
ngenes=13000
ntargets=1000
ncells=50000
perturb=probability # fit model using perturbation probabilities 
pseudocount=0.01 # add pseudocount of 0.01

# fit models - set --targeting flag to only fit models to genes with ground truth perturbation
Rscript $PWD/fit_GLiMMIRS-base_sim.R --h5 $h5 \
--cells $ncells --genes $ngenes \
--targeting --out $outdir --perturb $perturb \
--pseudocount $pseudocount

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "Fitting GLiMMIRS-base on simulated data for production finished successfully" &>/dev/null
	else
		slack "Fitting GLiMMIRS-base on simulated data exited with error code $exit_code"
	fi
fi
exit "$exit_code"
