#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

### define power analysis parameters
#outdir=/iblm/netapp/data1/jezhou/crisprQTL/test_power_analysis_03-May-2023_NEG
outdir=out/simulated/GLiMMIRS-int_power_analysis_NEG
#simdata=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_interactions_2guides_power_analysis_01-09-2023/sim.h5
simdata=data/simulated/int_power_analysis/sim.h5
ncells=50000
ngenes=13000

mkdir -p $outdir

Rscript src/simulated/models/GLiMMIRS-int_power_analysis.R --h5 $simdata \
--cells $ncells \
--genes $ngenes \
--out $outdir --neg


## message the user on slack if possible
#exit_code="$?"
#if command -v 'slack' &>/dev/null; then
#    if [ "$exit_code" -eq 0 ]; then
#		slack "power analysis for negative pairs finished successfully" &>/dev/null
#	else
#		slack "power analysis for negative pairs exited with error code $exit_code"
#	fi
#fi
#exit "$exit_code"
