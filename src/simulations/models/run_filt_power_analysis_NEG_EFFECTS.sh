#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err
#$ -t 1-4

wd=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_interactions_2guides_power_analysis_NEG_EFFECT_28-11-2023
# lambda=`awk -v line=$SGE_TASK_ID 'NR==line' $HOME/crisprQTL_dump/scripts/lambdas_list.txt`
lambda=`awk -v line=$SGE_TASK_ID 'NR==line' $wd/lambdas_list.txt`

echo "power analysis for lambda=$lambda, positive cases (with interaction, NEG EFFECT), filtering for min 10 cells"

outdir=/iblm/netapp/data1/jezhou/crisprQTL/power_analysis_NEG_EFFECTS_28-11-2023_filt10
simdata=$wd/sim.h5

mkdir -p $outdir

Rscript $PWD/filt_power_analysis.R --h5 $simdata --out $outdir --lambda $lambda --min.cells 10 --neg_effect --neg


# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "filtered power analysis (28-11-2023 simdata) for lambda=$lambda (negative cases) with neg. interactions & min. 10 cells finished successfully" &>/dev/null
	else
		slack "filtered power analysis (28-11-2023 simdata) for lambda=$lambda (negative cases) with neg. interactions & min. 10 cells exited with error code $exit_code"
	fi
fi
exit "$exit_code"
