#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err
#$ -t 1-5

lambda=`awk -v line=$SGE_TASK_ID 'NR==line' $HOME/crisprQTL/scripts/lambdas_list.txt`

echo "power analysis for lambda=$lambda, positive cases (with interaction)"

outdir=/iblm/netapp/data1/jezhou/crisprQTL/power_analysis_01-09-2023
simdata=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_interactions_2guides_power_analysis_01-09-2023/sim.h5

mkdir -p $outdir

Rscript $PWD/power_analysis.R --h5 $simdata --out $outdir --lambda $lambda


# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "power analysis for lambda=$lambda (positive cases) 09-Jan-2023 finished successfully" &>/dev/null
	else
		slack "power analysis for lambda=$lambda (positive cases) 09-Jan-2023 exited with error code $exit_code"
	fi
fi
exit "$exit_code"
