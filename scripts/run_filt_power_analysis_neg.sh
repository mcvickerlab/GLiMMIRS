#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err
#$ -t 1-5

lambda=`awk -v line=$SGE_TASK_ID 'NR==line' $HOME/crisprQTL/scripts/lambdas_list.txt`

echo "power analysis for lambda=$lambda, negative cases filtering for min. 20 cells"

outdir=/iblm/netapp/data1/jezhou/crisprQTL/power_analysis_01-09-2023_filt20
simdata=/iblm/netapp/data1/jezhou/crisprQTL/simulated_data_interactions_2guides_power_analysis_01-09-2023/sim.h5

mkdir -p $outdir

Rscript $PWD/filt_power_analysis.R --h5 $simdata --out $outdir --lambda $lambda --neg


# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "filtered power analysis for lambda=$lambda (negative cases) min. 20 cells finished successfully" &>/dev/null
	else
		slack "filtered power analysis for lambda=$lambda (negative cases) min. 20 cells exited with error code $exit_code"
	fi
fi
exit "$exit_code"
