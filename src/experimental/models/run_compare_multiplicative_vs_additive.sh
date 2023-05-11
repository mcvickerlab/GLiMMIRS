#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

outdir=/iblm/netapp/data1/jezhou/crisprQTL/multiplicative_vs_additive_330_pairs

mkdir -p $outdir

Rscript $PWD/run_model_experimental_suppl_data_table_2_enhancer_pairs_ADDITIVE.R \
--out $outdir 
# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "fitting multiplicative vs. additive models on 330 enhancer pairs finished successfully" &>/dev/null
	else
		slack "fitting multiplicative vs. additive models on 330 enhancer pairs exited with error code $exit_code"
	fi
fi
exit "$exit_code"
