#! /bin/bash
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

#outdir=/iblm/netapp/data1/jezhou/crisprQTL/multiplicative_vs_additive_330_pairs_11-May-2023
outdir=out/experimental/multiplicative_vs_additive

mkdir -p $outdir

Rscript src/experimental/models/compare_multiplicative_vs_additive.R \
--out $outdir

## message the user on slack if possible
#exit_code="$?"
#if command -v 'slack' &>/dev/null; then
#    if [ "$exit_code" -eq 0 ]; then
#		slack "fitting multiplicative vs. additive models on high confidence (330) pairs finished successfully" &>/dev/null
#	else
#		slack "fitting multiplicative vs. additive models on high confidence (330) exited with error code $exit_code"
#	fi
#fi
#exit "$exit_code"
