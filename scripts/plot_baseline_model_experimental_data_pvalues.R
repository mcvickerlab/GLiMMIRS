# This program plots the distribution of p-values obtained from running the baseline model on the
# 664 enhancer-gene pairs published in Gasperini et al. to check for any pileups of p-values,
# indicating optimization issues with the model.
#
# Author: Karthik Guruvayurappan

library(BoutrosLab.plotting.general)

baseline.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model.csv')
colnames(baseline.pvalues) <- c('enhancer', 'gene', 'pvalue')
create.histogram(
    x = baseline.pvalues$pvalue,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/baseline_model_p_value_distribution.tiff',
    resolution = 200
)
