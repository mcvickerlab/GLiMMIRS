# This program plots a qqplot of the p-values obtained by running the baseline model on the
# experimental data, and the two negative control sets. 
# This program was written by Karthik Guruvayurappan.

library(BoutrosLab.plotting.general)

# read in Gasperini paper p-values
published.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/gasperini_enhancer_gene_pairs_suppl_table_2.csv')
published.pvalues <- published.pvalues[, c('Target_Site', 'ENSG', 'Diff_expression_test_raw_pval')]
colnames(published.pvalues) <- c('enhancer', 'gene', 'pvalue')
print(head(published.pvalues, 20))
published.pvalues$set <- 'Gasperini'

# read in baseline model p-values
baseline.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model.csv')
colnames(baseline.pvalues) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues$set <- 'Baseline'

baseline.pvalues.scrambled.guide <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model_neg_scrambled_guides.csv')
colnames(baseline.pvalues.scrambled.guide) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues.scrambled.guide$set <- 'Scrambled Guide'

baseline.pvalues.mismatch.gene <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model_neg_mismatch_gene.csv')
colnames(baseline.pvalues.mismatch.gene) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues.mismatch.gene$set <- 'Mismatch Gene'

plot.df <- rbind(published.pvalues, baseline.pvalues, baseline.pvalues.scrambled.guide, baseline.pvalues.mismatch.gene)

create.qqplot.fit(
    x = ~ pvalue | set,
    data = plot.df,
    distribution = 'qunif',
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/baseline_model_experimental_data_qqplot.tiff',
    resolution = 200
)
