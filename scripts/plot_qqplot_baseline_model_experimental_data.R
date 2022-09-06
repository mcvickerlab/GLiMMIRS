# This program plots a qqplot of the p-values obtained by running the baseline model on the
# experimental data, and the two negative control sets. 
# This program was written by Karthik Guruvayurappan.

library(stats);
library(BoutrosLab.plotting.general)

# read in Gasperini paper p-values
published.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/gasperini_enhancer_gene_pairs_suppl_table_2.csv')
published.pvalues <- published.pvalues[, c('Target_Site', 'ENSG', 'Diff_expression_test_raw_pval')]
colnames(published.pvalues) <- c('enhancer', 'gene', 'pvalue')
published.pvalues <- published.pvalues[complete.cases(published.pvalues), ]
published.pvalues <- published.pvalues[order(published.pvalues$pvalue), ]
published.pvalues$unif <- seq(1, nrow(published.pvalues), length.out = nrow(published.pvalues)) / nrow(published.pvalues)
published.pvalues$set <- 'Gasperini'

# read in baseline model p-values
baseline.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model.csv')
colnames(baseline.pvalues) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues <- baseline.pvalues[complete.cases(baseline.pvalues), ]
baseline.pvalues <- baseline.pvalues[order(baseline.pvalues$pvalue), ]
baseline.pvalues$unif <- seq(1, nrow(baseline.pvalues), length.out = nrow(baseline.pvalues)) / nrow(baseline.pvalues)
baseline.pvalues$set <- 'Baseline'

baseline.pvalues.scrambled.guide <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model_neg_scrambled_guides.csv')
colnames(baseline.pvalues.scrambled.guide) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues.scrambled.guide <- baseline.pvalues.scrambled.guide[complete.cases(baseline.pvalues.scrambled.guide), ]
baseline.pvalues.scrambled.guide <- baseline.pvalues.scrambled.guide[order(baseline.pvalues.scrambled.guide$pvalue), ]
baseline.pvalues.scrambled.guide$unif <- seq(1, nrow(baseline.pvalues.scrambled.guide), length.out = nrow(baseline.pvalues.scrambled.guide)) / nrow(baseline.pvalues.scrambled.guide)
baseline.pvalues.scrambled.guide$set <- 'Scrambled Guide'

baseline.pvalues.mismatch.gene <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model_neg_mismatch_gene.csv')
colnames(baseline.pvalues.mismatch.gene) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues.mismatch.gene <- baseline.pvalues.mismatch.gene[complete.cases(baseline.pvalues.mismatch.gene), ]
baseline.pvalues.mismatch.gene <- baseline.pvalues.mismatch.gene[order(baseline.pvalues.mismatch.gene$pvalue), ]
baseline.pvalues.mismatch.gene$unif <- seq(1, nrow(baseline.pvalues.mismatch.gene), length.out = nrow(baseline.pvalues.mismatch.gene)) / nrow(baseline.pvalues.mismatch.gene)
baseline.pvalues.mismatch.gene$set <- 'Mismatch Gene'

plot.df <- rbind(published.pvalues, baseline.pvalues, baseline.pvalues.scrambled.guide, baseline.pvalues.mismatch.gene)
plot.df$pvalue[plot.df$pvalue == 0] <- 2.2e-308
plot.df$unif <- -log10(plot.df$unif)
plot.df$pvalue <- -log10(plot.df$pvalue)


create.scatterplot(
    formula = pvalue ~ unif,
    data = plot.df,
    groups = plot.df$set,
    col = default.colours(4),
    xlab.label = NULL,
    ylab.label = NULL,
    alpha  = 0.5,
    key = list(
        text = list(
            lab = c('Baseline','Gasperini','Mismatch Gene', 'Scrambled Guide'),
            cex = 1,
            col = 'black'
        ),
        points = list(
            pch = 19,
            col = default.colours(4),
            cex = 1
        ),
        x = 0.04,
        y = 0.95,
        padding.text = 2
    ),
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/baseline_model_experimental_data_qqplot.tiff',
    resolution = 200
)

plot.df <- rbind(baseline.pvalues.scrambled.guide, baseline.pvalues.mismatch.gene)
plot.df$pvalue[plot.df$pvalue == 0] <- 2.2e-308
plot.df$unif <- -log10(plot.df$unif)
plot.df$pvalue <- -log10(plot.df$pvalue)

create.scatterplot(
    formula = pvalue ~ unif,
    data = plot.df,
    groups = plot.df$set,
    col = default.colours(2),
    xlab.label = NULL,
    ylab.label = NULL,
    alpha  = 0.5,
    key = list(
        text = list(
            lab = c('Mismatch Gene', 'Scrambled Guide'),
            cex = 1,
            col = 'black'
        ),
        points = list(
            pch = 19,
            col = default.colours(2),
            cex = 1
        ),
        x = 0.04,
        y = 0.95,
        padding.text = 2
    ),
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/baseline_model_experimental_data_qqplot_neg_zoom_in.tiff',
    resolution = 200
)
