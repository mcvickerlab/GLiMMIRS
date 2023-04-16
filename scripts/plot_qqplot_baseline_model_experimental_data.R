# This program plots a qqplot of the p-values obtained by running the baseline model on the
# experimental data, and the two negative control sets. 
# This program was written by Karthik Guruvayurappan.

library(stats)
library(ggplot2)
library(RColorBrewer)

# read in Gasperini paper p-values
published.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/gasperini_enhancer_gene_pairs_suppl_table_2.csv')
published.pvalues <- published.pvalues[, c('Target_Site', 'ENSG', 'Diff_expression_test_raw_pval')]
colnames(published.pvalues) <- c('enhancer', 'gene', 'pvalue')
published.pvalues <- published.pvalues[complete.cases(published.pvalues), ]
published.pvalues <- published.pvalues[order(published.pvalues$pvalue), ]
published.pvalues$unif <- seq(1, nrow(published.pvalues), length.out = nrow(published.pvalues)) / nrow(published.pvalues)
published.pvalues$set <- 'Gasperini'

# read in baseline model p-values
baseline.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_gene_pairs_suppl_table_2_baseline_pseudocount_model.csv')
colnames(baseline.pvalues) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues <- baseline.pvalues[complete.cases(baseline.pvalues), ]
baseline.pvalues <- baseline.pvalues[order(baseline.pvalues$pvalue), ]
baseline.pvalues$unif <- seq(1, nrow(baseline.pvalues), length.out = nrow(baseline.pvalues)) / nrow(baseline.pvalues)
baseline.pvalues$set <- 'Baseline'

baseline.pvalues.scrambled.guide <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_gene_pairs_suppl_table_2_baseline_pseudocount_model_neg_scrambled_guides.csv')
colnames(baseline.pvalues.scrambled.guide) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues.scrambled.guide <- baseline.pvalues.scrambled.guide[complete.cases(baseline.pvalues.scrambled.guide), ]
baseline.pvalues.scrambled.guide <- baseline.pvalues.scrambled.guide[order(baseline.pvalues.scrambled.guide$pvalue), ]
baseline.pvalues.scrambled.guide$unif <- seq(1, nrow(baseline.pvalues.scrambled.guide), length.out = nrow(baseline.pvalues.scrambled.guide)) / nrow(baseline.pvalues.scrambled.guide)
baseline.pvalues.scrambled.guide$set <- 'Shuffled Guides'

baseline.pvalues.mismatch.gene <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_gene_pairs_suppl_table_2_baseline_pseudocount_model_neg_mismatch_gene.csv')
colnames(baseline.pvalues.mismatch.gene) <- c('enhancer', 'gene', 'pvalue')
baseline.pvalues.mismatch.gene <- baseline.pvalues.mismatch.gene[complete.cases(baseline.pvalues.mismatch.gene), ]
baseline.pvalues.mismatch.gene <- baseline.pvalues.mismatch.gene[order(baseline.pvalues.mismatch.gene$pvalue), ]
baseline.pvalues.mismatch.gene$unif <- seq(1, nrow(baseline.pvalues.mismatch.gene), length.out = nrow(baseline.pvalues.mismatch.gene)) / nrow(baseline.pvalues.mismatch.gene)
baseline.pvalues.mismatch.gene$set <- 'Mismatch Gene'

plot.df <- rbind(published.pvalues, baseline.pvalues, baseline.pvalues.scrambled.guide, baseline.pvalues.mismatch.gene)
plot.df$pvalue[plot.df$pvalue == 0] <- 2.2e-308
plot.df$unif <- -log10(plot.df$unif)
plot.df$pvalue <- -log10(plot.df$pvalue)

qq.plot <- ggplot(plot.df, aes(x = unif, y = pvalue, color = set)) + 
    geom_point(size = 5) +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Expected -log[10](italic(p)))) + 
    ylab(bquote(Observed -log[10](italic(p)))) +
    theme_classic() +
    theme(
        axis.line = element_line(linewidth = 1),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.text = element_text(size = 32, color = 'black'),
        axis.ticks = element_line(color = 'black', linewidth = 1),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = c(0.10, 0.89),
        legend.text = element_text(size = 24, color = 'black'),
        plot.margin = rep(unit(10, 'mm'), 4),
    ) +
    scale_colour_brewer(palette = 'Set1')

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_16_baseline_model_experimental_data_qqplot.pdf',
    device = 'pdf',
    plot = qq.plot,
    width = 46,
    height = 27,
    units = 'cm'
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_16_baseline_model_experimental_data_qqplot.png',
    device = 'png',
    plot = qq.plot,
    width = 46,
    height = 27,
    units = 'cm'
)

zoom.in.published <- published.pvalues[published.pvalues$pvalue > min(baseline.pvalues.mismatch.gene$pvalue), ]
zoom.in.baseline <- baseline.pvalues[baseline.pvalues$pvalue > min(baseline.pvalues.mismatch.gene$pvalue), ]
plot.df <- rbind(baseline.pvalues.scrambled.guide, baseline.pvalues.mismatch.gene, zoom.in.published, zoom.in.baseline)
plot.df$pvalue[plot.df$pvalue == 0] <- 2.2e-308
plot.df$unif <- -log10(plot.df$unif)
plot.df$pvalue <- -log10(plot.df$pvalue)

qq.plot <- ggplot(plot.df, aes(x = unif, y = pvalue, color = set)) + 
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Expected -log[10](italic(p)))) + 
    ylab(bquote(Observed -log[10](italic(p)))) +
    theme_classic() +
    theme(
        axis.line = element_line(linewidth = 1),
        axis.title.x = element_text(size = 20, color = 'black'),
        axis.title.y = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 20, color = 'black'),
        axis.ticks = element_line(color = 'black', linewidth = 1),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = c(0.18, 0.89),
        legend.text = element_text(size = 8, color = 'black'),
        plot.margin = rep(unit(10, 'mm'), 4),
    ) +
    scale_colour_brewer(palette = 'Set1')

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_28_baseline_model_experimental_data_neg_controls_qqplot.pdf',
    device = 'pdf',
    plot = qq.plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_28_baseline_model_experimental_data_neg_controls_qqplot.png',
    device = 'png',
    plot = qq.plot
)
