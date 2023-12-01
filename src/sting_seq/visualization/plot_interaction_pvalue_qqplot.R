library(ggplot2)
library(RColorBrewer)

# read in interaction p-values
interaction.pvalues <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_interaction_models.csv')
interaction.pvalues$color <- 'black'
interaction.pvalues$color[interaction.pvalues$interaction.fdr.pvalues < 0.1] <- 'red'
interaction.pvalues <- interaction.pvalues[, c('interaction.pvalues', 'color')]
colnames(interaction.pvalues) <- c('pvalue', 'color')
interaction.pvalues <- interaction.pvalues[order(interaction.pvalues$pvalue), ]
interaction.pvalues$unif <- seq(1, nrow(interaction.pvalues), length.out = nrow(interaction.pvalues)) / nrow(interaction.pvalues)
interaction.pvalues$unif <- -log10(interaction.pvalues$unif)
interaction.pvalues$pvalue <- -log10(interaction.pvalues$pvalue)
 
plot <- ggplot(interaction.pvalues, aes(x = unif, y = pvalue, col = color)) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0) +
    theme_classic() +
    scale_color_manual(values = c('black' = 'black', 'red' = 'red'))

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_11_29_interaction_pvalue_ptprc_qqplot.png',
    device = 'png',
    plot = plot
)