# This script plots the concordance between the GLiMMIRS baseline model and
# the interaction model to test if the coefficients are changing in magnitude
# in the interaction model.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# read in baseline model results
baseline.results <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_11_19_enhancer_gene_pairs_suppl_table_2_baseline_pseudocount_model.csv')
baseline.results <- baseline.results[, c('enhancer.list', 'gene.list', 'enhancer.effect.list')]

# read in results from the 330 enhancer-enhancer pairs
interaction.results <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_10_31_enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv')

# isolate results from enhancer 1 from interaction models and merge
enhancer.1.int.results <- interaction.results[, c('enhancer.1.list', 'gene.list', 'enhancer.1.coeff.list')]

colnames(baseline.results) <- c('enhancer', 'gene', 'baseline.effect')
colnames(enhancer.1.int.results) <- c('enhancer', 'gene', 'interaction.effect')
enhancer.1.concordance <- merge(enhancer.1.int.results, baseline.results, by = c('enhancer', 'gene'), all.x = TRUE)

enhancer.2.int.results <- interaction.results[, c('enhancer.2.list', 'gene.list', 'enhancer.2.coeff.list')]
colnames(enhancer.2.int.results) <- c('enhancer', 'gene', 'interaction.effect')
enhancer.2.concordance <- merge(enhancer.2.int.results, baseline.results, by = c('enhancer', 'gene'), all.x = TRUE)

plot <- ggplot(enhancer.1.concordance, aes(x = baseline.effect, y = interaction.effect)) +
    geom_point() +
    theme_classic() +
    xlab('Enhancer 1 Baseline Effect') +
    ylab('Enhancer 1 Interaction Effect')

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_11_20_enhancer_1_concordance.png',
    device = 'png', 
    plot = plot
)

plot <- ggplot(enhancer.2.concordance, aes(x = baseline.effect, y = interaction.effect)) +
    geom_point() +
    theme_classic() +
    xlab('Enhancer 2 Baseline Effect') +
    ylab('Enhancer 2 Interaction Effect')

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_11_20_enhancer_2_concordance.png',
    device = 'png', 
    plot = plot
)
