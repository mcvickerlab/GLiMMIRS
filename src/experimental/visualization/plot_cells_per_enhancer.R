# This script plots the numbers of cells with perturbations to each of the
# enhancers in an enhancer pair, in addition to the number of cells with both
# and no enhancers perturbed. This code is based on the at-scale analysis,
# which aimed to test 4,666 gene loci.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# read in enhancer pairs with counts of each perturbation
enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_08_03_at_scale_enhancer_enhancer_pairs_cells_counts.csv')
matched.enhancers <- enhancer.pairs$enhancer.1.list == enhancer.pairs$enhancer.2.list
enhancer.pairs <- enhancer.pairs[!matched.enhancers, ]
enhancer.pairs <- enhancer.pairs[enhancer.pairs$count.list > 20, ]

# separate counts for enhancer 1 perturbation
enhancer.1.perturbations <- data.frame(enhancer.pairs$enhancer.1.count.list)
enhancer.1.perturbations$perturb <- 'enhancer 1'
colnames(enhancer.1.perturbations) <- c('count', 'perturb')

# separate counts for enhancer 2 perturbation
enhancer.2.perturbations <- data.frame(enhancer.pairs$enhancer.2.count.list)
enhancer.2.perturbations$perturb <- 'enhancer 2'
colnames(enhancer.2.perturbations) <- c('count', 'perturb')

# separate counts for both perturbations
both.perturbations <- data.frame(enhancer.pairs$count.list)
both.perturbations$perturb <- 'enhancer 1 + enhancer 2'
colnames(both.perturbations) <- c('count', 'perturb')

# create a violin plot and save to output file
plot.df <- rbind(enhancer.1.perturbations, enhancer.2.perturbations, both.perturbations)

plot <- ggplot(plot.df, aes(x = perturb, y = count)) +
    geom_boxplot() +
    coord_flip() +
    theme_classic()

ggsave(
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/out/', '23_08_04_perturbation_violin.png'),
        device = 'png',
        plot = plot,
        width = 6,
        height = 6,
        units = 'in'
)
