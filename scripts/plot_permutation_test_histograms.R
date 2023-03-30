# To determine p-values associated with the significant interaction
# coefficients that were sensitive to outliers, we ran 100 permutation tests
# for each enhancer pair and recorded the interaction coefficients. This
# code plots the distribution of those coefficients, along with where the
# observed interaction coefficient falls on that distribution.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# read in output results from at-scale enhancer-enhancer analysis
at.scale.results <- read.csv(
    '/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv'
)
at.scale.results <- at.scale.results[complete.cases(at.scale.results), ]

# add FDR adjusted p-values and filter for FDR < 0.1
at.scale.results$adjusted.interaction.pvalue <- p.adjust(at.scale.results$interaction.pvalue.list, method = 'fdr')
significant.interactions <- at.scale.results[at.scale.results$adjusted.interaction.pvalue < 0.1, ]

# filter for necessary columns in significant interactions
significant.interactions <- significant.interactions[, c('enhancer.1.list', 'enhancer.2.list', 'gene.list', 'interaction.coeff.list')]
colnames(significant.interactions) <- c('enhancer.1', 'enhancer.2', 'gene', 'interaction.coeff')

permutation.plots <- vector(mode = 'list', length = nrow(significant.interactions))

for (i in 1:nrow(significant.interactions)) {

    # read in enhancer and gene names
    enhancer.1 <- significant.interactions[i, 'enhancer.1']
    enhancer.2 <- significant.interactions[i, 'enhancer.2']
    gene <- significant.interactions[i, 'gene']
    interaction.coefficient <- significant.interactions[i, 'interaction.coeff']
    coeff.df <- data.frame(interaction.coefficient)
    colnames(coeff.df) <- c('coeff')

    # read in file with null distribution of interaction coefficients
    null.coeffs <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/gasperini_data/23_03_29_chr1.12443_chr1.12449_ENSG00000116903_null_interaction_coefficient_estimates.csv')
    colnames(null.coeffs) <- c('coeffs')

    # plot null distribution of coefficients as a histogram
    plot <- ggplot(null.coeffs, aes(x = coeffs)) +
        geom_histogram(color = 'black') +
        geom_vline(xintercept = interaction.coefficient, color = 'red', linetype = 'dashed') +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab(bquote("Interaction Coefficient")) + 
        ylab(bquote(Count)) +
        theme_classic() +
        theme(
            axis.line = element_line(linewidth = 1),
            axis.title.x = element_text(size = 20, color = 'black'),
            axis.title.y = element_text(size = 20, color = 'black'),
            axis.text = element_text(size = 20, color = 'black'),
            axis.ticks = element_line(color = 'black', linewidth = 1),
            axis.ticks.length = unit(2, 'mm'),
            plot.margin = rep(unit(10, 'mm'), 4),
            legend.position = 'none'
    )

    permutation.plots[[i]] <- plot
}

combined.plot <- grid.arrange(grobs = permutation.plots,
             nrow = 2,
             ncol = 2
)

ggsave(
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/plots/', '23_03_29_permutation_histograms.pdf'),
        device = 'pdf',
        plot = combined.plot,
        width = 12,
        height = 12,
        units = 'in'
)

ggsave(
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/plots/', '23_03_29_permutation_histograms.png'),
        device = 'png',
        plot = combined.plot,
        width = 12,
        height = 12,
        units = 'in'
)
