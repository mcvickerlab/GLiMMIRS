# To determine p-values associated with the significant interaction
# coefficients that were sensitive to outliers, we ran 100 permutation tests
# for each enhancer pair and recorded the interaction coefficients. This
# code plots the distribution of those coefficients, along with where the
# observed interaction coefficient falls on that distribution.
#
# Author: Karthik Guruvayurappan

source('src/experimental/models/analysis_helpers.R')

library(ggplot2)

# create directory to hold output files
dir.create('out/adaptive_permutation_histograms/')

# get significant interactions from at-scale analysis
significant_interactions <- get_significant_results()

# iterate through significant interactions and plot permutation histograms
for (i in 1:nrow(significant_interactions)) {

  # get enhancer 1, enhancer 2, and gene name
  enhancer_1 <- significant_interactions[i, 'enhancer.1.list']
  enhancer_2 <- significant_interactions[i, 'enhancer.2.list']
  gene <- significant_interactions[i, 'gene.list']

  # store value of observed interaction coefficient
  observed_interaction <- significant_interactions[i, 'interaction.effects']

  # store filename for permutations and read in permutation results
  permutation_results_filename <- paste0(
    'data/experimental/processed/adaptive_permutation_test_results/',
    enhancer_1,
    '_',
    enhancer_2,
    '_',
    gene,
    '_adaptive_permutations.csv'
  )

  # read in file containing permutation test results
  permutation_results <- read.csv(permutation_results_filename)

  # plot histograms of permutation coefficients
  plot <- ggplot(permutation_results, aes(x = interaction.effects)) +
    geom_histogram(color = 'black') +
    theme_classic() +
    geom_vline(
      xintercept = observed_interaction,
      color = 'red',
      linetype = 'dashed'
    ) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab('Interaction Coefficient') +
    ylab('Count') +
    ggtitle(paste(enhancer_1, enhancer_2, gene)) +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.title.x = element_text(size = 14, color = 'black'),
      axis.title.y = element_text(size = 14, color = 'black'),
      axis.text = element_text(size = 10, color = 'black'),
      axis.ticks = element_line(color = 'black', linewidth = 1),
      axis.ticks.length = unit(2, 'mm'),
      plot.margin = rep(unit(10, 'mm'), 4),
      plot.title = element_text(size = 8.5, color = 'black', hjust = 0.5),
    )

  ggsave(
    filename = paste0(
      'out/adaptive_permutation_histograms/',
      enhancer_1,
      '_',
      enhancer_2,
      '_',
      gene,
      '.png'
    ),
    device = 'png',
    plot = plot
  )
  

}

















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
    null.coeffs <- read.csv(paste0('/iblm/netapp/home/karthik/GLiMMIRS/gasperini_data/23_03_31_', enhancer.1, '_', enhancer.2, '_', gene, '_null_interaction_coefficient_estimates.csv'))
    colnames(null.coeffs) <- c('coeffs')

    print(paste0('pvalue: ', mean(abs(null.coeffs) > interaction.coefficient)))

    # plot null distribution of coefficients as a histogram
    plot <- ggplot(null.coeffs, aes(x = coeffs)) +
        geom_histogram(color = 'black') +
        geom_vline(xintercept = interaction.coefficient, color = 'red', linetype = 'dashed') +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab(bquote("Interaction Coefficient")) + 
        ylab(bquote(Count)) +
        ggtitle(paste(enhancer.1, enhancer.2, gene)) +
        theme_classic() +
        theme(
            axis.line = element_line(linewidth = 1),
            axis.title.x = element_text(size = 14, color = 'black'),
            axis.title.y = element_text(size = 14, color = 'black'),
            axis.text = element_text(size = 10, color = 'black'),
            axis.ticks = element_line(color = 'black', linewidth = 1),
            axis.ticks.length = unit(2, 'mm'),
            plot.margin = rep(unit(10, 'mm'), 4),
            legend.position = 'none',
            plot.title = element_text(size = 8.5, color = 'black', hjust = 0.5),
    )

    permutation.plots[[i]] <- plot
}

combined.plot <- grid.arrange(grobs = permutation.plots,
             nrow = 2,
             ncol = 2
)

ggsave(
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/plots/', '23_04_24_permutation_histograms.pdf'),
        device = 'pdf',
        plot = combined.plot,
        width = 6,
        height = 6,
        units = 'in'
)

ggsave(
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/plots/', '23_04_24_permutation_histograms.png'),
        device = 'png',
        plot = combined.plot,
        width = 6,
        height = 6,
        units = 'in'
)
