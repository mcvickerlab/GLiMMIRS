# This program plots dotplots with error bars for the bootstrap confidnece
# interval estimates computed for the significant enhancer-enhancer interaction
# terms.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(gridExtra)

# read in output results from at-scale enhancer-enhancer analysis
at.scale.results <- read.csv(
    '/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv'
)
at.scale.results <- at.scale.results[complete.cases(at.scale.results), ]

# add FDR adjusted p-values and filter for FDR < 0.1
at.scale.results$adjusted.interaction.pvalue <- p.adjust(at.scale.results$interaction.pvalue.list, method = 'fdr')
significant.interactions <- at.scale.results[at.scale.results$adjusted.interaction.pvalue < 0.1, ]

# filter for necessary columns in significant interactions
significant.interactions <- significant.interactions[, c('enhancer.1.list', 'enhancer.2.list', 'gene.list')]
colnames(significant.interactions) <- c('enhancer.1', 'enhancer.2', 'gene')

# keep vector of bootstrap plots
bootstrap.plots <- vector(mode = 'list', length = nrow(significant.interactions))

for (i in 1:nrow(significant.interactions)) {

    # define enhancers and gene
    enhancer.1 <- significant.interactions[i, 'enhancer.1']
    enhancer.2 <- significant.interactions[i, 'enhancer.2']
    gene <- significant.interactions[i, 'gene']

    # read in bootstrap interaction coefficient estimates
    bootstrap.estimates.filename <- paste0('/iblm/netapp/data1/external/Gasperini2019/processed/22_12_21_', enhancer.1, '_', enhancer.2, '_', gene, '_bootstrap_coefficient_estimates.csv')
    bootstrap.interaction.estimates <- read.csv(
        bootstrap.estimates.filename
    )
    colnames(bootstrap.interaction.estimates) <- c('coefficient')
    bootstrap.interaction.estimates$name <- ''

    bootstrap.dotplot <- ggplot(bootstrap.interaction.estimates, aes(x=name, y=coefficient)) +
                         geom_dotplot(binaxis='y', stackdir='center') +
                         theme_bw() +
                         theme(
                             plot.title=element_text(hjust=0.5, size=11),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank(),
                             panel.border = element_blank(),
                             axis.line = element_line(colour='black'),
                         ) + 
                         stat_summary(fun.data=mean_sdl, fun.args = list(mult=2), geom="pointrange", color="blue") +
                         ggtitle(paste(enhancer.1, enhancer.2, gene)) +
                         xlab('') +
                         ylab('Interaction Coefficient')

    bootstrap.plots[[i]] <- bootstrap.dotplot
}

combined.plot <- grid.arrange(grobs = bootstrap.plots,
             nrow = 2,
             ncol = 3
)

ggsave(
        paste0('/iblm/netapp/home/karthik/crisprQTL/plots/', 'bootstrap_dotplot.tiff'),
        device = 'tiff',
        plot = combined.plot,
        width = 12,
        height = 12,
        units = 'in'
)
