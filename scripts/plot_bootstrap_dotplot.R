# This program plots dotplots with error bars for the bootstrap confidnece
# interval estimates computed for the significant enhancer-enhancer interaction
# terms.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
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
significant.interactions <- significant.interactions[, c('enhancer.1.list', 'enhancer.2.list', 'gene.list')]
colnames(significant.interactions) <- c('enhancer.1', 'enhancer.2', 'gene')

# keep vector of bootstrap plots
bootstrap.plots <- vector(mode = 'list', length = nrow(significant.interactions))

print(significant.interactions)
significant.interactions$gene.hgnc <- c('EXOC8', 'BABAM2', 'H2BC12', 'ZBED9')

for (i in 1:nrow(significant.interactions)) {

    # define enhancers and gene
    enhancer.1 <- significant.interactions[i, 'enhancer.1']
    enhancer.2 <- significant.interactions[i, 'enhancer.2']
    gene <- significant.interactions[i, 'gene']
    gene.hgnc <- significant.interactions[i, 'gene.hgnc']

    # read in bootstrap interaction coefficient estimates
    bootstrap.estimates.filename <- paste0('/iblm/netapp/data1/external/Gasperini2019/processed/23_02_23_', enhancer.1, '_', enhancer.2, '_', gene, '_bootstrap_coefficient_estimates.csv')
    bootstrap.interaction.estimates <- read.csv(
        bootstrap.estimates.filename
    )
    colnames(bootstrap.interaction.estimates) <- c('coefficient')
    bootstrap.interaction.estimates$name <- ''
    print(quantile(bootstrap.interaction.estimates$coefficient, probs = c(0.005, 0.995)))

    low.percentile <- function(x) {
        return(quantile(x, probs = c(0.005)))
    }

    high.percentile <- function(x) {
        return(quantile(x, probs = c(0.995)))
    }

    mid.percentile <- function(x) {
        return(quantile(x, probs = c(0.5)))
    }

    bootstrap.dotplot <- ggplot(bootstrap.interaction.estimates, aes(x=name, y=coefficient)) +
                         geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.75, color = 'black', fill = 'black') +
                         # geom_jitter() +
                         theme_classic() +
                         theme(
                            axis.line = element_line(linewidth = 1),
                            axis.title.x = element_text(size = 14, color = 'black'),
                            axis.title.y = element_text(size = 10, color = 'black'),
                            axis.text = element_text(size = 14, color = 'black'),
                            axis.ticks = element_line(color = 'black', linewidth = 1),
                            axis.ticks.length = unit(2, 'mm'),
                            plot.margin = rep(unit(5, 'mm'), 4),
                            plot.title = element_text(size = 10, color = 'black', hjust = 0.5),
                        ) +
                         geom_hline(yintercept = 0, linetype = 'dashed') +
                         stat_summary(fun.min = low.percentile, fun = mid.percentile, fun.max = high.percentile, geom="pointrange", color="red", size = 0.75, linewidth = 1.2) +
                         ggtitle(paste(enhancer.1, enhancer.2, gene.hgnc)) +
                         xlab('') +
                         ylab('Interaction Coefficient')

    bootstrap.plots[[i]] <- bootstrap.dotplot
}

combined.plot <- grid.arrange(grobs = bootstrap.plots,
             nrow = 2,
             ncol = 2
)

ggsave(
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/plots/', '23_04_20_bootstrap_dotplot.pdf'),
        device = 'pdf',
        plot = combined.plot,
        width = 6,
        height = 6,
        units = 'in'
)

ggsave(
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/plots/', '23_04_20_bootstrap_dotplot.png'),
        device = 'png',
        plot = combined.plot,
        width = 6,
        height = 6,
        units = 'in'
)
