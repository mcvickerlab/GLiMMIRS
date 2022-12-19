# plot_outlier_cell.R
#
# When analyzing significant enhancer-enhancer interactions, the significant terms seemed to have
# results dominated by a single cell with higher expression when analyzing the enhancer pair that
# impacted the histone genes. This script plots the expression of some common housekeeping genes
# and highlights the expression value of the outlier cell to see if this cell has higher overall
# expression.
#
# Author: Karthik Guruvayurappan

library(rhdf5)
library(MASS)
library(ggplot2)
library(gridExtra)

# read in cell-guide matrix
print('reading in cell-guide matrix!')
cell.guide.matrix <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'cell.guide.matrix')
guide.spacers <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'guide.spacers')
colnames(cell.guide.matrix) <- guide.spacers

# read in counts matrix
print('reading in counts matrix!')
counts.matrix <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.counts')
gene.names <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.names')
rownames(counts.matrix) <- gene.names

# compute scaling factors based on count matrix
print('computing scaling factors!')
scaling.factors <- colSums(counts.matrix) / 1e6

# variable to hold number of high expression cell
high.expression.cell <- 62383

# list of housekeeping genes published by Eisenberg et al. (2013): C1orf43, CHMP2A, EMC7, GPI, PSMB2, PSMB4, RAB7A, REEP5, SNRPD3, VCP, VPS29
# Ensembl Gene IDs: ENSG00000143612, ENSG00000130724, ENSG00000134153, ENSG00000105220, ENSG00000126067, ENSG00000159377, ENSG00000075785, ENSG00000129625, ENSG00000100028, ENSG00000165280, ENSG00000111237

# define set of housekeeping genes

housekeeping.genes.common.names <- c(
    'C1orf43', 'CHMP2A', 'EMC7', 'GPI',
    'PSMB2', 'PSMB4', 'RAB7A', 'REEP5',
    'SNRPD3', 'VCP', 'VPS29'
)

housekeeping.genes <- c(
    'ENSG00000143612', 'ENSG00000130724', 'ENSG00000134153', 'ENSG00000105220',
    'ENSG00000126067', 'ENSG00000159377', 'ENSG00000075785', 'ENSG00000129625',
    'ENSG00000100028', 'ENSG00000165280', 'ENSG00000111237'
)

# loop through housekeeping genes
gene.plots <- vector(mode = 'list', length = length(housekeeping.genes))

for (i in 1:length(housekeeping.genes)) {

    # get counts for gene and format in data frame for plotting
    gene.counts <- counts.matrix[housekeeping.genes[i], ]
    gene.counts <- data.frame(gene.counts)

    # plotting tutorial: https://felixfan.github.io/ggplot2-remove-grid-background-margin/
    # plotting tutorial: http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization 
    gene.plot <- ggplot(gene.counts, aes(x=gene.counts)) + 
                    geom_histogram(fill='black') +
                    geom_vline(
                        aes(xintercept=gene.counts[high.expression.cell]),
                        color='red',
                        linetype='dashed',
                        size=1
                    ) +
                 ggtitle(housekeeping.genes.common.names[i]) +
                 xlab('Gene UMI Count') +
                 ylab('Count') +
                 theme_bw() +
                 theme(
                     plot.title=element_text(hjust=0.5, size=16),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(colour='black'),
                 )

    gene.plots[[i]] <- gene.plot
}

combined.plot <- grid.arrange(grobs = gene.plots,
             nrow = 3,
             ncol = 4
)

ggsave(
        paste0('/iblm/netapp/home/karthik/crisprQTL/plots/', 'outlier_cell_housekeeping_counts.tiff'),
        device = 'tiff',
        plot = combined.plot,
        width = 12,
        height = 12,
        units = 'in'
)
