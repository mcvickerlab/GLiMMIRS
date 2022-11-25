# This script plots zoom-in violin plots for the 5 significant enhancer-enhancer interactions from
# the at-scale enhancer-enhancer pair data. This program was written by Karthik Guruvayurappan

library(stats)
library(rhdf5)
library(BoutrosLab.plotting.general)

# read in enhancer-enhancer p-values
enhancer.enhancer.pvalues <- read.csv(
    '/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_at_scale_20_cells_model.csv'
)
enhancer.enhancer.pvalues <- enhancer.enhancer.pvalues[complete.cases(enhancer.enhancer.pvalues), ]

# adjust interaction p-values for multiple testing correction
enhancer.enhancer.pvalues$adjusted.interaction.pvalue <- p.adjust(
    enhancer.enhancer.pvalues$interaction.pvalue.list,
    method = 'fdr'
)

# filter to only include those with low FDR
fdr.rate <- 0.1
enhancer.enhancer.pvalues <- enhancer.enhancer.pvalues[enhancer.enhancer.pvalues$adjusted.interaction.pvalue < fdr.rate, ]

# filter to only include enhancer names and gene names
enhancer.pairs <- enhancer.enhancer.pvalues[, c('enhancer.1.list', 'enhancer.2.list', 'gene.list')]
colnames(enhancer.pairs) <- c('enhancer.1', 'enhancer.2', 'gene')

# read in table mapping enhancers to spacers and reformat enhancer names
print('reading in enhancer-to-spacer table!')
enhancer.to.spacer.table <- read.table(
    '/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_grna_groups.at_scale.txt',
    sep = '\t'
)
colnames(enhancer.to.spacer.table) <- c('target.site', 'spacer.sequence')
enhancer.to.spacer.table$target.site <- sapply(enhancer.to.spacer.table$target.site, FUN = function(x) {
    if (startsWith(x, 'chr')) {
        return (strsplit(x, '_')[[1]][1])
    }
    else {
        return (x)
    }
})

guide.efficiencies.table <- h5read(
    '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5',
    'guidescan.output'
)
guide.efficiencies.table$spacer <- substring(
    guide.efficiencies.table$gRNA,
    1,
    nchar(guide.efficiencies.table$gRNA) - 3
)

# read in cell-guide matrix
print('reading in cell-guide matrix!')
cell.guide.matrix <- h5read(
    '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5',
    'cell.guide.matrix'
)
guide.spacers <- h5read(
    '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5',
    'guide.spacers'
)
colnames(cell.guide.matrix) <- guide.spacers

# read in counts matrix
print('reading in counts matrix!')
counts.matrix <- h5read(
    '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 
    'gene.counts'
)
gene.names <- h5read(
    '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5',
    'gene.names'
)
rownames(counts.matrix) <- gene.names

# compute scaling factors based on count matrix
print('computing scaling factors!')
scaling.factors <- colSums(counts.matrix) / 1e6

for (i in 1:nrow(enhancer.pairs)) {
    
    # get enhancer 1, enhancer 2, and gene name
    enhancer.1 <- enhancer.pairs[i, 'enhancer.1']
    enhancer.2 <- enhancer.pairs[i, 'enhancer.2']
    gene <- enhancer.pairs[i, 'gene']

    # get spacer sequences for each enhancer
    enhancer.1.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.1, ]$spacer.sequence
    enhancer.2.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.2, ]$spacer.sequence
    
    # get guide effiencies corresponding to spacers
    enhancer.1.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.1.spacers, c('spacer', 'Cutting.Efficiency')]
    enhancer.2.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.2.spacers, c('spacer', 'Cutting.Efficiency')]

    enhancer.1.spacers.efficiencies[is.na(enhancer.1.spacers.efficiencies)] <- 0
    enhancer.2.spacers.efficiencies[is.na(enhancer.2.spacers.efficiencies)] <- 0

    # define guide efficiency probabilities for enhancer 1
    enhancer.1.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.1.spacers.efficiencies)) {
        guide.spacer <- enhancer.1.spacers.efficiencies$spacer[j]
        guide.efficiency <- enhancer.1.spacers.efficiencies$Cutting.Efficiency[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.1.indicator.probs <- enhancer.1.indicator.probs * guide.probs
    }

    enhancer.1.indicator.probs <- 1 - enhancer.1.indicator.probs
    enhancer.1.indicator.vector <- enhancer.1.indicator.probs

    # define guide efficiency probabilities for enhancer 2
    enhancer.2.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.2.spacers.efficiencies)) {
        guide.spacer <- enhancer.2.spacers.efficiencies$spacer[j]
        guide.efficiency <- enhancer.2.spacers.efficiencies$Cutting.Efficiency[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.2.indicator.probs <- enhancer.2.indicator.probs * guide.probs
    }

    enhancer.2.indicator.probs <- 1 - enhancer.2.indicator.probs
    enhancer.2.indicator.vector <- enhancer.2.indicator.probs

    # get gene counts for the gene and normalize
    gene.counts <- counts.matrix[gene, ]
    normalized.gene.counts <- gene.counts * scaling.factors

    # create different indicator vectors for different enhancer targets
    enhancer.12.indicator.vector <- enhancer.1.indicator.vector * enhancer.2.indicator.vector
    enhancer.12.indicator.vector[enhancer.12.indicator.vector > 0] <- 1

    enhancer.1.indicator.vector[enhancer.1.indicator.vector > 0] <- 1
    enhancer.1.indicator.vector[enhancer.12.indicator.vector > 0] <- 0

    enhancer.2.indicator.vector[enhancer.2.indicator.vector > 0] <- 1
    enhancer.2.indicator.vector[enhancer.12.indicator.vector > 0] <- 0

    enhancer.all.indicator.vector <- enhancer.1.indicator.vector + enhancer.2.indicator.vector + enhancer.12.indicator.vector
    enhancer.none.indicator.vector <- 1 - enhancer.all.indicator.vector
    
    # get normalized gene counts based on indicator vector
    enhancer.1.counts <- normalized.gene.counts[enhancer.1.indicator.vector > 0]
    enhancer.2.counts <- normalized.gene.counts[enhancer.2.indicator.vector > 0]
    enhancer.12.counts <- normalized.gene.counts[enhancer.12.indicator.vector > 0]
    enhancer.none.counts <- normalized.gene.counts[enhancer.none.indicator.vector > 0]

    enhancer.1.df <- data.frame(enhancer.1.counts)
    enhancer.2.df <- data.frame(enhancer.2.counts)
    enhancer.12.df <- data.frame(enhancer.12.counts)
    enhancer.none.df <- data.frame(enhancer.none.counts)

    enhancer.1.df$group <- 'enhancer.1'
    enhancer.2.df$group <- 'enhancer.2'
    enhancer.12.df$group <- 'enhancer.12'
    enhancer.none.df$group <- 'none'

    colnames(enhancer.1.df) <- c('normalized.count', 'group')
    colnames(enhancer.2.df) <- c('normalized.count', 'group')
    colnames(enhancer.12.df) <- c('normalized.count', 'group')
    colnames(enhancer.none.df) <- c('normalized.count', 'group')

    counts.df <- rbind(enhancer.1.df, enhancer.2.df, enhancer.12.df, enhancer.none.df)
    print(min(counts.df$normalized.count))
    print(head(counts.df))
    
    create.boxplot(
        formula = normalized.count ~ group,
        data = counts.df, 
        filename = paste0('/iblm/netapp/home/karthik/crisprQTL/plots/', enhancer.1, 
                          '_', enhancer.2, '_', gene, '_zoom_in_boxplot.tiff'
                        ),
        resolution = 200,
        xlab.label = NULL,
        ylab.label = 'Normalized Count',
        main = paste(enhancer.1, enhancer.2, gene, 'expression'),
        main.cex = 1,
        ylab.cex = 1.5,
        yaxis.cex = 1,
        xaxis.cex = 0.8,
        add.stripplot = TRUE,
        points.alpha = 0.5
    )
}




