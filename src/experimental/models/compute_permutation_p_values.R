# This script uses the FastQTL p-value strategy by fitting a beta distribution
# to p-values obtained from 100 permutation tests, and then estimating a more
# stringent p-value based on that.
#
# Author: Karthik Guruvayurappan

library(fitdistrplus)
library(ggplot2)

# get list of files from permutation tests
permutation.files <- list.files(
    'data/experimental/processed/permutation_test_results/'
)

enhancer.1.list <- rep(NA, length(permutation.files))
enhancer.2.list <- rep(NA, length(permutation.files))
gene.list <- rep(NA, length(permutation.files))
shape.1.list <- rep(NA, length(permutation.files))
shape.2.list <- rep(NA, length(permutation.files))

# iterate through files and fit beta distributions
for (i in 1:length(permutation.files)) {

    # get permutations file
    file.name <- permutation.files[i]

    # read in file and p-values
    p.values <- read.csv(
        paste0(
            'data/experimental/processed/permutation_test_results/',
            file.name
        )
    )
    p.values <- p.values$permutation.interaction.pvalues

    # fit a beta distribution
    fit <- fitdist(
        p.values,
        'beta'
    )
    print(fit$convergence)

    # store relevant info
    enhancer.1 <- strsplit(file.name, '_')[[1]][1]
    enhancer.2 <- strsplit(file.name, '_')[[1]][2]
    gene <- strsplit(file.name, '_')[[1]][3]
    shape.1 <- fit$estimate[[1]]
    shape.2 <- fit$estimate[[2]]

    enhancer.1.list[i] <- enhancer.1
    enhancer.2.list[i] <- enhancer.2
    gene.list[i] <- gene
    shape.1.list[i] <- shape.1
    shape.2.list[i] <- shape.2
}

# merge results into data frame
beta.df <- data.frame(cbind(
    enhancer.1.list,
    enhancer.2.list,
    gene.list,
    shape.1.list,
    shape.2.list
))

# get observed p-value estimate
at.scale.results <- data.frame()

for (i in 1:32) {
    batch.results <- read.csv(
        paste0('data/experimental/processed/enhancer_pairs_at_scale_', i, '.csv')
    )
    at.scale.results <- rbind(at.scale.results, batch.results)
}

# merge beta distribution results with at-scale results
beta.df <- merge(beta.df, at.scale.results)
beta.df$shape.1.list <- as.numeric(beta.df$shape.1.list)
beta.df$shape.2.list <- as.numeric(beta.df$shape.2.list)

# compute adjusted p-values
beta.p.values <- rep(NA, nrow(beta.df))

for (i in 1:nrow(beta.df)) {
    shape.1 <- beta.df$shape.1.list[i]
    shape.2 <- beta.df$shape.2.list[i]
    true.pvalue <- beta.df$interaction.pvalues[i]
    beta.p.value <- pbeta(true.pvalue, shape.1, shape.2)
    beta.p.values[i] <- beta.p.value
}

beta.df$permutation.pvalues <- beta.p.values

beta.df <- beta.df[
    ,
    c(
        'enhancer.1.list',
        'enhancer.2.list',
        'gene.list',
        'interaction.pvalues',
        'permutation.pvalues'
    )
]

# take negative log10 of p-values
beta.df$log.interaction.pvalues <- -log10(beta.df$interaction.pvalues)
beta.df$log.permutation.pvalues <- -log10(beta.df$permutation.pvalues)

# scatterplot
plot <- ggplot(
    beta.df,
    aes(x = log.interaction.pvalues, y = log.permutation.pvalues)
    ) +
    geom_point() +
    theme_classic()

ggsave(
    filename = 'out/24_01_31_permutation_pvalues_scatter.png',
    device = 'png',
    plot = plot
)

