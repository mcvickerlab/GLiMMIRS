# This script splits the 46 significant interactions from the at-scale analysis
# into 32 different files for parallelization in the permutation testing
# procedure.
#
# Author: Karthik Guruvayurappan

# read in at-scale enhancer pair analysis results
model.results <- data.frame()

for (i in 1:32) {
    batch.file <- paste0(
        'data/experimental/processed/enhancer_pairs_at_scale_',
        i,
        '.csv'
    )
    batch.results <- read.csv(batch.file)
    model.results <- rbind(model.results, batch.results)
}

# filter for results with valid guide efficiency info
model.results <- model.results[complete.cases(model.results), ]

# compute FDR-adjusted p-values and filter
model.results$adj.interaction.pvalues <- p.adjust(
    model.results$interaction.pvalues,
    method = 'fdr'
)
significant.results <- model.results[
    model.results$adj.interaction.pvalues < 0.1,

]

# split significant interactions into 32 batches
rownames(significant.results) <- NULL
batches <- as.integer(rownames(significant.results)) %% 32
significant.results <- split(significant.results, batches)

# write out batch significant interaction files
for (i in 1:32) {
    write.csv(
        significant.results[[i]],
        paste0('data/experimental/interim/significant_results_', i, '.csv'),
        row.names = FALSE,
        quote = FALSE
    )
}
