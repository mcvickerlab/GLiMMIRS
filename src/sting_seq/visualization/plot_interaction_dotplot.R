# There are two significant interactions from the STING-seq data. Here, we
# plot expression estimates derived from fitted models to see how the true
# data compares to our baseline multiplicative expectation.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# create list of PTPRC pairs
# read in supplementary table S1E
table.s1e <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s1e.csv'
)
ptprc.snps <- table.s1e$rs.ID

ptprc.pairs <- data.frame(t(combn(ptprc.snps, 2)))
colnames(ptprc.pairs) <- c('snp.1', 'snp.2')

# iterate through PTPRC pairs
for (i in 1:nrow(ptprc.pairs)) {

    snp.1 <- ptprc.pairs$snp.1[i]
    snp.2 <- ptprc.pairs$snp.2[i]
    print(paste('plotting', snp.1, snp.2))


    # read in file with bootstrap estimates
    boot.file <- paste0(
        '/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_bootstrap/24_01_05_sting_seq_bootstrap_interactions_',
        snp.1,
        '_',
        snp.2,
        '.csv'
    )
    bootstrap.estimates <- read.csv(boot.file)

    # compute scaled linear predictor results and expectations
    intercept <- bootstrap.estimates$true.null.estimates
    enhancer.1 <- bootstrap.estimates$true.enhancer.1.estimates
    enhancer.2 <- bootstrap.estimates$true.enhancer.2.estimates
    interaction <- bootstrap.estimates$true.interaction.estimates

    adj.intercept <- exp(intercept)
    adj.enhancer.1 <- exp(intercept + enhancer.1)
    adj.enhancer.2 <- exp(intercept + enhancer.2)
    adj.interaction <- exp(intercept + enhancer.1 + enhancer.2 + interaction)
    adj.expected <- exp(intercept + enhancer.1 + enhancer.2)
    adj.additive <- exp(intercept) - (exp(intercept) - exp(intercept + enhancer.1)) - (exp(intercept) - exp(intercept + enhancer.2))

    plot.df <- data.frame(cbind(adj.intercept, adj.enhancer.1, adj.enhancer.2, adj.interaction, adj.expected, adj.additive))

    # create plot
    plot.means <- colMeans(plot.df)
    plot.lower <- apply(plot.df, 2, quantile, probs = 0.05)
    plot.upper <- apply(plot.df, 2, quantile, probs = 0.95)
    plot.df <- data.frame(t(rbind(plot.means, plot.lower, plot.upper)))
    rownames(plot.df) <- c(
        'No Perturbation',
        'E1',
        'E2',
        'E1 + E2',
        'E1 + E2 (mult)',
        'E1 + E2 (add)'
    )
    plot.df$perturbation <- rownames(plot.df)
    plot.df$perturbation <- factor(plot.df$perturbation, levels=c('No Perturbation', 'E1', 'E2', 'E1 + E2', 'E1 + E2 (mult)', 'E1 + E2 (add)'))

    plot <- ggplot(plot.df, aes(x = perturbation, y = plot.means)) +
        geom_point() +
        geom_errorbar(aes(ymin = plot.lower, ymax = plot.upper)) +
        ylab('Scaled Expression') +
        theme_classic()

    # write to output file
    ggsave(
        paste0(
            '/iblm/netapp/home/karthik/GLiMMIRS/out/sting_seq_bootstrap/',
            snp.1,
            '_',
            snp.2,
            '.png'
        ),
        plot = plot,
        device = 'png'
    )
}



