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
        '/iblm/netapp/home/karthik/GLiMMIRS/data/gasperini/processed/sting_seq_bootstrap/24_01_05_sting_seq_bootstrap_interactions_',
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

    plot.df <- data.frame(cbind(adj.intercept, adj.enhancer.1, adj.enhancer.2, adj.interaction))

    # create plot
    plot.means <- colMeans(plot.df)
    plot.lower <- apply(plot.df, 2, quantile, probs = 0.05)
    plot.upper <- apply(plot.df, 2, quantile, probs = 0.95)
    plot.df <- data.frame(t(rbind(plot.means, plot.lower, plot.upper)))
    rownames(plot.df) <- c(
        'None',
        'E1',
        'E2',
        'E1 + E2'
    )
    plot.df$perturbation <- rownames(plot.df)
    plot.df$perturbation <- factor(plot.df$perturbation, levels=c('None', 'E1', 'E2', 'E1 + E2'))

    plot <- ggplot(plot.df, aes(x = perturbation, y = plot.means)) +
      geom_rect(
          aes(xmin = -Inf, xmax = Inf, ymin = quantile(adj.expected, 0.05), ymax = quantile(adj.expected, 0.95)),
          fill = 'darkgray',
          alpha = 0.6,
          color = 'black',
          linetype = 'dashed'
        ) +
      geom_pointrange(aes(ymin = plot.lower, ymax = plot.upper), fatten = 10, linewidth = 1) +
      ylab('Scaled Expression') +
      theme_classic() +
      xlab('Perturbation') +
      theme(
        axis.line = element_line(linewidth = 1),
        axis.title.x = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        axis.text = element_text(size = 20, color = 'black', family = 'Helvetica'),
        axis.ticks = element_line(color = 'black', linewidth = 1),
        axis.ticks.length = unit(2, 'mm'),
        plot.margin = rep(unit(10, 'mm'), 4)
      )

    # write to output file
    ggsave(
        paste0(
            '/iblm/netapp/home/karthik/GLiMMIRS/out/sting_seq_bootstrap/',
            snp.1,
            '_',
            snp.2,
            '.pdf'
        ),
        plot = plot,
        device = 'pdf',
        width = 7,
        height = 7,
        unit = 'in'
    )
}
