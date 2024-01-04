# There are two significant interactions from the STING-seq data. Here, we
# plot expression estimates derived from fitted models to see how the true
# data compares to our baseline multiplicative expectation.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# read in bootstrap estimates
bootstrap.estimates <- read.csv(
    '/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_bootstrap_interactions_rs1926231_rs1326279.csv'
)

intercept <- bootstrap.estimates$true.null.estimates
enhancer.1 <- bootstrap.estimates$true.enhancer.1.estimates
enhancer.2 <- bootstrap.estimates$true.enhancer.2.estimates
interaction <- bootstrap.estimates$true.interaction.estimates

adj.intercept <- exp(intercept)
adj.enhancer.1 <- exp(intercept + enhancer.1)
adj.enhancer.2 <- exp(intercept + enhancer.2)
adj.interaction <- exp(intercept + enhancer.1 + enhancer.2 + interaction)
adj.expected <- exp(intercept + enhancer.1 + enhancer.2)

plot.df <- data.frame(cbind(adj.intercept, adj.enhancer.1, adj.enhancer.2, adj.interaction, adj.expected))

plot.means <- colMeans(plot.df)
plot.sds <- apply(plot.df, 2, sd)
plot.ses <- plot.sds / sqrt(nrow(plot.df))
plot.lower <- apply(plot.df, 2, quantile, probs = 0.05)
plot.upper <- apply(plot.df, 2, quantile, probs = 0.95)



plot.df <- data.frame(t(rbind(plot.means, plot.sds, plot.lower, plot.upper)))
rownames(plot.df) <- c(
    'No Perturbation',
    'E1',
    'E2',
    'E1 + E2',
    'E1 + E2 (expected)'
)
plot.df$perturbation <- rownames(plot.df)
plot.df$perturbation <- factor(plot.df$perturbation, levels=c('No Perturbation', 'E1', 'E2', 'E1 + E2', 'E1 + E2 (expected)'))
plot <- ggplot(plot.df, aes(x = perturbation, y = plot.means)) +
    geom_point() +
    geom_errorbar(aes(ymin = plot.lower, ymax = plot.upper)) +
    ylab('Scaled Expression') +
    theme_classic()

ggsave(
    '/iblm/netapp/home/karthik/GLiMMIRS/out/test_expectation_plot.png',
    plot = plot,
    device = 'png'
)



