# This script plots estimates of the posterior interaction frequency given some
# prior interaction frequency and various degrees of power in our models,
# combined with the fact that we did not observe a significant interaction
# across all of our tests.
#
# Authors: Graham McVicker & Karthik Guruvayurappan

library(ggplot2)

# This function parameterizes the mean and standard deviation of a beta
# distribution into alpha and beta terms.
calc_alpha_beta <- function(mean_val, sd_val) {
  alpha <- ((1 - mean_val) / (sd_val**2) - (1/mean_val)) * mean_val**2
  beta <- alpha * (1/mean_val - 1)

  return (list(mean=mean_val, sd=sd_val, alpha=alpha, beta=beta))
}

# This function parameterizes the alpha and beta terms of a beta distribution
# into mean and standard deviation terms
calc_mean_sd <- function(alpha_val, beta_val) {
  mean_val <- alpha_val / (alpha_val + beta_val)
  sd_val <- sqrt((alpha_val * beta_val) / ((alpha_val + beta_val)**2) * (alpha_val + beta_val + 1))

  return(list(mean=mean_val, sd=sd_val, alpha=alpha_val, beta=beta_val))
}

# define prior interaction frequencies
prior_interaction_frequencies <- seq(0.05, 1, length.out = 20)
prior_interaction_sds <- rep(
  0.025,
  length(prior_interaction_frequencies)
)

# create data frame to hold all results
results_df <- data.frame()

for (i in 1:length(prior_interaction_frequencies)) {

  # get current prior interaction frequency + sd
  prior_interaction_frequency <- prior_interaction_frequencies[i]
  prior_interaction_sd <- prior_interaction_sds[i]

  # reparameterize in terms of alpha and beta
  prior_params <- calc_alpha_beta(
    prior_interaction_frequency,
    prior_interaction_sd
  )

  # create distribution of different powers
  prior_results_df <- data.frame(
    power = c(
      0.05,
      0.1,
      0.2, 
      0.4,
      0.6,
      0.8
    )
  )
  prior_results_df$prior <- prior_interaction_frequency

  # define number of pairs tested for interactions
  num_tests <- 264

  # define number of significant interactions
  num_signif <- 0

  # iterate through different powers
  for (j in 1:nrow(prior_results_df)) {

    # get current power
    power <- prior_results_df$power[j]

    # multiply prior of interaction by power of detection
    # this is the frequency with which you will detect true interactions
    prior_detect_mean <- power * prior_interaction_frequency
    prior_detect_sd <- power * prior_interaction_sd

    # reparameterize in terms of alpha and beta
    prior_detect_params <- calc_alpha_beta(prior_detect_mean, prior_detect_sd)

    # calculate posterior given observed interactions
    post_detect_alpha <- prior_detect_params$alpha + num_signif
    post_detect_beta <- num_tests - num_signif + prior_detect_params$beta
    post_detect_params <- calc_mean_sd(post_detect_alpha, post_detect_beta)

    # rescale mean of posterior by power
    prior_results_df$posterior_mean[j] <- post_detect_params$mean / power
  }

  # append to larger results data frame
  results_df <- rbind(results_df, prior_results_df)
}

# plot as line graph
results_df$power <- as.factor(results_df$power)

plot <- ggplot(results_df, aes(x = prior, y = posterior_mean, group = power, color = power)) +
  geom_line(linewidth = 1) +
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab('Prior Interaction Frequency') +
  ylab('Posterior Interaction Frequency') +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 24, color = 'black'),
    axis.title.y = element_text(size = 24, color = 'black'),
    axis.text = element_text(size = 24, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    legend.position = c(0.15, 0.80),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  labs(color = 'Power') +
  scale_color_brewer(palette = 'Dark2')

ggsave(
  plot = plot,
  filename = 'out/high_confidence_bayesian_interaction_detection.pdf',
  device = 'pdf'
)
