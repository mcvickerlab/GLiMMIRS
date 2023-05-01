# This program compares the distance between enhancer pairs and the interaction coefficients, for
# both the at-scale enhancer pairs and 330 enhancer pairs. (pseudocount model)
# This program was written by Karthik Guruvayurappan.

library(ggplot2)
library(stats)

# read in 330 interaction coefficients and enhancer distances
coeffs.distances <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_04_25_enhancer_distance_330_pairs.csv')
coeffs.distances <- coeffs.distances[complete.cases(coeffs.distances), ]

# fit loess curve
loess.curve <- loess(formula = interaction.coeff.list ~ distance, data = coeffs.distances)
line.df <- data.frame(cbind(coeffs.distances$distance, loess.curve$fitted))
line.df <- line.df[order(line.df$X1), ]

scatter.plot <- ggplot(coeffs.distances, aes(x = distance, y = interaction.coeff.list)) +
    geom_point(aes(size = 3), alpha = 0.2) +
    geom_smooth(method = 'loess', se = FALSE) + 
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Distance)) + 
    ylab(bquote(Coefficient)) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    legend.position = "none"
    )

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_25_distance_coefficient_330_pairs.pdf',
    device = 'pdf',
    plot = scatter.plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_25_distance_coefficient_330_pairs.png',
    device = 'png',
    plot = scatter.plot
)

# read in 3808 interaction coefficients and enhancer distances
coeffs.distances <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_04_25_enhancer_distance_3808_pairs.csv')
coeffs.distances <- coeffs.distances[complete.cases(coeffs.distances), ]

scatter.plot <- ggplot(coeffs.distances, aes(x = distance, y = interaction.coeff.list)) +
    geom_point(aes(size = 3), alpha = 0.2) +
    geom_smooth(method = 'loess', se = FALSE) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Distance)) + 
    ylab(bquote(Coefficient)) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    legend.position = "none"
    )

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_25_distance_coefficient_3808_pairs.pdf',
    device = 'pdf',
    plot = scatter.plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_25_distance_coefficient_3808_pairs.png',
    device = 'png',
    plot = scatter.plot
)
