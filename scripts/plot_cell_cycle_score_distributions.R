# This program plots the distribution of cell cycle scores output by Seurat which will be
# incorporated into the linear model.
# This program was written by Karthik Guruvayurappan.

library(BoutrosLab.plotting.general)

s.scores <- read.csv('/iblm/netapp/home/karthik/crisprQTL/gasperini_data/s_scores.csv')$S.Score
create.histogram(
    x = s.scores,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/s_score_distribution.tiff',
    resolution = 200,
    ylimits = c(0, 23),
    xlab.label = 'S Score'
)

g2m.scores <- read.csv('/iblm/netapp/home/karthik/crisprQTL/gasperini_data/g2m_scores.csv')$G2M.Score
create.histogram(
    x = g2m.scores,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/g2m_score_distribution.tiff',
    resolution = 200,
    ylimits = c(0, 28),
    xlab.label = 'G2M Score'
)
