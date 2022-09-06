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

s.scores <- read.csv('/iblm/netapp/home/karthik/crisprQTL/gasperini_data/s_scores.csv')
g2m.scores <- read.csv('/iblm/netapp/home/karthik/crisprQTL/gasperini_data/g2m_scores.csv')
cell.cycle.scores <- merge(s.scores, g2m.scores, by = 'X')

create.scatterplot(
    formula = G2M.Score ~ S.Score,
    data = cell.cycle.scores,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/cell_s_score_vs_g2m_score.tiff',
    resolution = 200,
    xlab.label = 'S Score',
    ylab.label = 'G2M Score',
    alpha = 0.2
)
