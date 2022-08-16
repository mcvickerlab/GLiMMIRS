# This program plots the distribution of guide efficiencies calculated using the GuideScan 2.0 tool.
# This program was written by Karthik Guruvayurappan.

library(BoutrosLab.plotting.general)

# read in guidescan output
guidescan.guide.info <- read.csv('/iblm/netapp/home/karthik/GuideScan/Gasperini2019/guidescan_output.csv')

create.histogram(
    x = guidescan.guide.info[, c('Cutting.Efficiency')],
    filename = file.path('/iblm/netapp/home/karthik/crisprQTL/plots/guide_efficiency_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Guide Efficiency',
    ylab.label = 'Count',
    type = 'count',
    xlimits = c(0, 1),
    ylimits = c(0, 2600),
    xat = seq(0, 1, 0.2),
    yat = seq(0, 2500, 500)
)
