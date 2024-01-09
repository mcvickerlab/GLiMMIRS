# This program computes the number of cells with each enhancer perturbation
# and both enhancers perturbed for the at-scale enhancer-enhancer pairs.
#
# Author: Karthik Guruvayurappan

library(rhdf5)

# define h5 file name as variable
h5.name <- 'data/experimental/processed/gasperini_data.h5'

# read in enhancer-enhancer pairs
enhancer.pairs <- h5read(
    h5.name,
    'enhancer_enhancer_at_scale'
)

# read in enhancer-guide table
enhancer.guide <- h5read(
    h5.name,
    'enhancer_guide'
)

# read in cell-guide matrix
guide.matrix <- h5read(
    h5.name,
    'grna/guide_matrix'
)
guide.names <- h5read(
    h5.name,
    'grna/guide_names'
)
rownames(guide.matrix) <- guide.names
barcodes <- h5read(
    h5.name,
    'grna/cell_barcodes'
)
colnames(guide.matrix) <- barcodes

# create vectors to hold output metrics
enh.1.list <- rep(NA, nrow(enhancer.pairs))
enh.2.list <- rep(NA, nrow(enhancer.pairs))
gene.list <- rep(NA, nrow(enhancer.pairs))
enh.1.count.list <- rep(NA, nrow(enhancer.pairs))
enh.2.count.list <- rep(NA, norw(enhancer.pairs))
double.count.list <- rep(NA, nrow(enhancer.pairs))

# iterate through enhancer-enhancer pairs
for (i in 1:nrow(enhancer.pairs)) {

    # get names of both enhancers
    enhancer.1 <- enhancer.pairs$enhancer_1[i]
    enhancer.2 <- enhancer.pairs$enhancer_2[i]

    # print statement (for progress)
    print(paste('running for', enhancer.1, enhancer.2))

    # get guides for both enhancers
    enh.1.guides <- enhancer.guide[enhancer.guide$target.site == enhancer.1, ]
    enh.1.guides <- enh.1.guides$spacer

    enh.2.guides <- enhancer.guide[enhancer.guide$target.site == enhancer.2, ]
    enh.2.guides <- enh.2.guides$spacer

    # filter for spacers in the guide matrix
    enh.1.guides <- enh.1.guides[enh.1.guides %in% guide.names]
    enh.2.guides <- enh.2.guides[enh.2.guides %in% guide.names]

    # compute total perturbation vectors for each enhancer
    enh.1.perturb <- colSums(guide.matrix[enh.1.guides, ])
    enh.2.perturb <- colSums(guide.matrix[enh.2.guides, ])

    # set max perturbation to 1
    enh.1.perturb[enh.1.perturb > 0] <- 1
    enh.2.perturb[enh.2.perturb > 0] <- 1

    # compute output metrics
    gene <- enhancer.pairs$gene[i]
    enh.1.count <- sum(enh.1.perturb)
    enh.2.count <- sum(enh.2.perturb)
    double.count <- sum(enh.1.perturb * enh.2.perturb)

    # store output metrics
    enh.1.list[i] <- enhancer.1
    enh.2.list[i] <- enhancer.2
    gene.list[i] <- gene
    enh.1.count.list[i] <- enh.1.count
    enh.2.count.list[i] <- enh.2.count
    double.count.list[i] <- double.count
}

# write to output file
output <- data.frame(
    cbind(
        enh.1.list,
        enh.2.list,
        gene.list,
        enh.1.count.list,
        enh.2.count.list,
        double.count.list
    )
)
write.csv(
    output,
    'data/experimental/processed/perturbation_counts.csv',
    row.names = FALSE,
    quote = FALSE
)
