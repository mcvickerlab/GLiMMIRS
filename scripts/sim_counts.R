library(argparse)
library(Matrix)
library(mvtnorm)
library(ggplot2)
library(hexbin)
library(rhdf5)

###### outputs
# simulated counts matrix
# x1, x2, x3
# beta0, beta1, beta2, beta3
# one hot encoding of guides x cells
# guide metadata (target gene, efficiency, effect size)
# scaling factors

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--genes', action = "store", type = "integer", 
	default = 13000,
	help = "number of genes in simulated dataset")
parser$add_argument("--targets", action = "store", type = "integer",
	default = 1000,
	help = "number of candidate enhancers targeted in simulation")
parser$add_argument("--cells", action = "store", type = 'integer',
	default = 50000,
	help = "number of cells in simualted dataset")
parser$add_argument("--d", action = "store", type = "integer",
	default = 2,
	help = "number of gRNAs for each target site (candidate enhancer)")
parser$add_argument("--lambda", action = "store", type = "double",
	default = 15,
	help = "lambda parameter for Poisson distribution to draw number of guides per cell")
parser$add_argument("--out", action = "store", type = "character",
	help = "where to save outputs")

args <- parser$parse_args()

nGuides = args$d * args$targets

##############################################################
#  simulated baseline (beta0)
#  one beta per gene
#  rnorm params obtained by fitting to experimental data
##############################################################
baselines <- rnorm(args$genes, mean = 2.24, sd = 1.8)

# plot simulated values of beta0 (baseline)
png(file.path(args$out, "hist_beta0.png"))
hist(baselines)
dev.off()

#######################################
#  assign guides to cells
#  store as one hot encoding
#  simulate guide efficiencies 
#######################################

# initialize one hot encoding
onehot.guides <- matrix(0, args$cells, nGuides)

# get # of guides per cell
guides.per.cell <- rpois(args$cells, args$lambda)

# plot # of guides per cell
png(file.path(args$out, "hist_guides_per_cell.png"))
hist(guides.per.cell, xlab = "guides per cell", ylab = "cells")  
dev.off()

# assign guides to each cell
guides.idx.list <- sapply(guides.per.cell, 
    function(x) {sample(1:nGuides, x, replace = FALSE)})

# update one hot encoding
for (i in 1:args$cells) {
    guides.idx <- guides.idx.list[[i]]
    onehot.guides[i, guides.idx] <- 1
}

# write cell by guides mapping to sparse matrix
writeMM(Matrix(onehot.guides, sparse = TRUE), 
    file.path(args$out, "guides_per_cell.mtx"))

# assign guide efficiencies to guides in library
efficiencies <- rbeta(nGuides, 6,3)

# plot guide efficiencies
png(file.path(args$out, "hist_guide_efficiencies.png"))
hist(efficiencies)
dev.off()

####################################################
#  assign target genes to gRNAs
####################################################
target.genes <- sample(1:args$genes, args$targets, replace = FALSE)
guide.gene.map <- rep(target.genes, args$d)

# get effect sizes of each enhancer on its target gene 
effect.sizes <- -1*(rgamma(args$targets, shape = 6, scale = 0.5))

####################################################
#  get vector of beta1 values
####################################################
beta1 <- numeric(args$genes) # initialize vector of zeros

for (i in 1:length(target.genes)) {
    # for each target gene, update its position in beta1 to effect size
    target <- target.genes[i]
    beta1[target] <- effect.sizes[i]
}

# plot beta1 values (expected pileup at zero)
png(file.path(args$out, "hist_beta1.png"))
hist(beta1)
dev.off()

############################################################################
#  define a function for calculating X1 (different values for each gene) 
############################################################################ 

combined_prob <- function(cell, gene, verbose = FALSE) {
    if (verbose) {
        cat(sprintf("calculating value of X1 for gene %d in cell %d\n", gene, cell))
    }
    # identify which gRNAs in our design target this gene
    guides <- which(guide.gene.map==gene)
    
    # check if any of these gRNAs are present in cell
    if (sum(onehot.guides[cell, guides]) > 0) {
        terms <- numeric(length(guides))
        for (i in 1:args$d) {
            if (onehot.guides[cell,guides[i]]!=0) {
                terms[i] <- 1-efficiencies[guides[i]]
            }
        }
        x1 <- rbinom(1,1,1-prod(1-terms))
        return(x1)
    } else {
        return(0)
    }
}

####################################################
#  write all guide metadata to file
####################################################

# row index = gRNA identifier 
guides.metadata <- data.frame(guide.gene.map, efficiencies, rep(effect.sizes,2))
colnames(guides.metadata) <- c("target.gene", "efficiency", "effect.size")
write.table(guides.metadata, file.path(args$out, "guides_metadata.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  use MVN to generate cell cycle scores (X2,X3)
####################################################
cov.mtx <- matrix(c(1, -0.8,-0.8,1), ncol = 2)

rmv.rand <- rmvnorm(n=args$cells, mean = c(0,0), sigma = cov.mtx)

s.scores <- rmv.rand[,1]
g2m.scores <- rmv.rand[,2]

# plot S/G2M scores as hexbin
png(file.path(args$out, "cell_cycle_scores_hexbin.png"))
bin<-hexbin(s.scores, g2m.scores, xbins=50)
plot(bin, main="Hexagonal Binning")
dev.off()

# write cell cycle scores (X2, X3) to file
# row index = cell identifier
cell.cycle.scores <- data.frame(s.scores, g2m.scores)
write.table(cell.cycle.scores, file.path(args$out, "cell_cycle_scores.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate beta2, beta3 (from gamma distr.)
####################################################
beta2 <- rgamma(args$genes, shape = 6, scale = 0.5)
beta3 <- rgamma(args$genes, shape =6, scale = 0.5)

# plot beta2
png(file.path(args$out, "beta2_hist.png"))
hist(beta2)
dev.off()

# plot beta3
png(file.path(args$out, "beta3_hist.png"))
hist(beta3)
dev.off()

######################################################
#  write coeffs to file (beta0, beta1, beta2, beta3)
######################################################

# row index = gene identifier
coeffs <- data.frame(baselines, beta1, beta2, beta3)
write.table(coeffs, file.path(args$out, "coeffs.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate scaling factors (poisson)
####################################################

# simulate total counts per cell
t.vec <- rpois(args$cells, 50000)

# plot total counts per cell
png(file.path(args$out, "t_pois_hist.png"))
hist(t.vec)
dev.off()

# get scaling factors
scaling.factors <- t.vec/1e6

# plot scaling factors 
png(file.path(args$out, "scaling_factors_hist.png"))
hist(scaling.factors)
dev.off()

# write scaling factors to file
# row index = cell identifier
write.table(data.frame(scaling.factors), file.path(args$out, "scaling_factors.txt"), 
    row.names = TRUE, quote = FALSE)

####################################################
#  simulate counts matrix 
####################################################

# initialize counts matrix
sim.counts <-  matrix(0, args$genes, args$cells)

# initialize matrix of x1 values (different set of values for each gene)
x1.mtx <- matrix(0, args$genes, args$cells)

# populate counts mtx by gene (by row)
for (gene in 1:args$genes) {
    # get coeffs
    b0 <- baselines[gene]
    b1 <- beta1[gene]
    b2 <- beta2[gene]
    b3 <- beta3[gene]
    
    # initialize x1 as vector of zeros (assume it is not affected by any gRNAs in library)
    x1 <- numeric(args$cells)
    
    # check if enhancer of gene is targeted by any gRNAs
    if (gene %in% target.genes) {
        for (i in 1:args$cells) {
            x1[i] <- combined_prob(i, gene)
        }
    }
    
    x1.mtx[gene,] <- x1

    # get cell cycle scores
    x2 <- s.scores
    x3 <- g2m.scores
    
    # calculate values of mu
    mu.vec <- scaling.factors*exp(b0 + b1*x1 + b2*x2 + b3*x3)
    
    # use rnbinom to generate counts of each cell for this gene and update counts matrix
#     counts <- sapply(mu.vec, function(x) {rnbinom(1, mu = x, size = 1.5)})
    counts <- rnbinom(length(mu.vec), mu = mu.vec, size = 1.5)
    sim.counts[gene,] <- counts
}

# write simulated counts to sparse matrix
writeMM(Matrix(sim.counts), file.path(args$out, "counts.mtx"))

# write matrix of x1 values to sparse matrix 
writeMM(Matrix(x1.mtx), file.path(args$out, "x1.mtx"))

####################################################
#  write all data to h5
####################################################
h5.path <- file.path(args$out, "sim.h5")
h5createFile(h5.path)

# write counts, chunked
h5createDataset(h5.path, "counts", dim(sim.counts),
    storage.mode = "integer", chunk=c(1000, 10000), level=5)

h5write(sim.counts, h5.path, "counts")

# write guide info
h5createGroup(h5.path, "guides")
h5createDataset(h5.path, "guides/one_hot", dim(onehot.guides),
    storage.mode = "integer", chunk=c(1000, 1), level=5)

h5write(onehot.guides, h5.path,"guides/one_hot")
h5write(guides.metadata, h5.path, "guides/metadata")

# write coeffs
h5write(coeffs, h5.path, "coeffs")

# write x1
h5createGroup(h5.path, "x")

h5createDataset(h5.path, "x/x1", dim(x1.mtx),
    storage.mode = "integer", chunk = c(1000,1000), level = 7)
h5write(x1.mtx, h5.path, "x/x1")

# write cell cycle scores
h5write(cell.cycle.scores, h5.path, "x/cell_cycle_scores")

# write scaling factors
h5write(scaling.factors, h5.path, "scaling_factors")
