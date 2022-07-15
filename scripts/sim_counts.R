library(argparse)
library(Matrix)
library(mvtnorm)
library(ggplot2)
library(hexbin)

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
	help = "number of gRNAs for each target")
parser$add_argument("--lambda", action = "store", type = "double",
	default = 15,
	help = "lambda parameter for Poisson distribution to draw number of guides per cell")
parser$add_argument("--out", action = "store", type = "character",
	help = "where to save outputs")

args <- parser$parse_args()

nGuides = args$d * args$targets

#######################################
#  simulated baseline beta0
#######################################
baselines <- rnorm(args$genes, mean = 2.24, sd = 1.8)

png(file.path(args$out, "hist_beta0.png"))
hist(baselines)
dev.off()

#######################################
#  assign guides to cells
#######################################

# initialize one hot encoding
onehot.guides <- matrix(0, args$cells, nGuides)

# get guides per cell
guides.per.cell <- rpois(args$cells, args$lambda)

png(file.path("hist_guides_per_cell.png"))
hist(guides.per.cell, xlab = "guides per cell", ylab = "cells")  
dev.off()

# assign guides to each cell
guides.idx.list <- sapply(guides.per.cell, function(x) {sample(1:nGuides, x, replace = FALSE)})

# update one hot encoding
for (i in 1:args$cells) {
    guides.idx <- guides.idx.list[[i]]
    onehot.guides[i, guides.idx] <- 1
}

# write cell by guides mapping to sparse matrix
writeMM(Matrix(onehot.guides, sparse = TRUE), file.path(args$out, "guides_per_cell.mtx"))

# assign guide efficiencies to guides in library
efficiencies <- rbeta(nGuides, 6,3)

png(file.path(args$out, "hist_guide_efficiencies.png"))
hist(efficiencies)
dev.off()

####################################################
#  assign target genes to candidate enhancers
####################################################
target.genes <- sample(1:args$genes, args$targets, replace = FALSE)
guide.gene.map <- rep(target.genes, args$d)

# get effect sizes of each enhancer
effect.sizes <- -1*(rgamma(args$targets, shape = 6, scale = 0.5))

####################################################
#  get vector of beta1 values
####################################################
beta1 <- numeric(args$genes)

for (i in 1:length(target.genes)) {
    target <- target.genes[i]
    beta1[target] <- effect.sizes[i]
}

png(file.path(args$out, "hist_beta1.png"))
hist(beta1)
dev.off()

####################################################
#  define function for calculating X1 (per cell)
####################################################

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
#  define all guide metadata to file
####################################################
write.table(data.frame(guide.gene.map, efficiencies, rep(effect.sizes,2)), 
	file.path(args$out, "guides_metadata.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  use MVN to generate cell cycle scores (X2,X3)
####################################################
cov.mtx <- matrix(c(1, -0.8,-0.8,1), ncol = 2)

rmv.rand <- rmvnorm(n=args$cells, mean = c(0,0), sigma = cov.mtx)

s.scores <- rmv.rand[,1]
g2m.scores <- rmv.rand[,2]

# plot as hex
png(file.path(args$out, "cell_cycle_scores_hexbin.png"))
bin<-hexbin(s.scores, g2m.scores, xbins=50)
plot(bin, main="Hexagonal Binning")
dev.off()

####################################################
#  select beta2, beta3
####################################################
beta2 <- rgamma(args$genes, shape = 6, scale = 0.5)
beta3 <- rgamma(args$genes, shape =6, scale = 0.5)

png(file.path(args$out, "beta2_hist.png"))
hist(beta2)
dev.off()

png(file.path(args$out, "beta3_hist.png"))
hist(beta3)
dev.off()

####################################################
#  scaling factors
####################################################
t.vec <- rpois(args$cells, 50000)

png(file.path(args$out, "t_pois_hist.png"))
hist(t.vec)
dev.off()

scaling.factors <- t.vec/1e6

png(file.path(args$out, "scaling_factors_hist.png"))
hist(scaling.factors)
dev.off()

####################################################
#  simulate data
####################################################

# initialize counts matrix
sim.counts <-  matrix(0, args$genes, args$cells)

# populate counts mtx by gene (by row)
for (gene in 1:args$genes) {
    # get coeffs
    b0 <- baselines[gene]
    b1 <- beta1[gene]
    b2 <- beta2[gene]
    b3 <- beta3[gene]
    
    # initialize x1 as vector of zeros (assume it will not be affected by any guides in library)
    x1 <- numeric(args$cells)
    
    # check if enhancer of gene is targeted by any guides in our library
    if (gene %in% target.genes) {
        for (i in 1:args$cells) {
            x1[i] <- combined_prob(i, gene)
        }
    }
    
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
