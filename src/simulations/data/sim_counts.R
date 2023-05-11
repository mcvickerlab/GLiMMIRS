# This program creates simulated count data for a Perturb-seq experiment by sampling counts from a
# negative binomial distribution, and including several experimental covariates. The outputs of
# this program include a simulated count matrix, values for the various covariates, model
# coefficients that were used for simulation, and guide metadata including efficiency and target
# gene.
# This program was written by Jessica Zhou.

library(argparse)
library(Matrix)
library(rhdf5)
library(sn)
library(dplyr)
library(tidyr)

# define parser to handle input arguments from command line
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
parser$add_argument("--guide_disp", action = "store", type = "integer", nargs = "+",
                    default = NULL,
                    help = "dispersion value(s) to use when simulating estimated guide efficiencies")
parser$add_argument("--out", action = "store", type = "character",
	                help = "where to save outputs")
args <- parser$parse_args()

##############################################################
#  write parameters of simulation to file for recordkeeping
##############################################################
args.df <- data.frame(args)
args.df$date <- Sys.Date()

write.csv(t(args.df), file.path(args$out, "simulation_params.csv"),
    row.names = TRUE, col.names = FALSE, quote = FALSE)

##############################################################
#  simulated baseline (beta0)
#  one beta per gene
#  rnorm params obtained by fitting to experimental data
##############################################################
baselines <- rnorm(args$genes, mean = 2.24, sd = 1.8)

#######################################
#  assign guides to cells
#  store as one hot encoding
#  simulate guide efficiencies 
#######################################

# compute number of guides
num.guides = args$d * args$targets

# initialize one hot encoding
onehot.guides <- matrix(0, args$cells, num.guides)

# get nr of guides per cell
guides.per.cell <- rpois(args$cells, args$lambda)

# assign guides to each cell
guides.idx.list <- sapply(guides.per.cell, function(x) {
    sample(1:num.guides, x, replace = FALSE)
    })

# update one hot encoding
for (i in 1:args$cells) {
    guides.idx <- guides.idx.list[[i]]
    onehot.guides[i, guides.idx] <- 1
}

# assign guide efficiencies to guides in library
efficiencies <- rbeta(num.guides, 6,3)

#######################################
#  simulated noisy (est.) guide efficiency
#######################################

calculate_noisy_efficiency <- function(true, D) {
    a = D*true
    b = a/true - a
    return(rbeta(1,a,b))
}

dispersions <- args$guide_disp
n.disps <- length(dispersions)

noisy.mtx <- matrix(0, length(efficiencies), n.disps)


for (i in 1:n.disps) {
    noisy.mtx[,i] <- sapply(efficiencies, calculate_noisy_efficiency, D = dispersions[i])
}

noisy.df <- data.frame(noisy.mtx)
colnames(noisy.df) <- as.character(dispersions)
noisy.df$guide <- 1:length(efficiencies)
noisy.df$true <- efficiencies

### write to file
write.table(noisy.df, file.path(args$out, "noisy_guide_efficiencies.csv"),
    row.names = TRUE, col.names = TRUE, quote = FALSE)

####################################################
#  assign target genes to gRNAs
####################################################
target.genes <- sample(1:args$genes, args$targets, replace = FALSE)

# create guide gene map which has the gene ID of the enhancer that gRNA is targeting
guide.gene.map <- rep(target.genes, args$d)

# get effect sizes of each enhancer on its target gene (beta1)
effect.sizes <- -1 * (rgamma(args$targets, shape = 6, scale = 0.5))

##############################################################
#  get vector of beta1 values (enhancer effect size)
##############################################################
beta1 <- rep(0, args$genes) # initialize vector of zeros

for (i in 1:length(target.genes)) {
    # for each target gene, update its position in beta1 to effect size
    target <- target.genes[i]
    beta1[target] <- effect.sizes[i]
}

####################################################
#  write all guide metadata to file
####################################################

# row index = gRNA identifier 
guides.metadata <- data.frame(guide.gene.map, efficiencies, rep(effect.sizes, args$d))
colnames(guides.metadata) <- c("target.gene", "efficiency", "effect.size")
write.table(guides.metadata, file.path(args$out, "guides_metadata.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate X2, X3 (S, G2M scores)
####################################################

# load cell cycle scores
covars <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5','covariates')
emp.s <- covars$s.score
emp.g2m <- covars$g2m.score

# define optimization fxn for normal
optim_norm <- function(par, data) {
    -sum(dnorm(data, mean = par[1], sd = par[2], log = TRUE))
}

# fit normal distribution to observed S scores
s.score.fit <- optim(c(1,1), optim_norm, data = emp.s, lower = c(-Inf, 1e-10), upper = c(Inf, Inf), method = "L-BFGS-B")

# simulate S scores
s.scores <- rnorm(args$cells, mean = s.score.fit$par[1], sd = s.score.fit$par[2])

# fit skew normal to observed G2M scores
g2m.fit <- selm(X ~ 1, data=data.frame(X=emp.g2m))
extractSECdistr(g2m.fit) #extract parameters from fitting skew normal to G2M scores

# simulate G2M scores (hardcoded based on values from fitting empirical data to skn in a notebook)
g2m.scores <- rsn(n=args$cells, xi=-0.2556359, omega=0.3124325, alpha=6.2932919, tau=0, dp=NULL)

# write cell cycle scores (X2, X3) to file
# row index = cell identifier
cell.cycle.scores <- data.frame(s.scores, g2m.scores)
write.table(cell.cycle.scores, file.path(args$out, "cell_cycle_scores.txt"), row.names = TRUE, quote = FALSE)

#################################################################
#  simulate beta2, beta3 (from gamma distr.) for S, G2M scores
#################################################################
beta2 <- rgamma(args$genes, shape = 6, scale = 0.5)
beta3 <- rgamma(args$genes, shape = 6, scale = 0.5)

####################################################
#  simulate beta4, x4 for percent.mito
####################################################
percent.mito <- rbeta(args$cells, shape1 = 3.3, shape2 = 81.48)
beta4 <- rgamma(args$genes , shape = 6 , scale = 0.5)

# write percent.mito to file
# row index = cell identifier
percent.mito.df <- data.frame(percent.mito)
write.table(percent.mito.df, file.path(args$out, "percent_mito.txt"), row.names = TRUE, quote = FALSE)

##############################################################
#  write coeffs to file (beta0, beta1, beta2, beta3, beta4)
##############################################################

# row index = gene identifier
coeffs <- data.frame(beta0=baselines, 
                    beta.enh=beta1, 
                    beta.s=beta2, 
                    beta.g2m=beta3, 
                    beta.mito=beta4)
write.table(coeffs, file.path(args$out, "coeffs.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate scaling factors (poisson)
####################################################

# simulate total counts per cell
t.vec <- rpois(args$cells, 50000)

# get scaling factors
scaling.factors <- t.vec / 1e6

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
# represents enhancer perturbation probability
x1.mtx <- matrix(0, args$genes, args$cells)

# initialize matrices for storing values of linear predictor and mu
lp.mtx <- matrix(0, args$genes, args$cells)

# populate counts mtx by gene (by row)
for (gene in 1:args$genes) {
    cat(sprintf("simulating counts for gene %d\n", gene))

    # get coeffs
    b0 <- baselines[gene]
    b1 <- beta1[gene]
    b2 <- beta2[gene]
    b3 <- beta3[gene]
    b4 <- beta4[gene]
    
    # initialize x1 as vector of zeros (assume it is not affected by any gRNAs in library)
    x1 <- numeric(args$cells)
    
    # check if enhancer of gene is targeted by any gRNAs
    if (gene %in% target.genes) {
        cat(sprintf("gene %d is a targeting gene\n", gene))
        # calculate perturbation probability (X1)
        guides.for.gene <- which(guide.gene.map==gene)
        temp.mtx <- t(efficiencies[guides.for.gene]*t(onehot.guides[,guides.for.gene]))
        x1 <- apply(temp.mtx, 1, function(x) {1-prod(1-x)})
        cat(sprintf("sum of perturbation probabilities = %.3f\n", sum(x1)))

        x1.mtx[gene,] <- x1
    }

    # get cell cycle scores
    x2 <- s.scores
    x3 <- g2m.scores
    # get percent.mito
    x4 <- percent.mito

    # calculate values of mu 
    lp.vec <- b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + log(scaling.factors)
    lp.mtx[gene,] <- lp.vec

    # use rnbinom to generate counts of each cell for this gene and update counts matrix
    counts <- rnbinom(length(mu.vec), mu = mu.vec, size = 1.5)
    sim.counts[gene,] <- counts
}


####################################################
#  write all data to h5
####################################################
print("writing data to h5")
h5.path <- file.path(args$out, "sim.h5")
h5createFile(h5.path)

# write counts, chunked
print('writing counts, chunked')
h5createGroup(h5.path, "counts")
h5createDataset(h5.path, "counts/counts", dim(sim.counts),
    storage.mode = "integer", chunk=c(1000, 10000), level=5)

h5write(sim.counts, h5.path, "counts/counts")

# linear predictor, chunked
print('writing linear predictor')
h5createGroup(h5.path, "linear_predictor")
h5createDataset(h5.path, "linear_predictor/linear_predictor", dim(lp.mtx),
    storage.mode = "double", chunk=c(1000, 10000), level=5)

h5write(lp.mtx, h5.path, "linear_predictor/linear_predictor")

# write guide info
print('writing guide info')
h5createGroup(h5.path, "guides")
h5createDataset(h5.path, "guides/one_hot", dim(onehot.guides),
    storage.mode = "integer", chunk=c(1000, 1), level=5)

h5write(onehot.guides, h5.path,"guides/one_hot")
h5write(guides.metadata, h5.path, "guides/metadata")

# write estimate guide efficiencies
if (!is.null(args$guide_disp)) {
    print('writing estimated guide efficiencies')
    h5write(noisy.df, h5.path, "guides/noisy_guide_efficiencies")
}

# write coeffs
print('writing coeffs')
h5write(coeffs, h5.path, "coeffs")

# write perturbation_probability
print('writing perturbation probability')
h5createGroup(h5.path, "x")
h5createDataset(h5.path, "x/perturbation_prob", dim(x1.mtx),
    storage.mode = "double", chunk = c(1000, 1000), level = 7)
h5write(x1.mtx, h5.path, "x/perturbation_prob")

# write cell cycle scores
print('writing cell cycle scores')
h5write(cell.cycle.scores, h5.path, "x/cell_cycle_scores")

# write percent.mito
print('writing percent.mito')
h5write(percent.mito, h5.path, "x/percent_mito")

# write scaling factors
print('writing scaling factors')
h5write(scaling.factors, h5.path, "scaling_factors")
