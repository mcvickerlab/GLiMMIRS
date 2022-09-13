# This program creates simulated count data for a Perturb-seq experiment by sampling counts from a
# negative binomial distribution, and including several experimental covariates. The outputs of
# this program include a simulated count matrix, values for the various covariates, model
# coefficients that were used for simulation, and guide metadata including efficiency and target
# gene.
# This program was written by Jessica Zhou.

library(argparse)
library(Matrix)
library(mvtnorm)
library(ggplot2)
library(hexbin)
library(rhdf5)
library(BoutrosLab.plotting.general)
library(sn)

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
# parser$add_argument("--data", action = "store", type = "character",
#                     help = "file with empirical data info")
args <- parser$parse_args()

##############################################################
#  write parameters of simulation to file for recordkeeping
##############################################################
args.df <- data.frame(args)
args.df$date <- Sys.Date()

write.csv(t(args.df), file.path(args$out, "simulation_params.csv"))

##############################################################
#  simulated baseline (beta0)
#  one beta per gene
#  rnorm params obtained by fitting to experimental data
##############################################################
baselines <- rnorm(args$genes, mean = 2.24, sd = 1.8)

# plot simulated values of beta0 (baseline)
create.histogram(
    x = baselines,
    filename = file.path(args$out, 'beta0_distribution.tiff'),
    resolution = 200,
    xlab.label = 'beta0 value',
    ylab.label = 'Percent'
)

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

# plot nr of guides per cell
create.histogram(
    x = guides.per.cell,
    filename = file.path(args$out, 'guides_per_cell_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Guides in cell',
    ylab.label = 'Percent'
)

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

# plot guide efficiencies
create.histogram(
    x = efficiencies,
    filename = file.path(args$out, 'guide_efficiency_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Guide efficiency',
    ylab.label = 'Percent of guides'
)

#######################################
#  simulated est. guide efficiency
#######################################
est.efficiencies.list <- list()
disps.list <- list()
i <- 1

if (!is.null(args$guide_disp)) {
    print("simulating estimated guide efficiency values")

    for (gd in args$guide_disp) {
        cat(sprintf("D=%d\n",gd))
        est.eff <- rbeta(num.guides, 
            efficiencies*gd, 
            efficiencies*gd)  
        est.efficiencies.list[[i]] <- est.eff 
        disps.list[[i]] <- rep(gd, num.guides)
        # writeLines(est.efficiencies.list, file.path(args$out, sprintf("est_efficiencies_D%d.txt", d)))

        # tiff(file.path(args$out, 
        #     sprintf("hist_est_guide_efficiences_D%d.tiff", gd)))
        png(file.path(args$out, sprintf("hist_est_guide_efficiencies_D%d.png", gd)))
        hist(est.eff, 
            main = sprintf("Histogram of estimated guide efficiencies, D=%d", gd))
        dev.off()

        i <- i + 1
    }
} 

est.efficiencies.df <- data.frame(est.efficiency = do.call(c, est.efficiencies.list), 
                                    D = do.call(c, disps.list))

# plot estimated efficiencies together on sme histogram
est.efficiencies.df$D <- as.factor(est.efficiencies.df$D)
est.efficiencies.p <- ggplot(est.efficiencies.df, 
                            aes(x = est.efficiency, fill = D, color = D)) + 
                        geom_histogram(position = "dodge") + 
                        theme_classic() + 
                        theme(text = element_text(size = 20))

# tiff(file.path(args$out, 
#     "hist_est_guide_efficiences.tiff"), 
#         res = 300, units = "in")
png(file.path(args$out, "hist_est_guide_efficiencies.png"))
print(est.efficiencies.p)
dev.off()

# head(est.efficiencies.df)
write.table(est.efficiencies.df, file.path(args$out, "est_guide_efficiencies.csv"),
    row.names = TRUE, col.names = TRUE, quote = FALSE)

####################################################
#  assign target genes to gRNAs
####################################################
target.genes <- sample(1:args$genes, args$targets, replace = FALSE)

# create guide gene map which has the gene ID of the enhancer that gRNA is targeting
guide.gene.map <- rep(target.genes, args$d)

# get effect sizes of each enhancer on its target gene (beta1)
effect.sizes <- -1 * (rgamma(args$targets, shape = 6, scale = 0.5))

####################################################
#  get vector of beta1 values
####################################################
beta1 <- rep(0, args$genes) # initialize vector of zeros

for (i in 1:length(target.genes)) {
    # for each target gene, update its position in beta1 to effect size
    target <- target.genes[i]
    beta1[target] <- effect.sizes[i]
}

# plot beta1 values (expected pileup at zero)
create.histogram(
    x = beta1,
    filename = file.path(args$out, 'beta1_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Beta 1 value',
    ylab.label = 'Percent'
)

####################################################
#  write all guide metadata to file
####################################################

# row index = gRNA identifier 
guides.metadata <- data.frame(guide.gene.map, efficiencies, rep(effect.sizes, 2))
colnames(guides.metadata) <- c("target.gene", "efficiency", "effect.size")
write.table(guides.metadata, file.path(args$out, "guides_metadata.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate X2, X3 (cell cycle scores)
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
extractSECdistr(g2m.fit)

# simulate G2M scores (hardcoded based on values in notebook)
g2m.scores <- rsn(n=args$cells, xi=-0.2556359, omega=0.3124325, alpha=6.2932919, tau=0, dp=NULL)

# Plot 
# tiff("../rnorm_s_scores.tiff",
#     res = 300, units = "in")
png(file.path(args$out, "rnorm_s_scores.png"))
hist(s.scores, 
    main = expression(
        paste("Simulated S scores, ", 
            mu, "=-1.296e-3, ", 
            sigma, "=0.11")))
dev.off()

# tiff("../rsn_g2m_scores.tiff",
#     res = 300, units = "in")
png(file.path(args$out, "rsn_g2m_scores.png"))
hist(g2m.scores,
    main = expression(
        paste("Simulated G2M scores\n", 
            xi, "=-0.256, ", 
            omega, "=0.312, ", 
            alpha, "=6.29, ", 
            tau, "=0")))
dev.off()


# write cell cycle scores (X2, X3) to file
# row index = cell identifier
cell.cycle.scores <- data.frame(s.scores, g2m.scores)
write.table(cell.cycle.scores, file.path(args$out, "cell_cycle_scores.txt"), row.names = TRUE, quote = FALSE)

# ####################################################
# #  use MVN to generate cell cycle scores (X2,X3)
# ####################################################
# cov.mtx <- matrix(c(1, -0.8,-0.8,1), ncol = 2)

# rmv.rand <- rmvnorm(n=args$cells, mean = c(0,0), sigma = cov.mtx)

# s.scores <- rmv.rand[,1]
# g2m.scores <- rmv.rand[,2]

# # plot S/G2M scores as hexbin
# png(file.path(args$out, "cell_cycle_scores_hexbin.png"))
# bin < -hexbin(s.scores, g2m.scores, xbins = 50)
# plot(bin, main="Hexagonal Binning")
# dev.off()

# # write cell cycle scores (X2, X3) to file
# # row index = cell identifier
# cell.cycle.scores <- data.frame(s.scores, g2m.scores)
# write.table(cell.cycle.scores, file.path(args$out, "cell_cycle_scores.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate beta2, beta3 (from gamma distr.)
####################################################
beta2 <- rgamma(args$genes, shape = 6, scale = 0.5)
beta3 <- rgamma(args$genes, shape = 6, scale = 0.5)

# plot beta2
create.histogram(
    x = beta2,
    filename = file.path(args$out, 'beta2_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Beta 2 value',
    ylab.label = 'Percent'
)

# plot beta3
create.histogram(
    x = beta3,
    filename = file.path(args$out, 'beta3_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Beta 3 value',
    ylab.label = 'Percent'
)

####################################################
#  simulate beta4, x4
####################################################
percent.mito <- rbeta(args$cells, shape1 = 3.3, shape2 = 81.48)

# plot distribution of percent mito
create.histogram(
    x = percent.mito,
    filename = file.path(args$out, 'percent_mito_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Percent mito value',
    ylab.label = 'Percent'
)

beta4 <- rgamma(args$genes , shape = 6 , scale = 0.5)

# plot
create.histogram(
    x = beta4,
    filename = file.path(args$out, 'beta4_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Beta 4 value',
    ylab.label = 'Percent'
)

# write percent.mito to file
# row index = cell identifier
percent.mito.df <- data.frame(percent.mito)
write.table(percent.mito.df, file.path(args$out, "percent_mito.txt"), row.names = TRUE, quote = FALSE)

######################################################
#  write coeffs to file (beta0, beta1, beta2, beta3)
######################################################

# row index = gene identifier
coeffs <- data.frame(baselines, beta1, beta2, beta3, beta4)
write.table(coeffs, file.path(args$out, "coeffs.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate scaling factors (poisson)
####################################################

# simulate total counts per cell
t.vec <- rpois(args$cells, 50000)

# plot total counts per cell
create.histogram(
    x = t.vec,
    filename = file.path(args$out, 'counts_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Counts per cell',
    ylab.label = 'Percent'
)

# get scaling factors
scaling.factors <- t.vec / 1e6

# plot distribution of scaling factors
create.histogram(
    x = scaling.factors,
    filename = file.path(args$out, 'scaling_factor_distribution.tiff'),
    resolution = 200,
    xlab.label = 'Scaling factor',
    ylab.label = 'Percent'

)

# write scaling factors to file
# row index = cell identifier
write.table(data.frame(scaling.factors), file.path(args$out, "scaling_factors.txt"), 
                       row.names = TRUE, quote = FALSE)


############################################################################
#  define a function for calculating X1 (different values for each gene) 
############################################################################ 

# combined_prob <- function(cell, gene, verbose = FALSE) {
#     # calculate X1 as a binary indicator variable
#     if (verbose) {
#         cat(sprintf("calculating value of X1 for gene %d in cell %d\n", gene, cell))
#     }

#     # identify which gRNAs in our design target this gene
#     guides <- which(guide.gene.map==gene)
    
#     # check if any of these gRNAs are present in cell
#     if (sum(onehot.guides[cell, guides]) > 0) {
#         terms <- numeric(length(guides))
#         for (i in 1:args$d) {
#             if (onehot.guides[cell,guides[i]]!=0) {
#                 terms[i] <- 1-efficiencies[guides[i]]
#             }
#         }
#         x1 <- rbinom(1,1,1-prod(1-terms))
#         return(x1)
#     } else {
#         return(0)
#     }
# }

# combined_prob <- function(cell, gene, verbose = FALSE) {
#     # calculate X1 as a combined probability 
#     if (verbose) {
#         cat(sprintf("calculating value of X1 for gene %d in cell %d\n", gene, cell))
#     }

#     # identify which gRNAs in our design target this gene
#     guides <- which(guide.gene.map==gene)
    
#     terms <- numeric(length(guides))
#     for (i in 1:length(guides)) {
#         cat(sprintf("i = %d\n", i))
#         if (onehot.guides[cell,guides[i]]!=0) {
#         # if (sum(h5read(args$h5, "guides/one_hot", index = list(cell, guides[i])))!=0) {
#             terms[i] <- efficiencies[guides[i]]
#         }
#     }
#     x1 <- 1-prod(1-terms)
#     return(x1)
# }


####################################################
#  simulate counts matrix 
####################################################

# initialize counts matrix
sim.counts.discrete <-  matrix(0, args$genes, args$cells)
sim.counts.continuous <-  matrix(0, args$genes, args$cells)

# initialize matrix of x1 values (different set of values for each gene)
x1.discrete.mtx <- matrix(0, args$genes, args$cells)
x1.continuous.mtx <- matrix(0, args$genes, args$cells)

# initialize matrices for storing values of linear predictor and mu
lp.discrete.mtx <- matrix(0, args$genes, args$cells)
lp.continuous.mtx <- matrix(0, args$genes, args$cells)
mu.discrete.mtx <- matrix(0, args$genes, args$cells)
mu.continuous.mtx <- matrix(0, args$genes, args$cells)

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
    x1.discrete <- numeric(args$cells)
    x1.continuous <- numeric(args$cells)
    
    # check if enhancer of gene is targeted by any gRNAs
    if (gene %in% target.genes) {
        cat(sprintf("gene %d is a targeting gene\n", gene))
        guides.for.gene <- which(guide.gene.map==gene)
        temp.mtx <- t(efficiencies[guides.for.gene]*t(onehot.guides[,guides.for.gene]))
        x1.continuous <- apply(temp.mtx, 1, function(x) {1-prod(1-x)})
        x1.discrete <- rbinom(args$cells,1,x1.continuous)
        # for (i in 1:args$cells) {
        #     x1[i] <- combined_prob(i, gene)
        # }
        # cat(sprintf("x1 total = %.3f\n", sum(x1)))

        x1.continuous.mtx[gene,] <- x1.continuous
        x1.discrete.mtx[gene,] <- x1.discrete
    }

    # get cell cycle scores
    x2 <- s.scores
    x3 <- g2m.scores
    x4 <- percent.mito
    
    # calculate values of mu - with discrete x
    lp.discrete.vec <- b0 + b1*x1.discrete + b2*x2 + b3*x3 + b4*x4 + log(scaling.factors)
    lp.discrete.mtx[gene,] <- lp.discrete.vec
    mu.discrete.vec <- exp(lp.discrete.vec)
    mu.discrete.mtx[gene,] <- mu.discrete.vec

    # calculate values of mu - with continuous x
    lp.continuous.vec <- b0 + b1*x1.continuous + b2*x2 + b3*x3 + b4*x4 + log(scaling.factors)
    lp.continuous.mtx[gene,] <- lp.continuous.vec
    mu.continuous.vec <- exp(lp.continuous.vec)
    mu.continuous.mtx[gene,] <- mu.continuous.vec

    # use rnbinom to generate counts of each cell for this gene and update counts matrix
#     counts <- sapply(mu.vec, function(x) {rnbinom(1, mu = x, size = 1.5)})
    counts.discrete <- rnbinom(length(mu.discrete.vec), mu = mu.discrete.vec, size = 1.5)
    sim.counts.discrete[gene,] <- counts.discrete
    counts.continuous <- rnbinom(length(mu.continuous.vec), mu = mu.continuous.vec, size = 1.5)
    sim.counts.continuous[gene,] <- counts.continuous
}

# # write simulated counts to sparse matrix
# print('writing simulated counts to sparse matrix')
# writeMM(Matrix(sim.counts), file.path(args$out, "counts.mtx"))

# # write matrix of x1 values to sparse matrix 
# print('writing x1 to sparse matrix')
# writeMM(Matrix(x1.mtx), file.path(args$out, "x1.mtx"))

# # write linear predictor/mu values to sparse matrix
# print('writing linear predictors to MM file')
# print(head(lp.mtx))
# writeMM(Matrix(lp.mtx), file.path(args$out, "linear_predictor.mtx"))

# print('writing mu to MM file')
# print(head(mu.mtx))
# writeMM(Matrix(mu.mtx), file.path(args$out, "mu.mtx"))

####################################################
#  write all data to h5
####################################################
print("writing data to h5")
h5.path <- file.path(args$out, "sim.h5")
h5createFile(h5.path)

# write counts, chunked
print('writing counts, chunked')
h5createGroup(h5.path, "counts")
h5createDataset(h5.path, "counts/discrete", dim(sim.counts.discrete),
    storage.mode = "integer", chunk=c(1000, 10000), level=5)

h5write(sim.counts.discrete, h5.path, "counts/discrete")

h5createDataset(h5.path, "counts/continuous", dim(sim.counts.continuous),
    storage.mode = "integer", chunk=c(1000, 10000), level=5)

h5write(sim.counts.continuous, h5.path, "counts/continuous")

# linear predictor, chunked
print('writing linear predictor')
h5createGroup(h5.path, "linear_predictor")
h5createDataset(h5.path, "linear_predictor/discrete", dim(lp.discrete.mtx),
    storage.mode = "double", chunk=c(1000, 10000), level=5)

h5write(lp.discrete.mtx, h5.path, "linear_predictor/discrete")

h5createDataset(h5.path, "linear_predictor/continuous", dim(lp.continuous.mtx),
    storage.mode = "double", chunk=c(1000, 10000), level=5)

h5write(lp.continuous.mtx, h5.path, "linear_predictor/continuous")

# write mu, chunked
print('writing mu')
h5createGroup(h5.path, "mu")
h5createDataset(h5.path, "mu/discrete", dim(mu.discrete.mtx),
    storage.mode = "double", chunk=c(1000, 10000), level=5)

h5write(mu.discrete.mtx, h5.path, "mu/discrete")

h5createDataset(h5.path, "mu/continuous", dim(mu.continuous.mtx),
    storage.mode = "double", chunk=c(1000, 10000), level=5)

h5write(mu.continuous.mtx, h5.path, "mu/continuous")

# write guide info
print('writing guide info')
h5createGroup(h5.path, "guides")
h5createDataset(h5.path, "guides/one_hot", dim(onehot.guides),
    storage.mode = "integer", chunk=c(1000, 1), level=5)

h5write(onehot.guides, h5.path,"guides/one_hot")
h5write(guides.metadata, h5.path, "guides/metadata")

# write estimate guide efficiencies
print('writing estimated guide efficiencies')
if (!is.null(args$guide_disp)) {
    for (i in 1:length(est.efficiencies.list)) {
        cat(sprintf("writing est efficiencies with dispersion D = %d\n",i))
        h5write(est.efficiencies.list[[i]], h5.path, sprintf("guides/est_efficiency_D%d", args$guide_disp[i]))
    }
}

# write coeffs
print('writing coeffs')
h5write(coeffs, h5.path, "coeffs")

# write x1
print('writing x1')
h5createGroup(h5.path, "x")

h5createDataset(h5.path, "x/x1_discrete", dim(x1.discrete.mtx),
    storage.mode = "integer", chunk = c(1000,1000), level = 7)
h5write(x1.discrete.mtx, h5.path, "x/x1_discrete")

h5createDataset(h5.path, "x/x1_continuous", dim(x1.continuous.mtx),
    storage.mode = "double", chunk = c(1000, 1000), level = 7)
h5write(x1.continuous.mtx, h5.path, "x/x1_continuous")

# write cell cycle scores
print('writing cell cycle scores')
h5write(cell.cycle.scores, h5.path, "x/cell_cycle_scores")

# write percent.mito
print('writing percent.mito')
h5write(percent.mito, h5.path, "x/percent_mito")

# write scaling factors
print('writing scaling factors')
h5write(scaling.factors, h5.path, "scaling_factors")
