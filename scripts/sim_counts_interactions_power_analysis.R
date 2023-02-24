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
library(dplyr)
library(tidyr)
library(patchwork)

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
parser$add_argument("--lambda", action = "store", type = "double", nargs = "+",
	                default = 15,
	                help = "lambda parameter(s) for Poisson distribution to draw number of guides per cell")
parser$add_argument("--effect", action = "store", type = "double", nargs = "+",
                    default = c(0.5, 1, 3, 5, 7),
                    help = "effect size of interaction terms (fixed value)")
parser$add_argument("--neg_effect", action = "store_true", 
                    help = "if true, use negative values for effect sizes")
parser$add_argument("--guide_disp", action = "store", type = "integer", nargs = "+",
                    default = NULL,
                    help = "dispersion value(s) to use when simulating estimated guide efficiencies")
parser$add_argument("--npairs", action = "store", type = "integer",
                    default = 500,
                    help = "number of ground truth pairs with interaction effects")
parser$add_argument("--nneg", action = "store", type = "integer",
                    default = NULL,
                    help = "number of negative controls - pairs targeting same gene w/o interaction")
parser$add_argument("--out", action = "store", type = "character",
	                help = "where to save outputs")
parser$add_argument("--png", action = "store_true", 
                    help = "if true, save PNGs of data metrics")
parser$add_argument("--tiff", action = "store_true",
                    help = "if true, save TIFFs of data metrics")
parser$add_argument("--pdf", action = "store_true", 
                    help = "if true, save PDF of plots")
# parser$add_argument("--data", action = "store", type = "character",
#                     help = "file with empirical data info")
args <- parser$parse_args()

##############################################################
#  write parameters of simulation to file for recordkeeping
##################################å############################
print('writing parameters of simulation job to file')
print(args)

##############################################################
#  simulated baseline (beta0)
#  one beta per gene
#  rnorm params obtained by fitting to experimental data
##############################################################
baselines <- rnorm(args$genes, mean = 2.24, sd = 1.8)

baselines.p <- hist(baselines, xlab = expression(beta[0]), 
                    main = expression(paste("N(", mu, "=2.24, ", sigma,
                    "=1.8)")))

if (args$png) {
    png(file.path(args$out, "baselines.png"))
    print(baselines.p)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "baselines.tiff"))
    print(baselines.p)
    dev.off()
}

if (args$pdf) {
    pdf(file.path(args$out, "baselines.pdf"))
    print(baselines.p)
    dev.off()
}

#######################################
#  assign guides to cells 
#  store as one hot encoding
#######################################
num.guides = args$d * args$targets

cat(sprintf("%d unique gRNAs with %d gRNAs per target (%d total targets)\n",
    num.guides, args$d, args$targets))

# get nr of gRNAs per cell - different set for each value of lambda
guides.per.cell.list <- lapply(args$lambda, function(x) {rpois(args$cells, x)})

### visualize 
guides.per.cell.plotdf <- data.frame(guides.per.cell.list)
colnames(guides.per.cell.plotdf) <- args$lambda
guides.per.cell.plotdf$cell <- 1:args$cells

guides.per.cell.plotdf <- guides.per.cell.plotdf %>% 
                            pivot_longer(-cell, names_to = "lambda", values_to = "nguides")
guides.per.cell.plotdf$lambda <- factor(guides.per.cell.plotdf$lambda, levels = sort(args$lambda))

guides.per.cell.hist <- ggplot(guides.per.cell.plotdf, aes(x = nguides, color = lambda, fill = lambda)) + 
                            geom_histogram(alpha = 0.5, position = "identity") + 
                            theme_classic() +
                            theme(text = element_text(size = 14)) +
                            ggtitle("gRNAs per cell")

if (args$png) {
    png(file.path(args$out, "guides_per_cell.png"))
    print(guides.per.cell.hist)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "guides_per_cell.tiff"))
    print(guides.per.cell.hist)
    dev.off()
}

if (args$pdf) {
    pdf(file.path(args$out, "guides_per_cell.pdf"))
    print(guides.per.cell.hist)
    dev.off()
} 

### assign gRNAs to cells
# initialize list of one hot encodings
onehot.matrices.list <- list()

for (i in 1:length(args$lambda)) {
    onehot.guides <- matrix(0, args$cells, num.guides)

    l <- args$lambda[i]
    cat(sprintf("lambda = %d\n", l))
    
    gpc <- guides.per.cell.list[[i]] 
    
    # assign guides to each cell
    guides.idx.list <- sapply(gpc, function(x) {
        sample(1:num.guides, x, replace = FALSE)
        })
    
    # update one hot encoding
    for (j in 1:args$cells) {
        guides.idx <- guides.idx.list[[j]]
        onehot.guides[j, guides.idx] <- 1
    }
    
    onehot.matrices.list[[i]] <- onehot.guides
    
}

#######################################
#  assign guide efficiencies
#######################################
efficiencies <- rbeta(num.guides, 6,3)
efficiencies.p <- hist(efficiencies,
                    main = expression(paste(beta, "(", a, "=6, ", b,
                    "=3)")))
# visualize
if (args$png) {
    png(file.path(args$out, "efficiencies_hist.png"))
    print(efficiencies.p)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "efficiencies_hist.tiff"))
    print(efficiencies.p)
    dev.off()
}

if (args$pdf) {
    pdf(file.path(args$out, "efficiencies_hist.pdf"))
    print(efficiencies.p)
    dev.off()
}

#######################################
#  assign guides to target sites
#######################################
guide.target.map <- data.frame(guides = 1:num.guides, 
                                target = rep(1:args$targets, args$d),
                                efficiencies = efficiencies)
write.csv(guide.target.map, 
    file.path(args$out, "guide_target_map.csv"), 
    quote = FALSE, row.names = FALSE)

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

### Plot histogram
noisy.df.hist <- noisy.df %>% pivot_longer(!guide, names_to = "D", values_to = "efficiency")

group.colors <- c(`1` = "#DABFFF", `10` = "#907AD6", `100` ="#4F518C", `true` = "#F0EC57")

est.efficiencies.hist.p <- ggplot(noisy.df.hist,
                                aes(x = efficiency, fill = D, color = D)) + 
                            geom_histogram(position = "dodge", alpha = 0.5) + 
                            theme_classic() + 
                            theme(text = element_text(size = 20)) + 
                            scale_fill_manual(values=group.colors) + 
                            scale_color_manual(values = group.colors)


# plot
if (args$png) {
    png(file.path(args$out, "hist_est_guide_efficiencies_all.png"))
    print(est.efficiencies.hist.p)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "hist_est_guide_efficiencies_all.tiff"))
    print(est.efficiencies.hist.p)
    dev.off()
}

if (args$pdf) {
    pdf(file.path(args$out, "hist_est_guide_efficiencies_all.pdf"))
    print(est.efficiencies.hist.p)
    dev.off()
}


### Plot scatterplot
noisy.df.scatter <- noisy.df %>% 
                    pivot_longer(cols = as.character(dispersions), 
                                    names_to = "D", values_to = "noisy")


est.efficiencies.scatter.p <- ggplot(noisy.df.scatter, aes(x = true, y = noisy, col = D)) + 
                                        geom_point() +                         
                                        theme_classic() + 
                                        theme(text = element_text(size = 20)) + 
                                        scale_fill_manual(values=group.colors[1:n.disps]) + 
                                        scale_color_manual(values = group.colors[1:n.disps])


# plot
if (args$png) {
    png(file.path(args$out, "scatter_est_guide_efficiencies_all.png"))
    print(est.efficiencies.scatter.p)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "scatter_est_guide_efficiencies_all.tiff"))
    print(est.efficiencies.scatter.p)
    dev.off()

}

if (args$pdf) {
    pdf(file.path(args$out, "scatter_est_guide_efficiencies_all.pdf"))
    print(est.efficiencies.scatter.p)
    dev.off()

}

### write to file
write.csv(noisy.df, file.path(args$out, "noisy_guide_efficiencies.csv"),
    row.names = TRUE, col.names = TRUE, quote = FALSE)

####################################################
#  determine target sites of interest
#  target sites = putative enhancers
####################################################
# determine total number of target sites
cat(sprintf("number of target sites = %d\n", args$targets))

# total possible target pairs
targets <- 1:args$targets
possible.target.pairs <- combn(targets,2)
cat(sprintf("total possible target pairs = %d\n", dim(possible.target.pairs)[2]))

###########################################################
#  identify target site pairs WITH INTERACTIONS (POS)
#  aka "ground truth" target site pairs
###########################################################
# get ix of pairs from possible pairs
ts.pairs.ix <- sample(1:dim(possible.target.pairs)[2], args$npairs, replace = FALSE)

# subset for ix of pairs and convert to df
ts.pairs <- as.data.frame(t(possible.target.pairs[,ts.pairs.ix]))

# convert to df of target site pairs with interactions 
colnames(ts.pairs) <- c("tsA","tsB")
row.names(ts.pairs) <- c(1:args$npairs)
head(ts.pairs)

#########################################################################
#  identify "negative control" pairs (NEG)
#  aka pairs that act on the same gene without interaction term
#########################################################################
# set number of negative control pairs
n.neg <- args$nneg
if (is.null(args$nneg)) {
    n.neg <- args$npairs
}

# select neg ctrl cases
possible.neg.ix <- setdiff(1:dim(possible.target.pairs)[2], ts.pairs.ix)
neg.pairs.ix <- sample(possible.neg.ix, n.neg, replace = FALSE)

# convert to df 
neg.ctrl.pairs <- as.data.frame(t(possible.target.pairs[,neg.pairs.ix]))
colnames(neg.ctrl.pairs) <- c("neg.tsA", "neg.tsB")
head(neg.ctrl.pairs)

####################################################
#  assign target genes to target site pairs
####################################################

# determine target genes for POS pairs (with interaction)
target.genes <- sample(1:args$genes, args$npairs, replace = FALSE)

# determine target genes for NEG pairs (no interaction)
target.genes.neg <- sample(setdiff(1:args$genes, target.genes), n.neg, replace = FALSE)

# add data to data frames of pos/neg pairs
ts.pairs$target.genes <- target.genes
neg.ctrl.pairs$target.genes <- target.genes.neg

#####################################################################
#  check num. cells where each target site pair is perturbed
#####################################################################
ncells.per.pos.pair.list <- list()

for (i in 1:length(args$lambda)) {
    
    one.hot.guides <- onehot.matrices.list[[i]]
    
    ncells <- c()

    for (j in 1:nrow(ts.pairs)) {
        tsA <- ts.pairs[j,]$tsA
        tsB <- ts.pairs[j,]$tsB
        tsA.guides <- guide.target.map %>% filter(target == tsA) %>% pull(guides)
        tsB.guides <- guide.target.map %>% filter(target == tsB) %>% pull(guides)
        cells.with.tsA <- rowSums(one.hot.guides[,tsA.guides])>0
        cells.with.tsB <- rowSums(one.hot.guides[,tsB.guides])>0
        n <- sum(cells.with.tsA & cells.with.tsB)
        ncells[j] <- n
    }
    
    ncells.per.pos.pair.list[[i]] <- ncells
}

### visualize
ncells.per.pair.df <- data.frame(ncells.per.pos.pair.list)
colnames(ncells.per.pair.df) <- args$lambda
ncells.per.pair.df$pair <- 1:nrow(ts.pairs)

# write to table
write.csv(ncells.per.pair.df, file.path(args$out, "ncells_per_pair.csv"),
            quote = FALSE, row.names = FALSE)

# coerce df for plotting
ncells.per.pair.plotdf <- ncells.per.pair.df %>% 
                            pivot_longer(!pair, names_to = "lambda", values_to = "ncells")

ncells.per.pair.plotdf$lambda <- factor(ncells.per.pair.plotdf$lambda, levels = sort(args$lambda))

# histogram
ts.pair.freq.p <- ggplot(ncells.per.pair.plotdf, aes(x = ncells, color = lambda, fill = lambda)) + 
                    geom_histogram(alpha = 0.5, position = "identity") + 
                    theme_classic() +
                    theme(text = element_text(size = 14)) +
                    ggtitle("Frequency of target site pairs") 

# plot
if (args$png) {
    png(file.path(args$out, "ts_pair_freq_hist.png"))
    print(ts.pair.freq.p)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "ts_pair_freq_hist.tiff"))
    print(ts.pair.freq.p)
    dev.off()

}

if (args$pdf) {
    pdf(file.path(args$out, "ts_pair_freq_hist.pdf"))
    print(ts.pair.freq.p)
    dev.off()

}

# boxplot
ts.pair.freq.box <- ggplot(ncells.per.pair.plotdf, aes(x = lambda, y = ncells, color = lambda)) + 
                    geom_boxplot() + 
                    theme_classic() +
                    theme(text = element_text(size = 14)) +
                    ggtitle("Frequency of target site pairs")

# plot
if (args$png) {
    png(file.path(args$out, "ts_pair_freq_boxplot.png"))
    print(ts.pair.freq.box)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "ts_pair_freq_boxplot.tiff"))
    print(ts.pair.freq.box)
    dev.off()

}

if (args$pdf) {
    pdf(file.path(args$out, "ts_pair_freq_boxplot.pdf"))
    print(ts.pair.freq.box)
    dev.off()

}                    

# plot histogram and boxplot together
if (args$png) {
    png(file.path(args$out, "ts_pair_freq_patchwork.png"))
    print(ts.pair.freq.p  | ts.pair.freq.box)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "ts_pair_freq_patchwork.tiff"))
    print(ts.pair.freq.p  | ts.pair.freq.box)
    dev.off()

}

if (args$pdf) {
    pdf(file.path(args$out, "ts_pair_freq_patchwork.pdf"))
    print(ts.pair.freq.p  | ts.pair.freq.box)
    dev.off()
}


#####################################################################
#  assign effect sizes for individual target sites on target gene
#  betaA, betaB, betaAB
#####################################################################
ts.pairs$betaA <- -1*rgamma(dim(ts.pairs)[1], shape = 6, scale = 0.5)
ts.pairs$betaB <- -1*rgamma(dim(ts.pairs)[1], shape = 6, scale = 0.5)

neg.ctrl.pairs$betaA <- -1*rgamma(dim(neg.ctrl.pairs)[1], shape = 6, scale = 0.5)
neg.ctrl.pairs$betaB <- -1*rgamma(dim(neg.ctrl.pairs)[1], shape = 6, scale = 0.5)

# combine pos + neg pairs into one table
neg.ctrl.pairs.tmp <- neg.ctrl.pairs
colnames(neg.ctrl.pairs.tmp)[1:2] <- c("tsA","tsB")
all.ts.pairs <- rbind(ts.pairs, neg.ctrl.pairs.tmp)
all.ts.pairs$set <- c(rep("positive", args$npairs), rep("negative", n.neg))

####################################################
#  simulate X_S, X_G2M (cell cycle scores)
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

# collect data in data frame
cell.cycle.scores <- data.frame(cell = 1:args$cells, s = s.scores, g2m = g2m.scores)
cell.cycle.scores.plotdf <- cell.cycle.scores %>% 
                            tidyr::pivot_longer(!cell, names_to = "cell cycle", values_to = "score")


if (args$png) {
    png(file.path(args$out, "cell_cycle_scores_hist.png"))
    print(ggplot(cell.cycle.scores.plotdf, 
        aes(x = score, fill = `cell cycle`, color = `cell cycle`)) + 
            geom_histogram(position = "dodge", alpha = 0.5) + 
            theme_classic() + theme(text = element_text(size = 20)) 
            )
    dev.off()
}

if (args$pdf) {
    pdf(file.path(args$out, "cell_cycle_scores_hist.pdf"))
    print(ggplot(cell.cycle.scores.plotdf, 
        aes(x = score, fill = `cell cycle`, color = `cell cycle`)) + 
            geom_histogram(position = "dodge", alpha = 0.5) + 
            theme_classic() + theme(text = element_text(size = 20)) 
            )
    dev.off()
}

# write cell cycle scores (X2, X3) to file
# row index = cell identifier
cell.cycle.scores <- data.frame(s.scores, g2m.scores)
write.table(cell.cycle.scores, file.path(args$out, "cell_cycle_scores.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate beta_S, beta_G2M (from gamma distr.)
####################################################
beta.s <- rgamma(args$genes, shape = 6, scale = 0.5)
beta.g2m <- rgamma(args$genes, shape = 6, scale = 0.5)

if (args$png) {
    png(file.path(args$out, "beta_S_hist.png"))
    hist(beta.s)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "beta_S_hist.tiff"))
    hist(beta.s)
    dev.off()
}

if (args$png) {
    png(file.path(args$out, "beta_G2M_hist.png"))
    hist(beta.g2m)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "beta_G2M_hist.tiff"))
    hist(beta.g2m)
    dev.off()
}

####################################################
#  simulate beta_pct.mito, x_pct.mito
####################################################
percent.mito <- rbeta(args$cells, shape1 = 3.3, shape2 = 81.48)
beta.pct.mito <- rgamma(args$genes, shape = 6 , scale = 0.5)

if (args$png) {
    png(file.path(args$out, "percent_mito_hist.png"))
    hist(percent.mito)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "percent_mito_hist.tiff"))
    hist(percent.mito)
    dev.off()
}

if (args$png) {
    png(file.path(args$out, "beta_percent_mito_hist.png"))
    hist(beta.pct.mito)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "beta_percent_mito_hist.tiff"))
    hist(beta.pct.mito)
    dev.off()
}

# write percent.mito to file
# row index = cell identifier
percent.mito.df <- data.frame(percent.mito)
write.table(percent.mito.df, file.path(args$out, "percent_mito.txt"), row.names = TRUE, quote = FALSE)

######################################################
#  write all "true" coeff values to table
######################################################
all.target.genes <- all.ts.pairs$target.gene
nt.genes <- setdiff(1:args$genes, all.target.genes)

# get table of betaA, betaB, and interaction coeffs for NT genes (0)
ab.coeffs.nt <- data.frame(gene = nt.genes, betaA = 0, betaB = 0)

# subset coeffs for pos and neg pairs and rename columns to match table for NT genes
ab.coeffs.targets <- all.ts.pairs %>% select(target.genes, betaA, betaB)
colnames(ab.coeffs.targets)[1] <- "gene"

# bind tables
ab.coeffs.all <- rbind(ab.coeffs.targets, ab.coeffs.nt) %>% arrange(gene)

# add all other coeffs to table
all.coeffs <- ab.coeffs.all
all.coeffs$beta0 <- baselines
all.coeffs$beta.s <- beta.s
all.coeffs$beta.g2m <- beta.g2m
all.coeffs$percent.mito <- beta.pct.mito

print('all coeffs:')
head(all.coeffs)

# row index = gene identifier
write.table(all.coeffs, file.path(args$out, "coeffs.txt"), row.names = TRUE, quote = FALSE)

####################################################
#  simulate scaling factors (poisson)
####################################################

# simulate total counts per cell
t.vec <- rpois(args$cells, 50000)

# get scaling factors
scaling.factors <- t.vec / 1e6

# plot distribution of scaling factors
if (args$png) {
    png(file.path(args$out, "scaling_factors_hist.png"))
    hist(scaling.factors)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "scaling_factors_hist.tiff"))
    hist(scaling.factors)
    dev.off()
}

if (args$pdf) {
    pdf(file.path(args$out, "scaling_factors_hist.pdf"))
    hist(scaling.factors)
    dev.off()
}
# write scaling factors to file
# row index = cell identifier
write.table(data.frame(scaling.factors), file.path(args$out, "scaling_factors.txt"), 
                       row.names = TRUE, quote = FALSE)

####################################################
#  initialize H5 for writing simulated data
####################################################
print("writing data to h5")
h5.path <- file.path(args$out, "sim.h5")
h5createFile(h5.path)

# create groups
print('creating groups')
h5createGroup(h5.path, "counts")
h5createGroup(h5.path, "linear_predictor")
h5createGroup(h5.path, "mu")
h5createGroup(h5.path, "x")
h5createGroup(h5.path, "x/x_a")
h5createGroup(h5.path, "x/x_b")
h5createGroup(h5.path, "x/x_ab")
h5createGroup(h5.path, "guides")
h5createGroup(h5.path, "guides/one_hot")

####################################################
#  simulate counts matrix 
####################################################

# initialize X_A, X_B, X_AB as zero 
# one matrix per value of lambda
print('initializing matrices for XA, XB, XAB')
xa.mtx.list <- replicate(length(args$lambda),
                        matrix(0,args$genes, args$cells), simplify = FALSE)
xb.mtx.list <- replicate(length(args$lambda),
                        matrix(0,args$genes, args$cells), simplify = FALSE)
xab.mtx.list <- replicate(length(args$lambda),
                        matrix(0,args$genes, args$cells), simplify = FALSE)


# initialize counts matrices list
# one matrix per value of lambda and interaction size
# total = length(lambda)*length(effect sizes)
# sim.counts.list <- replicate(length(args$lambda)*length(args$effect),
#                             matrix(0,args$genes, args$cells), simplify = FALSE)
# sim.counts.list <- lapply(1:length(args$lambda), function(x) {
#         lapply(1:length(args$effect), function(x) {
#                 matrix(0, args$genes, args$cells)
#             })
#     })
print('initializing simulated counts matrices')
sim.counts.list <- replicate(length(args$lambda), 
                            replicate(length(args$effect), 
                                matrix(0,args$genes, args$cells), 
                                simplify = FALSE), 
                            simplify = FALSE)

# initialize matrix lists for storing LP/mu
# lp.mtx.list <- replicate(length(args$lambda)*length(args$effect),
#                          matrix(0,args$genes, args$cells), simplify = FALSE)
print('initializing LP matrices')
lp.mtx.list <- replicate(length(args$lambda), 
                            replicate(length(args$effect), 
                                matrix(0,args$genes, args$cells), 
                                simplify = FALSE), 
                            simplify = FALSE)
# mu.mtx.list <- replicate(length(args$lambda)*length(args$effect),
#                          matrix(0,args$genes, args$cells), simplify = FALSE)
print('initializing mu matrices')
mu.mtx.list <- replicate(length(args$lambda), 
                            replicate(length(args$effect), 
                                matrix(0,args$genes, args$cells), 
                                simplify = FALSE), 
                            simplify = FALSE)


### test with genes targeted by pair with interaction term
pos.test.genes <- sample(ts.pairs$target.genes, 10, replace = FALSE)

# iterate through genes
for (gene in 1:args$genes) {
# for (gene in pos.test.genes) {
    cat(sprintf("simulating counts for gene %d\n", gene))
    
    # get coeffs 
    b0 <- baselines[gene]
    beta.a <- all.coeffs$betaA[gene]
    beta.b <- all.coeffs$betaB[gene]
    beta.ab <- 0
    beta.s.score <- all.coeffs$beta.s[gene]
    beta.g2m.score <- all.coeffs$beta.g2m[gene]
    beta.mito <- all.coeffs$percent.mito[gene]
    
    # get cell cycle scores + percent.mito 
    # these variables are fixed for all values of lambda
    x.s <- s.scores
    x.g2m <- g2m.scores
    x.mito <- percent.mito
    
    # check if gene is targeted by any gRNAs (is it in pos/neg pair)
    if (gene %in% all.ts.pairs$target.genes) {
        cat(sprintf("gene %d is targeted by gRNAs in this design\n", gene))
        
        # identify tsA and tsB
        tmp <- all.ts.pairs %>% filter(target.genes == gene) 
        tsA <- tmp$tsA
        tsB <- tmp$tsB
        
        # get guides targeting tsA/tsB
        guides.A <- guide.target.map %>% filter(target == tsA) %>% pull(guides)
        guides.B <- guide.target.map %>% filter(target == tsB) %>% pull(guides)
        
        # calculate X_A and X_B for different values of lambda
        for (i in 1:length(args$lambda)) {
            l <- args$lambda[i]
            cat(sprintf("calculating X_A and X_B for lambda = %d\n", l))
            
            # calculate X_A (probability of tsA perturbation)
            temp.mtx.a <- t(efficiencies[guides.A]*t(onehot.matrices.list[[i]][,guides.A]))
            xa <- apply(temp.mtx.a, 1, function(x) {1-prod(1-x)})

            # calculate X_B (probability of tsB perturbation)
            temp.mtx.b <- t(efficiencies[guides.B]*t(onehot.matrices.list[[i]][,guides.B]))
            xb <- apply(temp.mtx.b, 1, function(x) {1-prod(1-x)})
            
            xa.mtx.list[[i]][gene,] <- xa
            xb.mtx.list[[i]][gene,] <- xb

            xab <- xa*xb
            xab.mtx.list[[i]][gene,] <- xab
            
            # check if gene is targeted by an interaction effect (pos pair)
            if (gene %in% ts.pairs$target.genes) {
                print('gene has interaction effect')
                for (j in 1:length(args$effect)) {
                    size <- args$effect[j]
                    if (args$neg_effect) {
                        size <- -1*size
                    }
                    cat(sprintf("interaction size = %.2f\n", size))
                    beta.ab <- size
                    # calculate lp and mu 
                    lp <- b0 + beta.a*xa + beta.b*xb + beta.ab*xab + 
                                beta.s.score*x.s + beta.g2m.score*x.g2m + beta.mito*x.mito + 
                                log(scaling.factors)
                    # lp.mtx.list[[i*j]][gene,] <- lp
                    lp.mtx.list[[i]][[j]][gene,] <- lp
                    mu <- exp(lp)
                    # mu.mtx.list[[i*j]][gene,] <- mu
                    mu.mtx.list[[i]][[j]][gene,] <- mu

                    # use rnbinom to generate counts for gene for each cell and update counts matrix 
                    counts <- rnbinom(length(mu), mu = mu, size = 1.5)
                    # sim.counts.list[[i*j]][gene,] <- counts
                    sim.counts.list[[i]][[j]][gene,] <- counts
                }
            } else { 
                # gene is not targeted by interaction effect (neg ctrl pair)
                print('gene has no interaction effect')
                # calculate lp and mu 
                lp <- b0 + beta.a*xa + beta.b*xb +
                            beta.s.score*x.s + 
                            beta.g2m.score*x.g2m + 
                            beta.mito*x.mito + 
                            log(scaling.factors)
                mu <- exp(lp)
                counts <- rnbinom(length(mu), mu = mu, size = 1.5)
                
                for (j in 1:length(args$effect)) {
                    # lp.mtx.list[[i*j]][gene,] <- lp
                    # mu.mtx.list[[i*j]][gene,] <- mu
                    # sim.counts.list[[i*j]][gene,] <- counts

                    lp.mtx.list[[i]][[j]][gene,] <- lp
                    mu.mtx.list[[i]][[j]][gene,] <- mu
                    sim.counts.list[[i]][[j]][gene,] <- counts
                } 
            }
        }
    } else {
        # gene is not targeted by any guides
        print('gene is not targeted by any gRNAs in this design')
        lp <- b0 + beta.s.score*x.s + beta.g2m.score*x.g2m + beta.mito*x.mito + 
                                log(scaling.factors)
        mu <- exp(lp)
        counts <- rnbinom(length(mu), mu = mu, size = 1.5)
        for (i in 1:length(args$lambda)) {
            for (j in 1:length(args$effect)) {
                # lp.mtx.list[[i*j]][gene,] <- lp
                # mu.mtx.list[[i*j]][gene,] <- mu
                # sim.counts.list[[i*j]][gene,] <- counts
                lp.mtx.list[[i]][[j]][gene,] <- lp
                mu.mtx.list[[i]][[j]][gene,] <- mu
                sim.counts.list[[i]][[j]][gene,] <- counts
            }
        }
    }
}

####################################################
#  write data to h5
####################################################
# write everything to h5
for (i in 1:length(args$lambda)) {
    l <- args$lambda[i]
    
    # write x values pertaining to perturbations
    print('writing X_A, X_B, X_AB')

    h5createDataset(h5.path, paste0("x/x_a/", l), dim(xa.mtx.list[[i]]),
        storage.mode = "double", chunk = c(1000,1000), level = 7)
    h5write(xa.mtx.list[[i]], h5.path, paste0("x/x_a/", l))

    h5createDataset(h5.path, paste0("x/x_b/", l), dim(xb.mtx.list[[i]]),
        storage.mode = "double", chunk = c(1000,1000), level = 7)
    h5write(xb.mtx.list[[i]], h5.path, paste0("x/x_b/", l))

    h5createDataset(h5.path, paste0("x/x_ab/", l), dim(xab.mtx.list[[i]]),
        storage.mode = "double", chunk = c(1000,1000), level = 7)
    h5write(xab.mtx.list[[i]], h5.path, paste0("x/x_ab/", l))

    # write one-hot encoded matrices
    print('writing one-hot encoded matrices for gRNAs in cells')
    h5createDataset(h5.path, paste0("guides/one_hot/",l), dim(onehot.matrices.list[[i]]),
    storage.mode = "integer", chunk=c(1000, 1), level=5)
    h5write(onehot.matrices.list[[i]], h5.path, paste0("guides/one_hot/",l))

    for (j in 1:length(args$effect)) {
        size <- args$effect[j]
        set.name <- paste0("lambda",l, "_size", size)
        if (args$neg_effect) {
            set.name <- paste0("lambda",l,"_size-NEG", size)
        }
        print(set.name)
        
        # write counts and relevant info for this set of params to h5
        print(paste('writing to h5 simulations for', set.name))
        print('writing counts, chunked')
        h5createDataset(h5.path, paste0("counts/", set.name), 
                        dim(sim.counts.list[[i]][[j]]), 
                        storage.mode = "integer", chunk = c(1000,10000), level = 5)
        # h5write(sim.counts.list[[i*j]], h5.path, paste0("counts/",set.name))
        h5write(sim.counts.list[[i]][[j]], h5.path, paste0("counts/",set.name))

        
        print('writing LP')
        h5createDataset(h5.path, paste0("linear_predictor/", set.name), 
                        dim(lp.mtx.list[[i]][[j]]),
                        storage.mode = "double", chunk=c(1000, 10000), level=5)
        # h5write(lp.mtx.list[[i*j]], h5.path, paste0("linear_predictor/", set.name))
        h5write(lp.mtx.list[[i]][[j]], h5.path, paste0("linear_predictor/", set.name))
        
        print('writing mu')
        h5createDataset(h5.path, paste0("mu/", set.name), dim(mu.mtx.list[[i]][[j]]),
                        storage.mode = "double", chunk=c(1000, 10000), level=5)
        # h5write(mu.mtx.list[[i*j]], h5.path, paste0("mu/", set.name))
        h5write(mu.mtx.list[[i]][[j]], h5.path, paste0("mu/", set.name))
        

    }
}

# write table with info about pairs targeting genes (POS & NEG)
print('writing info about target site pairs')
h5write(all.ts.pairs, h5.path, "pairs")

# write guide info
print('writing guide info')
h5write(guide.target.map, h5.path, "guides/guide_target_map")

if (!is.null(args$guide_disp)) {
    print('writing estimated guide efficiencies')
    h5write(noisy.df, h5.path, "guides/noisy_guide_efficiencies")
}

# write coeffs
print('writing coeffs')
h5write(all.coeffs, h5.path, "coeffs")

# write cell cycle scores
print('writing cell cycle scores')
h5write(cell.cycle.scores, h5.path, "x/cell_cycle_scores")

# write percent.mito
print('writing percent.mito')
h5write(percent.mito, h5.path, "x/percent_mito")

# write scaling factors
print('writing scaling factors')
h5write(scaling.factors, h5.path, "scaling_factors")

# writing lambda and effect size values used 
print('writing lambda and effect sizes')
h5write(args$lambda, h5.path, "lambda")
h5write(args$effect, h5.path, "effect.sizes")
