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
parser$add_argument("--pct", action = "store", type = "double",
                    help = "pick percent of targets to define as true interactions",
                    default = 0.5)
parser$add_argument("--out", action = "store", type = "character",
	                help = "where to save outputs")
parser$add_argument("--png", action = "store_true", 
                    help = "if true, save PNGs of data metrics")
parser$add_argument("--tiff", action = "store_true",
                    help = "if true, save TIFFs of data metrics")
# parser$add_argument("--data", action = "store", type = "character",
#                     help = "file with empirical data info")
args <- parser$parse_args()

##############################################################
#  write parameters of simulation to file for recordkeeping
##################################å############################
print('writing parameters of simulation job to file')
print(args)

args.df <- data.frame(args)
args.df$date <- Sys.Date()

print(args.df)

print(data.frame(t(args.df)))

write.csv(data.frame(t(args.df)), file.path(args$out, "simulation_params.csv"),
    row.names = TRUE, col.names = FALSE, quote = FALSE)

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
num.guides = args$d * args$targets

cat(sprintf("%d unique gRNAs with %d gRNAs per target (%d total targets)\n",
    num.guides, args$d, args$targets))

# initialize one hot encoding  
onehot.guides <- matrix(0, args$cells, num.guides)

# get nr of gRNAs per cell
guides.per.cell <- rpois(args$cells, 15)

# visualize 
if (args$png) {
    png(file.path(args$out, "guides_per_cell.png"))
    hist(guides.per.cell)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "guides_per_cell.tiff"))
    hist(guides.per.cell)
    dev.off()
}

# assign guides to cells
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

if (args$png) {
    png(file.path(args$out, "efficiencies_hist.png"))
    print(hist(efficiencies,
                main = expression(paste(beta, "(", a, "=6, ", b,
                "=3)"))))
    dev.off()
}

if (args$tiff) {
    png(file.path(args$out, "efficiencies_hist.tiff"))
    print(hist(efficiencies,
                main = expression(paste(beta, "(", a, "=6, ", b,
                "=3)"))))
    dev.off()
}

#######################################
#  assign guides to target sites
#######################################
guide.target.map <- data.frame(guides = 1:num.guides, target = rep(1:args$targets, args$d))
write.csv(guide.target.map, file.path(args$out, "guide_target_map.csv"), quote = FALSE, row.names = FALSE)

#######################################
#  make guides metadata table
#######################################
guides.metadata <- guide.target.map
guides.metadata$efficiencies <- efficiencies
write.csv(guides.metadata, file.path(args$out, "guides_metadata.csv"), quote = FALSE, row.names = FALSE)

####################################################
#  assign guide efficiencies to guides in library
####################################################
efficiencies <- rbeta(num.guides, 6,3)

# visualize
if (args$png) {
    png(file.path(args$out, "guide_efficiencies_hist.png"))
    hist(efficiencies)
    dev.off()
}

if (args$tiff) {
    tiff(file.path(args$out, "guide_efficiencies_hist.tiff"))
    hist(efficiencies)
    dev.off()
}

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


### Plot scatterplot
noisy.df.scatter <- noisy.df %>% pivot_longer(cols = as.character(dispersions), names_to = "D", values_to = "noisy")


est.efficiencies.scatter.p <- ggplot(noisy.df.scatter, aes(x = true, y = noisy, col = D)) + geom_point() +                         
                                        theme_classic() + 
                                        theme(text = element_text(size = 20)) + 
                                        scale_fill_manual(values=group.colors[1:3]) + 
                                        scale_color_manual(values = group.colors[1:3])


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

### write to file
write.table(noisy.df, file.path(args$out, "noisy_guide_efficiencies.csv"),
    row.names = TRUE, col.names = TRUE, quote = FALSE)


####################################################
#  determine target sites of interest
#  target sites = putative enhancers
####################################################
# determine total number of targets 
targets <- 1:args$targets
cat(sprintf("number of target sites = %d\n", args$targets))

# total possible target pairs
possible.target.pairs <- combn(targets,2)

cat(sprintf("total possible target pairs = %d\n", dim(possible.target.pairs)[2]))

# # determine percent of possible pairs to define as true interactions
# pct.true.interactions <- args$pct

# # see number of pairs for each percentage
# n.pairs <- floor(ceiling(pct.true.interactions*dim(possible.target.pairs)[2]), args$targets)

# set number of pairs with interactions based on percentage argument
n.pairs <- ceiling(args$pct*args$targets)

# # make df and visualize
# interactions.df <- data.frame(pct.pairs = pct.true.interactions, n.pairs = n.pairs)

###########################################################
#  identify target site pairs WITH INTERACTIONS (POS)
###########################################################
# get ix of pairs from possible pairs
ts.pairs.ix <- sample(1:dim(possible.target.pairs)[2], n.pairs, replace = FALSE)

# subset for ix of pairs and convert to df
ts.pairs <- as.data.frame(t(possible.target.pairs[,ts.pairs.ix]))

# convert to df of target site pairs with interactions 
colnames(ts.pairs) <- c("tsA","tsB")
row.names(ts.pairs) <- c(1:n.pairs)

#########################################################################
#  identify "negative control" pairs (NEG)
#  defined as pairs that act on the same gene without interaction term
#########################################################################
# set number of negative control pairs
n.neg <- args$targets - n.pairs

# select neg ctrl cases
possible.neg.ix <- setdiff(1:dim(possible.target.pairs)[2], ts.pairs.ix)
neg.pairs.ix <- sample(possible.neg.ix, n.neg, replace = FALSE)

# convert to df 
neg.ctrl.pairs <- as.data.frame(t(possible.target.pairs[,neg.pairs.ix]))
colnames(neg.ctrl.pairs) <- c("neg.tsA", "neg.tsB")

####################################################
#  assign target genes to target site pairs
####################################################

# determine target genes for POS pairs (with interaction)
target.genes <- sample(1:args$genes, n.pairs, replace = FALSE)

# determine target genes for NEG pairs (no interaction)
target.genes.neg <- sample(setdiff(1:args$genes, target.genes), n.neg, replace = FALSE)

# add data to data frames of pos/neg pairs
ts.pairs$target.genes <- target.genes
neg.ctrl.pairs$target.genes <- target.genes.neg

#####################################################################
#  assign effect sizes for individual target sites on target gene
#  betaA, betaB, betaAB
#####################################################################
assign_ts_effect <- function(df, shape = 6, scale = 0.5, interaction = TRUE) {
    n <- dim(df)[1]
    beta.a <- -1 * (rgamma(n, shape = shape, scale = scale))
    beta.b <- -1 * (rgamma(n, shape = shape, scale = scale))
    df$beta.A <- beta.a
    df$beta.B <- beta.b
    if (interaction) {
        interaction.pos <- 1*rgamma(ceiling(n/2), shape = shape, scale = scale)
        interaction.neg <- rgamma(n-ceiling(n/2), shape = shape, scale = scale)
        interactions <- c(interaction.pos, interaction.neg)
        df$interaction <- interactions 
    }
    return(df)
}

# assign target site effect sizes for pairs in POS/NEG pairs
ts.pairs <- assign_ts_effect(ts.pairs)
neg.ctrl.pairs <- assign_ts_effect(neg.ctrl.pairs, interaction = FALSE)
neg.ctrl.pairs$interaction <- 0

# combine tables of POS/NEG target sites pairs
neg.ctrl.pairs.tmp <- neg.ctrl.pairs
colnames(neg.ctrl.pairs.tmp)[1:2] <- c("tsA","tsB")
all.pairs.pos.and.neg <- rbind(ts.pairs, neg.ctrl.pairs.tmp)
all.pairs.pos.and.neg$set <- c(rep("positive", n.pairs), rep("negative", n.neg))

head(all.pairs.pos.and.neg)

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

# if (args$pdf) {
#     png(file.path(args$out, "cell_cycle_scores_hist.pdf"))
#     print(ggplot(cell.cycle.scores.plotdf, 
#         aes(x = score, fill = `cell cycle`, color = `cell cycle`)) + 
#             geom_histogram(position = "dodge", alpha = 0.5) + 
#             theme_classic() + theme(text = element_text(size = 20)) 
#             )
#     dev.off()
# }

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
all.target.genes <- all.pairs.pos.and.neg$target.gene
nt.genes <- setdiff(1:args$genes, all.target.genes)

# get table of betaA, betaB, and interaction coeffs for NT genes (0)
ab.coeffs.nt <- data.frame(gene = nt.genes, beta.A = 0, beta.B = 0, interaction = 0)

# subset coeffs for pos and neg pairs and rename columns to match table for NT genes
ab.coeffs.targets <- all.pairs.pos.and.neg %>% select(target.genes, beta.A, beta.B, interaction)
colnames(ab.coeffs.targets)[1] <- "gene"

# bind tables
ab.coeffs.all <- rbind(ab.coeffs.targets, ab.coeffs.nt) %>% arrange(gene)

# add all other coeffs to table
all.coeffs <- ab.coeffs.all
all.coeffs$beta0 <- baselines
all.coeffs$beta.s <- beta.s
all.coeffs$beta.g2m <- beta.g2m
all.coeffs$percent.mito <- beta.pct.mito

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
# write scaling factors to file
# row index = cell identifier
write.table(data.frame(scaling.factors), file.path(args$out, "scaling_factors.txt"), 
                       row.names = TRUE, quote = FALSE)

####################################################
#  simulate counts matrix 
####################################################

# initialize counts matrix
sim.counts <- matrix(0, args$genes, args$cells)

# initialize matrix of x values pertaining to perturbations 
xa.mtx <- matrix(0, args$genes, args$cells)
xb.mtx <- matrix(0, args$genes, args$cells)
xab.mtx <- matrix(0, args$genes, args$cells)

# initialize matrices for storing values of linear predictor and mu
lp.mtx <- matrix(0, args$genes, args$cells)
mu.mtx <- matrix(0, args$genes, args$cells)

# populate counts mtx by gene (by row)
for (gene in 1:args$genes) {
# for (gene in sample(args$genes, 10)) {
    cat(sprintf("simulating counts for gene %d\n", gene))

    # get coeffs
    b0 <- baselines[gene]
    beta.a <- all.coeffs$beta.A[gene]
    beta.b <- all.coeffs$beta.B[gene]
    beta.ab <- all.coeffs$interaction[gene]
    beta.s.score <- beta.s[gene]
    beta.g2m.score <- beta.g2m[gene]
    beta.mito <- beta.pct.mito[gene]
    
    # initialize X_A, X_B, X_AB as zero 
    xa <- numeric(args$cells)
    xb <- numeric(args$cells)
    xab <- numeric(args$cells)
    
    # check if enhancer of gene is targeted by any gRNAs
    if (gene %in% ab.coeffs.targets$gene) {
        cat(sprintf("gene %d is targeted by gRNAs in this design\n", gene))

        # identify tsA and tsB
        tmp <- all.pairs.pos.and.neg %>% filter(target.genes == gene) 
        tsA <- tmp$tsA
        tsB <- tmp$tsB

        # get guides targeting tsA 
        guides.A <- guide.target.map %>% filter(target == tsA) %>% pull(guides)
        guides.B <- guide.target.map %>% filter(target == tsB) %>% pull(guides)
        
        # calculate X_A (probability of tsA perturbation)
        temp.mtx.a <- t(efficiencies[guides.A]*t(onehot.guides[,guides.A]))
        xa <- apply(temp.mtx.a, 1, function(x) {1-prod(1-x)})
        
        # calculate X_B (probability of tsB perturbation)
        temp.mtx.b <- t(efficiencies[guides.B]*t(onehot.guides[,guides.B]))
        xb <- apply(temp.mtx.b, 1, function(x) {1-prod(1-x)})
        
        # calculate X_AB (probability of both perturbed)
        xab <- xa*xb

        # store values
        xa.mtx[gene,] <- xa
        xb.mtx[gene,] <- xb
        xab.mtx[gene,] <- xab
    }

    # get cell cycle scores + percent.mito
    x.s <- s.scores
    x.g2m <- g2m.scores
    x.pct.mito <- percent.mito
    
    # calculate linear predictor/mu
    lp <- b0 + beta.a*xa + beta.b*xb + beta.ab*xab + 
                beta.s.score*x.s + beta.g2m.score*x.g2m + beta.mito*x.pct.mito + 
                log(scaling.factors)
    lp.mtx[gene,] <- lp
    
    mu <- exp(lp)
    mu.mtx[gene,] <- mu

    # use rnbinom to generate counts for gene for each cell and update counts matrix 
    counts <- rnbinom(length(mu), mu = mu, size = 1.5)
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

# write mu, chunked
print('writing mu')
h5createGroup(h5.path, "mu")
h5createDataset(h5.path, "mu/mu", dim(mu.mtx),
    storage.mode = "double", chunk=c(1000, 10000), level=5)

h5write(mu.mtx, h5.path, "mu/mu")

# write table with info about pairs targeting genes (POS & NEG)
print('writing info about pairs')
h5write(all.pairs.pos.and.neg, h5.path,"pairs")

# write guide info
print('writing guide info')
h5createGroup(h5.path, "guides")
h5createDataset(h5.path, "guides/one_hot", dim(onehot.guides),
    storage.mode = "integer", chunk=c(1000, 1), level=5)

h5write(onehot.guides, h5.path,"guides/one_hot")
h5write(guide.target.map, h5.path, "guides/guide_target_map")
h5write(guides.metadata, h5.path, "guides/metadata")

if (!is.null(args$guide_disp)) {
    print('writing estimated guide efficiencies')
    h5write(noisy.df, h5.path, "guides/noisy_guide_efficiencies")
}

# write coeffs
print('writing coeffs')
h5write(all.coeffs, h5.path, "coeffs")

# write x values pertaining to perturbaitons
print('writing X_A, X_B, X_AB')
h5createGroup(h5.path, "x")

h5createDataset(h5.path, "x/x_a", dim(xa.mtx),
    storage.mode = "double", chunk = c(1000,1000), level = 7)
h5write(xa.mtx, h5.path, "x/x_a")

h5createDataset(h5.path, "x/x_b", dim(xb.mtx),
    storage.mode = "double", chunk = c(1000,1000), level = 7)
h5write(xb.mtx, h5.path, "x/x_b")

h5createDataset(h5.path, "x/x_ab", dim(xab.mtx),
    storage.mode = "double", chunk = c(1000,1000), level = 7)
h5write(xab.mtx, h5.path, "x/x_ab")

# write cell cycle scores
print('writing cell cycle scores')
h5write(cell.cycle.scores, h5.path, "x/cell_cycle_scores")

# write percent.mito
print('writing percent.mito')
h5write(percent.mito, h5.path, "x/percent_mito")

# write scaling factors
print('writing scaling factors')
h5write(scaling.factors, h5.path, "scaling_factors")
