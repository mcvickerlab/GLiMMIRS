### this code runs a power analysis on simulated dataset for interactions between regulatory elements

library(argparse)
# library(Matrix)
# library(ggplot2)
library(MASS)
library(rhdf5)
library(dplyr)
library(broom)

# define parser to handle input arguments from command line
parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--h5', action = "store", type = "character", 
                    help = "path to h5 file with simulated data")
parser$add_argument("--cells", action = "store", type = 'integer',
	                default = 50000,
	                help = "number of cells in simualted dataset")
parser$add_argument('--genes', action = "store", type = "integer", 
	                default = 13000,
                    help = "number of genes in simulated dataset")
parser$add_argument("--out", action = "store", type = "character",
	                help = "where to save outputs")
parser$add_argument("--pseudocount", action = "store", type = "double",
					default = NULL,
					help = "pseudocount to add when modeling")
parser$add_argument("--lambda", action = "store", type = "integer",
					# default = NULL,
					help = "specify which MOI to do power analysis for; if NULL test all")
# parser$add_argument("--pos", action = "store_true",
# 					help = "if true, test positive pairs, e.g. enhancers with interaction")
parser$add_argument("--neg", action = "store_true",
					help = "if true, test negative pairs, e.g. enhancers without interactions (DEFAULT=test positive pairs")
parser$add_argument('--test', action = 'store_true',
					help = "if true, just test five genes")
parser$add_argument("--min.cells", action = "store", type = "integer",
					default = 20,
					help = "minimum number of cells containing perturbations for both enhancers")
args <- parser$parse_args()

print(args)

# show contents of h5 file
h5ls(args$h5)

### define helper functions 
load_counts <- function(lambda, effect.sizes) {
	# lambda should be fixed value
	# effect.sizes is a list of effect sizes
	sim.counts.list <- list()

	for (j in 1:length(effect.sizes)) {
		set.name <- paste0("lambda", args$lambda, "_size", effect.sizes.list[j])
	    counts <- h5read(args$h5, paste0("counts/", set.name))
	    sim.counts.list[[j]] <- counts
	}

	return(sim.counts.list)
}

### load important info
all.ts.pairs <- h5read(args$h5, "pairs")
lambda.list <- h5read(args$h5, "lambda")
effect.sizes.list <- h5read(args$h5, "effect.sizes")
guide.target.map <- h5read(args$h5, "guides/guide_target_map")
coeffs <- h5read(args$h5, "coeffs")
cell.cycle.scores <- h5read(args$h5, "x/cell_cycle_scores")
percent.mito <- h5read(args$h5, "x/percent_mito")
scaling.factors <- h5read(args$h5, "scaling_factors")

# get list of positive/negative control cases of interactions 
genes.to.test <- c()
if (args$neg) {
	genes.to.test <- all.ts.pairs %>% filter(set == "negative") %>% pull(target.genes)
} else {
	genes.to.test <- all.ts.pairs %>% filter(set == "positive") %>% pull(target.genes)
}

# if test flag is set, iterate through just a subset of genes
if (args$test) {
	genes.to.test <- genes.to.test[1:200]
	# genes.to.test <- tail(genes.to.test)
}
# pos.genes.to.test <- all.ts.pairs %>% filter(set=='positive') %>% pull(target.genes)
# neg.genes.to.test <- all.ts.pairs %>% filter(set == "negative") %>% pull(target.genes)

# list for collecting results of model fitting
broom.results <- list()

# if running power analysis for just one value of lambda
if (!is.null(args$lambda)) {

	# load one hot mtx + x_a, x_b, x_ab matrices for lambda
	one.hot.mtx <- h5read(args$h5, paste0("guides/one_hot/", args$lambda))
	xa <- h5read(args$h5, paste0("x/x_a/", args$lambda))
	xb <- h5read(args$h5, paste0("x/x_b/", args$lambda))
	xab <- h5read(args$h5, paste0("x/x_ab/", args$lambda))

	# load list of simulated counts for this value of lambda and all effect sizes
	sim.counts.list <- load_counts(args$lambda, effect.sizes.list) 

	### POWER ANALYSIS 
	pct.correc.vec <- c()
	coeffs.compare.list <- list()

	# iterate through interaction effect sizes
	for (j in 1:length(effect.sizes.list)) {
	    print(paste('effect size =', effect.sizes.list[j]))
	    
	    n.correct <- 0

	    coeffs.compare.list.sub <- list()

	    broom.sublist <- list()

	    k <- 1

	    # iterate through genes to test
	    for (tg in genes.to.test) {
	        print(paste('gene', tg))

	        # get observed counts for gene 
	        obs.counts <- sim.counts.list[[j]][tg,]

	        # if using pseudocount, add pseudocount
	        if (!is.null(args$pseudocount)) {
	            obs.counts <- obs.counts + args$pseudocount
	        }

	        # get X_A, X_B, X_AB for gene
	        xa.gene <- xa[tg,]
	        xb.gene <- xb[tg,]
	        xab.gene <- xab[tg,]
	        
	        # check if min number of cells contains guides for both sites 
	        if (sum(xab.gene !=0)>args$min.cells) {

	        	# collect data in df for model fitting
		        gene.data <- data.frame(tsA = xa.gene, 
		        						tsB = xb.gene,
				                       tsAB = xab.gene,
				                       s.score = cell.cycle.scores$s.scores,
				                       g2m.score = cell.cycle.scores$g2m.scores,
				                       percent.mito = percent.mito,
				                       counts = obs.counts, 
				                       scaling.factor = scaling.factors)
	        
	        	# fit model
		        alt <- glm.nb(counts ~ tsA * tsB + s.score + g2m.score + percent.mito + offset(log(scaling.factor)),
		                      data = gene.data)
				
				# tidy fitted model results		        
		        alt.brm <- broom::tidy(alt)

		        alt.brm$gene <- tg
		        alt.brm$effect <- effect.sizes.list[j]

		        # store list of fitted model results 
		        broom.sublist[[k]] <- alt.brm 

		        # store coefficient data compared to actual data
		        size <- effect.sizes.list[j]
		        if (args$neg) {
		        	size <- 0
		        }

		        coeffs.compare.df <- bind_cols(alt.brm %>% dplyr::select(term, estimate), 
		                                       true = t(coeffs %>% filter(gene == tg) %>%
		                                                rename("tsA" = "betaA",
		                                                       "tsB" = "betaB",
		                                                       "(Intercept)" = "beta0",
		                                                       "s.score" = "beta.s",
		                                                       "g2m.score" = "beta.g2m") %>% 
		                                                dplyr::select("(Intercept)",
		                                                              "tsA",
		                                                              "tsB",
		                                                              "s.score",
		                                                              "g2m.score",
		                                                              "percent.mito") %>%
		                                                mutate(`tsA:tsB`= size)), 
		                                       gene = tg)

	        
		        coeffs.compare.list.sub[[k]] <- coeffs.compare.df
		        k <- k + 1

		        # count number correctly called 
		        if (args$neg) {
		        	# check for cases called as non-significant interaction
		        	if (broom::tidy(alt) %>% filter(term=="tsA:tsB") %>% pull(p.value)>=0.05) {
			            n.correct <- n.correct + 1
			        }
		        } else {
		        	# check for cases called as significant interaction 
		        	if (broom::tidy(alt) %>% filter(term=="tsA:tsB") %>% pull(p.value)<0.05) {
			            n.correct <- n.correct + 1
			        }
		        }
		    }
        }	

        # calculate percent of genes tested that were correctly called 
	    pct.correc.vec[j] <- n.correct/(k-1)

	    # rbind df lists
	    coeffs.compare.list[[j]] <- do.call(rbind,coeffs.compare.list.sub)
	    broom.results[[j]] <- do.call(rbind, broom.sublist)
	}

	### plot power curves based on percentage of correct hits 
	print('calculating pct detected and making plot df for power curves')
	# pct.correct <- sapply(n.correct.list, function(x) {x/length(genes.to.test)})
	power.plotdf <- data.frame(effect.size = effect.sizes.list, 
								pct.detected = pct.correc.vec,
								lambda = args$lambda)

	# power.curve <- ggplot(power.plotdf, aes(x = effect.size, y = pct.detected)) + 
	# 					geom_point() + geom_line() +
	# 					theme_clasic() + 
	# 					theme(text = element_text(size = 16))


	# # save plot as pdf
	# pdf(file.path(args$out, paste0("lambda", args$lambda, "_power_curves.pdf")))
	# print(power.curve)
	# dev.off()

	# # save plot as png
	# png(file.path(args$out, paste0("lambda", args$lambda,"_power_curves.png")), 
	# 	res = 300, units = "in", width = 10, height = 6)
	# print(power.curve)
	# dev.off()

	if (args$neg) {
		print('saving results (NEG)')
		# save power.plotdf (in case we want to work with it later...)
		write.csv(power.plotdf, file.path(args$out, 
					paste0("lambda", args$lambda, "_power_plotdf_NEG_filt", args$min.cells, ".csv")), 
			quote=F, row.names = F)

		# save df of broom results
		saveRDS(broom.results, file.path(args$out, 
					paste0("lambda", args$lambda, "_broom_list_NEG_filt", args$min.cells, ".rds")))

		### export table of est vs. true coeffs for downstream analysis
		saveRDS(coeffs.compare.list, file.path(args$out, 
					paste0("lambda", args$lambda,"_coeffs_list_NEG_filt", args$min.cells, ".rds")))
	} else {
		print('saving results (POS)')
		# save power.plotdf (in case we want to work with it later...)
		write.csv(power.plotdf, file.path(args$out, 
					paste0("lambda", args$lambda, "_power_plotdf_filt", args$min.cells, ".csv")), 
			quote=F, row.names = F)

		# save df of broom results
		saveRDS(broom.results, file.path(args$out, 
					paste0("lambda", args$lambda, "_broom_list_filt", args$min.cells, ".rds")))

		### export table of est vs. true coeffs for downstream analysis
		saveRDS(coeffs.compare.list, file.path(args$out, 
					paste0("lambda", args$lambda,"_coeffs_list_filt", args$min.cells, ".rds")))
	}

}
