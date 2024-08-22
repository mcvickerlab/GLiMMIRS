### this code runs a power analysis fitting GLiMMIRS-int to simulated dataset 

library(argparse)
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
					help = "specify which MOI to do power analysis for; if NULL analyze all")
parser$add_argument("--neg", action = "store_true",
					help = "if true, test negative pairs, e.g. enhancers without interactions (DEFAULT=test positive pairs")
parser$add_argument("--neg_effect", action = "store_true",
					help = "if true, interaction effect sizes are negative")
args <- parser$parse_args()

print(args)

# show contents of h5 file
h5ls(args$h5)

### define helper functions 
load_counts <- function(lambda, effect.sizes, neg = FALSE) {
	# lambda should be fixed value
	# effect.sizes is a list of effect sizes
	sim.counts.list <- list()

	for (j in 1:length(effect.sizes)) {
		set.name <- paste0("lambda", lambda, "_size", effect.sizes[j])
		if (neg) {
			set.name <- paste0("lambda", lambda, "_size-NEG", effect.sizes[j])
		}
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

# if running power analysis for just one value of lambda
if (!is.null(args$lambda)) {
	one.hot.mtx <- h5read(args$h5, paste0("guides/one_hot/", args$lambda))
	xa <- h5read(args$h5, paste0("x/x_a/", args$lambda))
	xb <- h5read(args$h5, paste0("x/x_b/", args$lambda))
	xab <- h5read(args$h5, paste0("x/x_ab/", args$lambda))

	# load list of simulated counts for this value of lambda and all effect sizes
	sim.counts.list <- load_counts(args$lambda, effect.sizes.list, args$neg_effect) 

	### POWER ANALYSIS 
	n.correct.list <- c()
	coeffs.compare.list <- list()
	broom.list <- list()

	for (j in 1:length(effect.sizes.list)) {
	    print(paste('effect size =', effect.sizes.list[j]))
	    
	    n.correct <- 0
	    
	    coeffs.compare.list.sub <- list()
	    broom.list.sub <- list()
	    k <- 1

	    effect.size <- effect.sizes.list[j]
	    if (args$neg_effect) {
	    	effect.size <- -1*effect.size
	    }

	    for (tg in genes.to.test) {
	        print(paste('gene', tg))
	        obs.counts <- sim.counts.list[[j]][tg,]
	        if (!is.null(args$pseudocount)) {
	            obs.counts <- obs.counts + args$pseudocount
	        }

	        xa.gene <- xa[tg,]
	        xb.gene <- xb[tg,]
	        xab.gene <- xab[tg,]
	        
	        gene.data <- data.frame(tsA = xa.gene, 
	                       tsB = xb.gene,
	                       tsAB = xab.gene,
	                       s.score = cell.cycle.scores$s.scores,
	                       g2m.score = cell.cycle.scores$g2m.scores,
	                       percent.mito = percent.mito,
	                       counts = obs.counts, 
	                       scaling.factor = scaling.factors)
	        
	        alt <- glm.nb(counts ~ tsA * tsB + s.score + g2m.score + percent.mito + offset(log(scaling.factor)),
	                      data = gene.data)
	        
	        alt.brm <- broom::tidy(alt)
	        alt.brm$gene <- tg
			broom.list.sub[[k]] <- alt.brm

	        coeffs.compare.df <- NULL
	        if (args$neg) {
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
	                                                mutate(`tsA:tsB`= 0)), 
	                                       gene = tg)
	        } else {
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
			                                        mutate(`tsA:tsB`= effect.size)), 
			                               gene = tg)
	        }

	        
	        coeffs.compare.list.sub[[k]] <- coeffs.compare.df
	        k <- k + 1

	        if (args$neg) {
	        	if (broom::tidy(alt) %>% filter(term=="tsA:tsB") %>% pull(p.value)>=0.05) {
		            n.correct <- n.correct + 1
		        }
	        } else {
	        	if (broom::tidy(alt) %>% filter(term=="tsA:tsB") %>% pull(p.value)<0.05) {
		            n.correct <- n.correct + 1
		        }
	        }

	    }
	    n.correct.list[j] <- n.correct
	    coeffs.compare.list[[j]] <- do.call(rbind, coeffs.compare.list.sub)
	    broom.list[[j]] <- do.call(rbind, broom.list.sub)

	}

	### plot power curves based on percentage of correct hits 
	print('calculating pct detected and making plot df for power curves')
	pct.correct <- sapply(n.correct.list, function(x) {x/length(genes.to.test)})
	power.plotdf <- data.frame(effect.size = effect.sizes.list, 
								pct.detected = pct.correct,
								lambda = args$lambda)
	if (args$neg_effect) {
		power.plotdf <- data.frame(effect.size = -1*effect.sizes.list, 
							pct.detected = pct.correct,
							lambda = args$lambda)
	}

	if (args$neg) {
		print('saving results (NEG)')
		# save power.plotdf (in case we want to work with it later...)
		write.csv(power.plotdf, file.path(args$out, paste0("lambda", args$lambda, "_power_plotdf_NEG.csv")), 
			quote=F, row.names = F)
		
		# save list of tidy model results (with p-vals)
		saveRDS(broom.list, file.path(args$out, paste0("lambda", args$lambda, "_tidy_alt_NEG.rds")))

		### export table of est vs. true coeffs for downstream analysis
		saveRDS(coeffs.compare.list, file.path(args$out, paste0("lambda", args$lambda,"_coeffs_list_NEG.rds")))
	} else {
		print('saving results (POS)')
		# save power.plotdf (in case we want to work with it later...)
		write.csv(power.plotdf, file.path(args$out, paste0("lambda", args$lambda, "_power_plotdf.csv")), 
			quote=F, row.names = F)
		
		# save list of tidy model results (with p-vals)
		saveRDS(broom.list, file.path(args$out, paste0("lambda", args$lambda, "_tidy_alt.rds")))
		
		### export table of est vs. true coeffs for downstream analysis
		saveRDS(coeffs.compare.list, file.path(args$out, paste0("lambda", args$lambda,"_coeffs_list.rds")))
	}

} else {
	# store df for plotting power curves for each value of lambda
	power.plotdf.list <- list()
	i <- 1

	for (lambda in lambda.list) {
		print(paste("lambda=", lambda))
		one.hot.mtx <- h5read(args$h5, paste0("guides/one_hot/", lambda))
		xa <- h5read(args$h5, paste0("x/x_a/", lambda))
		xb <- h5read(args$h5, paste0("x/x_b/", lambda))
		xab <- h5read(args$h5, paste0("x/x_ab/", lambda))

		# load list of simulated counts for this value of lambda and all effect sizes
		sim.counts.list <- load_counts(lambda, effect.sizes.list, args$neg_effect) 

		### POWER ANALYSIS 
		n.correct.list <- c()
		coeffs.compare.list <- list()
		broom.list <- list()

		for (j in 1:length(effect.sizes.list)) {
		    print(paste('effect size =', effect.sizes.list[j]))
		    
		    n.correct <- 0
		    
		    coeffs.compare.list.sub <- list()
		    broom.list.sub <- list()
		    k <- 1

		    effect.size <- effect.sizes.list[j]
		    if (args$neg_effect) {
		    	effect.size <- -1*effect.size
		    }

		    for (tg in genes.to.test) {
		        print(paste('gene', tg))
		        obs.counts <- sim.counts.list[[j]][tg,]
		        if (!is.null(args$pseudocount)) {
		            obs.counts <- obs.counts + args$pseudocount
		        }

		        xa.gene <- xa[tg,]
		        xb.gene <- xb[tg,]
		        xab.gene <- xab[tg,]
		        
		        gene.data <- data.frame(tsA = xa.gene, 
		                       tsB = xb.gene,
		                       tsAB = xab.gene,
		                       s.score = cell.cycle.scores$s.scores,
		                       g2m.score = cell.cycle.scores$g2m.scores,
		                       percent.mito = percent.mito,
		                       counts = obs.counts, 
		                       scaling.factor = scaling.factors)
		        
		        alt <- glm.nb(counts ~ tsA * tsB + s.score + g2m.score + percent.mito + offset(log(scaling.factor)),
		                      data = gene.data)
		        
		        alt.brm <- broom::tidy(alt)
		        alt.brm$gene <- tg
				broom.list.sub[[k]] <- alt.brm

		        coeffs.compare.df <- NULL
		        if (args$neg) {
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
		                                                mutate(`tsA:tsB`= 0)), 
		                                       gene = tg)
		        } else {
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
				                                        mutate(`tsA:tsB`= effect.size)), 
				                               gene = tg)
		        }

		        
		        coeffs.compare.list.sub[[k]] <- coeffs.compare.df
		        k <- k + 1

		        if (args$neg) {
		        	if (broom::tidy(alt) %>% filter(term=="tsA:tsB") %>% pull(p.value)>=0.05) {
			            n.correct <- n.correct + 1
			        }
		        } else {
		        	if (broom::tidy(alt) %>% filter(term=="tsA:tsB") %>% pull(p.value)<0.05) {
			            n.correct <- n.correct + 1
			        }
		        }

		    }
		    n.correct.list[j] <- n.correct
		    coeffs.compare.list[[j]] <- do.call(rbind, coeffs.compare.list.sub)
		    broom.list[[j]] <- do.call(rbind, broom.list.sub)

		}

		### plot power curves based on percentage of correct hits 
		print('calculating pct detected and making plot df for power curves')
		pct.correct <- sapply(n.correct.list, function(x) {x/length(genes.to.test)})
		power.plotdf <- data.frame(effect.size = effect.sizes.list, 
									pct.detected = pct.correct,
									lambda = lambda)
		if (args$neg_effect) {
			power.plotdf <- data.frame(effect.size = -1*effect.sizes.list, 
								pct.detected = pct.correct,
								lambda = lambda)
		}

		power.plotdf.list[[i]] <- power.plotdf
		i <- i + 1

		if (args$neg) {
			print('saving results (NEG)')		
			# save list of tidy model results (with p-vals)
			saveRDS(broom.list, file.path(args$out, paste0("lambda", args$lambda, "_tidy_alt_NEG.rds")))

			### export table of est vs. true coeffs for downstream analysis
			saveRDS(coeffs.compare.list, file.path(args$out, paste0("lambda", args$lambda,"_coeffs_list_NEG.rds")))
		} else {
			print('saving results (POS)')			
			# save list of tidy model results (with p-vals)
			saveRDS(broom.list, file.path(args$out, paste0("lambda", args$lambda, "_tidy_alt.rds")))
			
			### export table of est vs. true coeffs for downstream analysis
			saveRDS(coeffs.compare.list, file.path(args$out, paste0("lambda", args$lambda,"_coeffs_list.rds")))
		}
	}

	# save power.plotdf.list as a big dataframe for plotting power curves later
	if (args$neg) {
		print('writing power curve data to file (NEG)')
		write.csv(do.call(rbind, power.plotdf.list), 
			file.path(args$out, "power_analysis_plot_data_NEG.csv"), 
			quote = FALSE, row.names = FALSE)
	} else {
		print('writing power curve data to file (POS)')
		write.csv(do.call(rbind, power.plotdf.list), 
			file.path(args$out, "power_analysis_plot_data.csv"), 
			quote = FALSE, row.names = FALSE)		
	}
}
