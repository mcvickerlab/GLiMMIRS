# crisprQTL
Modeling enhancer-enhancer interactions using Generalized Linear Mixed Models (GLMMs)

The goal of this research project is to reanalyze an existing published Perturb-seq dataset from Gasperini et al. 
In their paper, they determined enhancer-gene pairs using a negative binomail generalized linear model with
several covariates, including gRNA count [TODO]. The goal of our project is to build on this existing model, which
is implemented in Monocle2, by modeling guide efficiency as a random effect, making this model a generalized linear
mixed model. Additionally, we aim to add several covariates. Using this improved model, we hope to also model 
enhancer-enhancer interactions to determine their nature as either synergistic, logistic, or additive. 

