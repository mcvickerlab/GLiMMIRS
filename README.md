# GLiMMIRS (**G**eneralized **Li**near **M**odels for **M**easuring **I**nteractions between **R**egulatory **S**equences) ✨

Our paper **Genome-wide analysis of CRISPR perturbations indicates that enhancers act multiplicatively and without epistatic-like interactions** is now available as a preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.04.26.538501v1)! Check it out!

This repo contains scripts that demonstrate the data simulation and model fitting functionalities of GLiMMIRS. We have provided tutorials for running our scripts as well as Jupyter notebooks that visualize simulated data and walk the user through the outputs of our scripts and how to interact with them. 

In this repository, you will find two directories:
- [`src/`](https://github.com/mcvickerlab/GLiMMIRS/tree/simulations/src): contains scripts used for dealing with experimental ([`src/experimental`](https://github.com/mcvickerlab/GLiMMIRS/tree/simulations/src/experimental)) and simulated ([`src/simulations`](https://github.com/mcvickerlab/GLiMMIRS/tree/simulations/src/simulations)) data. 
- [`notebooks/`](https://github.com/mcvickerlab/GLiMMIRS/tree/simulations/notebooks): contains Jupyter notebooks that demonstrate how to interact with outputs of scripts and visualize data

## Tutorial 
### Generating simulated data 
##### For GLiMMIRS-base
[`https://github.com/mcvickerlab/GLiMMIRS/blob/simulations/run_sim_counts.sh`](https://github.com/mcvickerlab/GLiMMIRS/blob/simulations/run_sim_counts.sh) is example shell script for generating simulated data for GLiMMIRS-base, which evaluates the effect of single enhancers acting on a single gene. This example runs [`sim_counts.R`](https://github.com/mcvickerlab/GLiMMIRS/blob/simulations/src/simulations/data/sim_counts.R) with the following parameters:
- 13000 genes (`--genes 13000`)
- 50000 cells (`--cells 50000`)
- 1000 target sites, or putative enhancers ('--targets 1000`)
- 2 gRNAs per target site (`--d 2`)
- Sampling number of unique gRNAs per cells from a Poisson distribution parameterized by $\lambda=15$ (https://github.com/mcvickerlab/GLiMMIRS/blob/e9ee30714f1caa88806069b0ee3c7172bc0bb4db/run_sim_counts.sh#L14) (`--lambda 15`)
- Generating corresponding noisy gRNA efficiency estimates at values of $D=1, 10, 100$ (`--guide_disp 1 10 100`)
- saving outputs to `data/simulated/base` (`--out data/simulated/base`)


Please refer to README files in each directory for further details on scripts and notebooks. 

<!-- Here is a comprehensive overview of the directory structure:
```
notebooks/
├── experimental
│   ├── compare_additive_vs_multiplicative_model.ipynb
│   ├── explore_cell_cycle_scores.ipynb
│   ├── explore_enhancer_distance.ipynb
│   └── explore_guidescan_guide_efficiency_missingness.ipynb
└── simulations
    ├── eval_GLiMMIRS-base_on_sim.ipynb
    ├── plot_power_curves.ipynb
    ├── visualize_sim_base_data.ipynb
    └── visualize_sim_interactions_data.ipynb
src/
├── experimental
│   ├── features
│   │   ├── compute_cell_cycle_scores.R
│   │   ├── compute_multiple_enhancer_guide_count_330_published.R
│   │   ├── compute_multiple_enhancer_guide_count_at_scale.R
│   │   ├── create_cell_guide_matrix.py
│   │   ├── create_guidescan_query.py
│   │   └── filter_guidescan_output.py
│   ├── models
│   │   ├── compare_multiplicative_vs_additive.R
│   │   ├── run_baseline_model_experimental_data_neg_mismatch_gene.R
│   │   ├── run_baseline_model_experimental_data_neg_scrambled_guides.R
│   │   ├── run_baseline_model_experimental_data.R
│   │   ├── run_bootstrap_significant_interactions.R
│   │   ├── run_compare_multiplicative_vs_additive.sh
│   │   ├── run_model_at_scale_enhancer_pairs.R
│   │   ├── run_model_experimental_suppl_data_table_2_enhancer_pairs_no_pseudocount.R
│   │   ├── run_model_experimental_suppl_data_table_2_enhancer_pairs.R
│   │   └── run_permutation_null_interaction_coefficients.R
│   ├── README.md
│   └── visualization
│       ├── plot_bootstrap_dotplot.R
│       ├── plot_cell_cycle_score_distributions.R
│       ├── plot_cell_grna_count_distribution.R
│       ├── plot_enhancer_interaction_distance.R
│       ├── plot_enhancer_pair_count_at_scale.R
│       ├── plot_enhancer_pair_distances.R
│       ├── plot_guide_efficiency_distribution.R
│       ├── plot_outlier_dotplot.R
│       ├── plot_permutation_test_histograms.R
│       ├── plot_pseudocount_coefficients.R
│       ├── plot_qqplot_baseline_model_experimental_data.R
│       ├── plot_qqplot_interaction_term_pvalues_at_scale.R
│       └── plot_volcano_interaction_coefficients.R
└── simulations
    ├── data
    │   ├── run_sim_counts_power_analysis.sh
    │   ├── run_sim_counts.sh
    │   ├── sim_counts_interactions_power_analysis.R
    │   └── sim_counts.R
    ├── models
    │   ├── fit_GLiMMIRS-base_sim.R
    │   ├── GLiMMIRS-int_power_analysis.R
    │   ├── run_fit_GLiMMIRS-base_sim.sh
    │   ├── run_GLiMMIRS-int_power_analysis_NEG.sh
    │   └── run_GLiMMIRS-int_power_analysis.sh
    └── README.md
```
 -->
If you have any questions, please feel free to reach out to one of the creators. We will be happy to do what we can to help! 
- Karthik Guruvayurappan: kag4019 [at] med [dot] cornell [dot] edu
- Jessica Zhou: jlz014 [at] eng [dot] ucsd [dot] edu
