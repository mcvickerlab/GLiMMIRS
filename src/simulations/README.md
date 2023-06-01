# [src/simulations]

This directory contains scripts related to data simulation with GLiMMIRS. 

In this directory, you will find two subdirectories:

## [`data/`](https://github.com/mcvickerlab/GLiMMIRS/tree/simulations/src/simulations/data)
This directory contains scripts related to data simulation:

### [`sim_counts.R`](https://github.com/mcvickerlab/GLiMMIRS/blob/simulations/src/simulations/data/sim_counts.R)
This script simulates data for GLiMMIRS-base. It generates data resembling a multiplexed, single-cell RNA-seq CRISPRi screen like the one in Gasperini et al. 2019. Specifically, it simulates a CRISPRi experiment in which putative enhancers are targeted with CRISPRi and gene expression is evaluated to determine whether a putative enhancer has an effect on a gene's expression or not, thereby mapping enhancers to their target genes. The simulated data is stored in an h5 file, which contains the simulated counts in a gene by cell matrix, a mapping of gRNAs to cells stored in a one hot encoded matrix, metadata for the gRNAs in the simulated experiment (e.g. ground truth target gene of putative enhancer targeted by the gRNA, gRNA efficiency), per cell variables for modeling, and ground truth coefficient values that can be used as reference when evaluating the coefficient estimates from GLiMMIRS-base. This script takes several input arguments:

Here is a comprehensive overview of the directory structure:
```
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

These scripts were created by Jessica Zhou. Feel free to reach out if you have any questions: jlz014 [at] eng [dot] ucsd [dot] edu
