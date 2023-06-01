# `src/simulations`

This directory contains scripts related to data simulation with GLiMMIRS. 

In this directory, you will find two subdirectories:

## [`data/`](https://github.com/mcvickerlab/GLiMMIRS/tree/simulations/src/simulations/data)
This directory contains scripts related to data simulation:

### [`sim_counts.R`](https://github.com/mcvickerlab/GLiMMIRS/blob/simulations/src/simulations/data/sim_counts.R)
This script simulates data for GLiMMIRS-base. It generates data resembling a multiplexed, single-cell RNA-seq CRISPRi screen like the one in Gasperini et al. 2019. Specifically, it simulates a CRISPRi experiment in which putative enhancers are targeted with CRISPRi and gene expression is evaluated to determine whether a putative enhancer has an effect on a gene's expression or not, thereby mapping enhancers to their target genes. The simulated data is stored in an h5 file, which contains the simulated counts in a gene by cell matrix, a mapping of gRNAs to cells stored in a one hot encoded matrix, metadata for the gRNAs in the simulated experiment (e.g. ground truth target gene of putative enhancer targeted by the gRNA, gRNA efficiency), per cell variables for modeling, and ground truth coefficient values that can be used as reference when evaluating the coefficient estimates from GLiMMIRS-base. This script takes several input arguments:
- `--genes`: provide the number of genes to be recorded in the simulated experiment. Default: 13000
- `--targets`: provide the number of candidate enhancers targeted in the simulation. Default: 1000
- `--cells`: provide the number of cells in the simulated experiment. Default: 50000
- `--d`: provide the number of gRNAs targeting each target site (candidate enhancer). Default: 2
- `--lambda`: provide the $\lambda$ parameter defining the Poisson distribution from which the number of guides per cell are randomly drawn. Default: 15
- `--guide_disp`: provide one or more dispersion values ($D$, refer to Methods section of manuscript) to use when simulating estimated guide efficiencies. Smaller values of $D$ yield more noisy estimates of dispersion. Default: NULL (no noisy guide efficiencies estimated)
- `--out`: provide path for saving outputs

Example run (from root directory): 
```bash
Rscript src/simulations/data/sim_counts.R \
--out [outdir] \ # specify where simulated data is saved
--genes 13000 \ # 13000 genes in simulation
--cells 50000 \ # 50000 cells in simulation
--targets 1000 \ # 1000 candidate enhancers targeted in CRISPRi experiment
--lambda 15 \ # defines Poisson distribution for simulating number of gRNAs per cell
--guide_disp 1 10 100 \ # simulate noisy guide efficiencies at D=1,10,100
--d 2 # 2 gRNAs targeting each site
```


Here is a comprehensive overview of the directory structure:
```
├── data
│   ├── sim_counts_interactions_power_analysis.R
│   └── sim_counts.R
├── models
│   ├── fit_GLiMMIRS-base_sim.R
│   ├── GLiMMIRS-int_power_analysis.R
└── README.md
```

These scripts were created by Jessica Zhou. Feel free to reach out if you have any questions: jlz014 [at] eng [dot] ucsd [dot] edu
