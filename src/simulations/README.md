# `src/simulations`

This directory contains scripts related to data simulation with GLiMMIRS. 

In this directory, you will find two subdirectories:

## [`data/`](https://github.com/mcvickerlab/GLiMMIRS/tree/simulations/src/simulations/data)
This directory contains scripts related to data simulation:

### [`sim_counts.R`](https://github.com/mcvickerlab/GLiMMIRS/blob/simulations/src/simulations/data/sim_counts.R)
This script simulates data for GLiMMIRS-base. It generates data resembling a multiplexed, single-cell RNA-seq CRISPRi screen like the one in Gasperini et al. 2019. Specifically, it simulates a CRISPRi experiment in which putative enhancers are targeted with CRISPRi and gene expression is evaluated to determine whether a putative enhancer has an effect on a gene's expression or not, thereby mapping enhancers to their target genes. The simulated data is stored in an h5 file, which contains the simulated counts in a gene by cell matrix, a mapping of gRNAs to cells stored in a one hot encoded matrix, metadata for the gRNAs in the simulated experiment (e.g. ground truth target gene of putative enhancer targeted by the gRNA, gRNA efficiency), per cell variables for modeling, and ground truth coefficient values that can be used as reference when evaluating the coefficient estimates from GLiMMIRS-base. Specifically, this script takes into account the following terms:
- intercept ($\beta_0$)
- gRNA efficiency ($X_{perturb}, \beta_{enhancer}$)
- cell cycle scores ($X_S, X_{G2M}, \beta_S, \beta_{G2M}$) 
- percentage of mitochondrial DNA ($X_{mito}, \beta_{mito}$)
 
This script takes several input arguments:
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
##### Outputs
This script will produce an h5 file (`sim.h5`) that contains all of the simulated data and relevant metadata for evaluating the performance of GLiMMIRS-base applied to this simulated dataset. 
###### Notes:
- line 149 contains hardcoded path to `/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5`, from which the cell cycle scores corresponding to the experimental data is located. We fit a distribution to the observed cell cycle scores and use it to generated simulated cell cycle scores resembling the empirical data. 

Throughout this script, we have hardcoded values to define the distributions from which coefficients and variables are sampled. These values have been determined via maximum likelihood estimation of the experimental data from Gasperini et al. 2019. Therefore, this script is designed to simulate data resembling the pilot experiment from that paper. In order to simulate data resembling another CRISPR regulatory screen, the user may need to modify the script to define the sampling distributions differently, or use different statistical distributions altogether as appropriate. Here are some of the places in the script where you might want to change: 

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
