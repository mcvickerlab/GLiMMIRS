# GLiMMIRS (**G**eneralized **Li**near **M**odels for **M**easuring **I**nteractions between **R**egulatory **S**equences) :star2:
Modeling enhancer-enhancer interactions using Generalized Linear Models (GLMMs)

```create_guidescan_query.py```: appends 'NGG' to guide RNA spacer sequences from Supplementary Table 2 and outputs file ready to input into GuideScan 2.0 gRNA Sequence Search Tool

```explore_guidescan_guide_efficiency_missingness.ipynb```: Jupyter Notebook to explore missingness in GuideScan 2.0 outputs for guide RNA sequences from the at-scale screen previously published by Gasperini et al.

```filter_guidescan_output.py```: filters guide RNAs from supplementary table 2 to only include enhancer-targeting guides and merges with guide efficiency information

```plot_guide_efficiency_distribution.R```: plots distribution of enhancer-targeting guide RNA efficiencies as a histogram

```explore_cell_cycle_scores.ipynb```: Jupyter notebook that runs cell cycle scoring functions and produces a PCA that separates cells in the at-scale screen based on cell cycle score

```compute_cell_cycle_scores.R```: computes cell cycle scores using Seurat single-cell RNA-sequencing package and generates PCA plot that shows cells separating based on cell cycle score

