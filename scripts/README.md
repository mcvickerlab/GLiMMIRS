I attempted to compute guide efficiencies using two different methods. 

The first approach involved using FlashFry, since this was the method initially used in the Gasperini et al. paper. I used the protospacer sequence provided in the dataset, created a FlashFry index for the hg19 genome, converted the guide library into a FASTA file, and then computed off targets and scores. However, only the off-target scoring metrics work, since the on-target scoring metrics (more similar to efficiency) also require flanking sequence. 

The second approach was to use GuideScan 2.0, and the new gRNA sequence search, which takes as input a list of protospacers + NGG (PAM), and outputs the coordinates, efficiency, and specificity of the guide sequence. Mapping to hg38 should not create significant issues, since the assemblies are mostly similar. This approach was smoother overall, but there are roughly 3000 guide sequences out of roughly 13000 that did not find a match in the GuideScan database. 

I plotted the distribution of guide efficiencies derived from the second method using BoutrosLab.plotting.general.

