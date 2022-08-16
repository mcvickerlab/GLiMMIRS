import pandas as pd

# read in file containing guide spacer sequences from Gasperini 2019 paper
guide_sequences = pd.read_csv('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_grna_groups.at_scale.txt',
                              sep = '\t',
                              names = ['guide_group', 'guide_spacer'])

guide_sequences['guide_spacer'] = guide_sequences['guide_spacer'] + 'NGG'
guide_sequences['guide_spacer'].to_csv('/iblm/netapp/home/karthik/GuideScan/Gasperini2019/guide_sequences.csv', index = False)