# This program creates a FASTA file which contains all of the guide spacer sequences used in the
# Gasperini et al. 2019 paper, to provide them as an input to FlashFry for guide efficiency
# scoring.
#
# Author: Karthik Guruvayurappan

import pandas as pd

# read in file containing guide spacer sequences from Gasperini 2019 paper
guide_sequences = pd.read_csv('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_grna_groups.at_scale.txt',
                              sep = '\t',
                              names = ['guide_group', 'guide_spacer'])

with open('/iblm/netapp/home/karthik/FlashFry/Gasperini2019/guide_sequences.fasta', 'w') as file:

    # write each guide sequence to FASTA file
    for guide_sequence in guide_sequences['guide_spacer']:
        sequence_label = 'guide_' + guide_sequence
        file.write('>' + sequence_label + '\n')
        file.write(guide_sequence + 'GGG' + '\n')
