# This script reads in the output of running the at-scale screen guide RNAs
# through GuideScan and filters those outputs to only include
# enhancer-targeting guides.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# load in output file from GuideScan and filter 'NGG' from end of each spacer sequence
guidescan_output = pd.read_csv('/iblm/netapp/home/karthik/GuideScan/Gasperini2019/guidescan_output.csv')
guidescan_output['guide_spacer'] = guidescan_output['gRNA'].apply(lambda x: x[:-3])

# load in guide RNA spacer sequences from supplementary table 2
guide_sequences = pd.read_csv('/iblm/netapp/data1/external/Gasperini2019/gasperini_2019_suppl_table_2.csv')

# merge data frames together on spacer sequence to connect efficiency with guide metadata
guide_sequences = guide_sequences.merge(guidescan_output, left_on = 'Spacer', right_on = 'guide_spacer')

# filter out TSS guides
guide_sequences = guide_sequences[guide_sequences['Category'] != 'TSS']

# filter out NTC guides
guide_sequences = guide_sequences[guide_sequences['Category'] != 'NTC']

# filter out globin locus targeting guides
guide_sequences = guide_sequences[guide_sequences['Category'] != 'Positive_control_to_globin_locus']

# write enhancer-targeting guide sequences to output file
guide_sequences.to_csv('/iblm/netapp/home/karthik/GuideScan/Gasperini2019/enhancer_guidescan_output.csv', index = False)

