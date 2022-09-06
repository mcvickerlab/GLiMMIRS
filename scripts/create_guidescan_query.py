import pandas as pd

# read in file containing guide spacer sequences from Gasperini 2019 paper
guide_sequences = pd.read_csv('/iblm/netapp/data1/external/Gasperini2019/gasperini_2019_suppl_table_2.csv')

guide_sequences['Spacer'] = guide_sequences['Spacer'] + 'NGG'
guide_sequences['Spacer'].to_csv('/iblm/netapp/home/karthik/GuideScan/Gasperini2019/guide_sequences.csv', index = False)
