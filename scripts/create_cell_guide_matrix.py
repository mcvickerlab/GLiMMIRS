# import data analysis packages
import numpy as np
import pandas as pd

# file and command line packages
import gzip 

# string parsing
import re

# read in phenodata file
colnames = open('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale.phenoData.colnames.txt') \
           .read().splitlines()
phenodata_df = pd.read_csv('/iblm/netapp/home/karthik/gasperini_project/data/phenodata.txt', sep=' ', names=colnames)

# drop nan values from dataframe (cannot parse strings)
phenodata_df = phenodata_df.dropna()

def get_guide_names(guides):
    '''regex function to separate guide names'''
    
    guide_list = re.findall(r"[A-Za-z0-9]*\_TSS|chr[A-Z0-9]{1,2}\.\w*top_two|chr[A-Z0-9]{1,2}\.\w*second_two", guides)
    return guide_list


# get guide names and split guide sequences for each cell 
phenodata_df['guides'] = phenodata_df['gene'].apply(get_guide_names)
phenodata_df['sequences'] = phenodata_df['barcode'].str.split('_')

# select necessary columns from dataframe
phenodata_df = phenodata_df[['cell', 'guides', 'sequences']]

# get list of guide sequences and cell names (slow step)
guide_sequences = pd.Series(phenodata_df['sequences'].sum()).unique()
cell_names = phenodata_df['cell']

def guides_present(cell_guide_list):
    '''helper function to determine which guides from guide sequence list are present in cell'''
    
    return pd.Series(guide_sequences).isin(cell_guide_list).astype(np.int64)


# generate cell guide matrix and transpose matrix so columns are cells
cell_guide_matrix = phenodata_df['sequences'].apply(guides_present)
cell_guide_matrix.index = cell_names
cell_guide_matrix.columns = guide_sequences
cell_guide_matrix = cell_guide_matrix.T

# write output to h5 file (for storing large data)
cell_guide_matrix.to_hdf('/iblm/netapp/home/karthik/crisprQTL/gasperini_data/cell_guide_matrix.h5', key='df', mode='w')
