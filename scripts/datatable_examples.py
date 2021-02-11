from datatable import dt, fread, iread, join, by, rbind, cbind, f
import re
import glob
import os
import numpy as np

# Check starting directory
if not bool(re.search('PharmacoDI$', os.getcwd())):
    try: 
        os.chdir('PharmacoDI')
    except: 
        FileNotFoundError('Please change into the top level PharmacoDI' 
            ' directory to run this script!')

# Path variables
data_dir = 'data/procdata'
output_dir = 'data/demo'
psets = ['CTRPv2', 'FIMM', 'gCSI', 'GDSC_v1',
         'GDSC_v2', 'GRAY', 'UHNBreast', 'CCLE']
pset_tables = {pset: os.listdir(os.path.join(data_dir, pset)) for pset in psets}
# pset_tables: ["dose_response", "drug", "datasets_cells", 
#     "dataset_statistics", "cell", "drug_annotation", "gene_drug", 
#     "profile", "dataset", "mol_cell", "gene_annotation", "dataset_cell", 
#     "experiment", "tissue", "gene"]'

pset_name = psets[3]  # GDSC_v1

# -- Read in a single .csv
experiment = fread(os.path.join(data_dir, pset_name, pset_tables[pset_name][-3], 
    f'*{pset_tables[pset_name][-3]}*.csv'))

# -- Read in multiple .csv files and make a single Frame
dose_response = rbind(*iread(os.path.join(data_dir, pset_name, pset_tables[pset_name][0], 
    f'*{pset_tables[pset_name][0]}*.csv')))

# Can use pattern matching to read in multiple files; ** will match any number of subdirectories
# Should make path parsing code much more compact
all_cell_tables = rbind(*iread(os.path.join(data_dir, '**', 'cell', '*cell.csv')))

# -- Write to csv
dose_response.to_csv(os.path.join(output_dir, f'{pset_tables[pset_name][0]}.csv'))

# -- Select (of the form df[filter, select, ...])
# f is for Frame and references variables within the Frame object (i.e., columns)
dose_response[:, [f.id, f.experiment_id]]

# -- Filter
dose_response[(f.id >= 500) & (f.id <= 1000), :]


# -- Filter using Group By (removing rows with duplicated values in one or more columns)

# -- Join
# Set the key column
# NOTE: Current limitation is that keys need to have the same name

# don't set key on dose_response due to
experiment.names = {'name': 'experiment_id'} # rename columns with dict mapping
experiment.key = 'experiment_id'

# Only does left outer join, but that should be fine for our use case
dose_response_experiment = dose_response[:, :, join(experiment)]
