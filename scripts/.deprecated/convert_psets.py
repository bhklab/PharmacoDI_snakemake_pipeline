import PharmacoDI as di
import os
import re
import glob
import pandas as pd
from rpy2.robjects import r, pandas2ri

# Change to the correct directory
if 'scripts' not in os.getcwd():
    os.chdir('PharmacoDI/scripts')

# Download the required PSets
if False:
    di.download_canonical_psets(save_dir='data/rawdata')

# Convert PSet to Python
pandas2ri.activate()

readRDS = r["readRDS"]

pset_files = glob.glob('../data/rawdata/*rds')

pset_file = pset_files[2]
pset = readRDS(pset_file)

# FIXME:: Boolean columns in R data.frame being converted to TRUE=1, FALSE=-2147483648
# Dealing with this in table creation currently ^
pset_py = di.convert_pset_to_py(pset)

# FIXME: Currently there are issues with serializing the python pset; seems like there are some R objects
#  hiding in there
pset_py

pset = pset_py

api_url = "https://www.orcestra.ca/api/psets/canonical"

canonical_names = pd.read_json(api_url).name
name = re.sub('_', '.*', pset.get('annotation').get('name')[0])

annot_dir = os.path.join("..", "data", "metadata")
save_dir = os.path.join("..", "data", "procdata")

home = os.environ['HOME']
data_dir = os.path.join(home, "Development", "DataIngestion", "PharmacoDI", "data", "rawdata")


