
#%%
# Aanalys patients from TCGA-BRCA-counts file
import pandas as pd

STUDY = "TCGA-BRCA"
file_path = f"../data/omics/{STUDY}/{STUDY}.star_counts.tsv.gz"
df = pd.read_csv(file_path, sep='\t', compression='gzip', low_memory=False, nrows=3)

sample_type_codes = pd.read_csv("../data/omics/metadata/TCGA_sample_type_codes.csv")

""" 
Code,Sample Type,Description,Short Code
01,Primary Solid Tumor,Primary solid tumor,TP
02,Recurrent Solid Tumor,Recurrent solid tumor,TR
"""


def patient_code_parser(code):
    """Extract patient code from TCGA barcode.
    example code 'TCGA-C8-A274-01A'
    I need just sample type '01' 
    """
    parts = code.split('-')
    if len(parts) >= 3:
        sample_type = parts[3][:2]  # Get the first two characters of the fourth part
        return sample_type
    return None

# get all sample type from df.columns
sample_types = {col: patient_code_parser(col) for col in df.columns if col.startswith("TCGA")}
# show unique sample types and counts

for i in set(sample_types.values()):
    descriptinion = sample_type_codes[sample_type_codes['Code'] == int(i)]['Sample Type'].values

    print(f"Sample Type {i}: {list(sample_types.values()).count(i)} samples, type: {descriptinion}")
