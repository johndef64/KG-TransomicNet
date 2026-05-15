"""
mondo_mapping_manual.py
=======================
PROVENANCE / EXPLORATORY SCRIPT — do NOT run as-is.

This is the working notebook used to build `data/mappings/mondo_to_primarydisease.tsv`
manually, by inspecting MONDO labels against the TCGA `disease_type.project` values.
The file is shipped with the repository; this script is kept only to document
how those 33 rows were assigned.

The original script connects to ArangoDB on import (so cannot be run in
isolation) and contains scratch exploration cells. It is intentionally NOT
refactored into a clean CLI — its only output is the curated TSV that is
already in `data/mappings/`.
"""

#%%
from sqlalchemy import null
from arangodb_utils import *
import pandas as pd

# NOTE: this line connects to ArangoDB at import time.
# Commented out to allow the file to be imported / inspected without a DB.
# db_connection = setup_arangodb_connection("PKT_test10000")

# funzioni per caricare multi omics datasets
TCGA_STUDIES = [
    "TCGA-LAML",  # Acute Myeloid Leukemia
    "TCGA-ACC",   # Adrenocortical Cancer
    "TCGA-CHOL",  # Bile Duct Cancer
    "TCGA-BLCA",  # Bladder Cancer
    "TCGA-BRCA",  # Breast Cancer
    "TCGA-CESC",  # Cervical Cancer
    "TCGA-COAD",  # Colon Cancer
    "TCGA-UCEC",  # Endometrioid Cancer
    "TCGA-ESCA",  # Esophageal Cancer
    "TCGA-GBM",   # Glioblastoma
    "TCGA-HNSC",  # Head and Neck Cancer
    "TCGA-KICH",  # Kidney Chromophobe
    "TCGA-KIRC",  # Kidney Clear Cell Carcinoma
    "TCGA-KIRP",  # Kidney Papillary Cell Carcinoma
    "TCGA-DLBC",  # Large B-cell Lymphoma
    "TCGA-LIHC",  # Liver Cancer
    "TCGA-LGG",   # Lower Grade Glioma
    "TCGA-LUAD",  # Lung Adenocarcinoma
    "TCGA-LUSC",  # Lung Squamous Cell Carcinoma
    "TCGA-SKCM",  # Melanoma
    "TCGA-MESO",  # Mesothelioma
    "TCGA-UVM",   # Ocular melanomas
    "TCGA-OV",    # Ovarian Cancer
    "TCGA-PAAD",  # Pancreatic Cancer
    "TCGA-PCPG",  # Pheochromocytoma & Paraganglioma
    "TCGA-PRAD",  # Prostate Cancer
    "TCGA-READ",  # Rectal Cancer
    "TCGA-SARC",  # Sarcoma
    "TCGA-STAD",  # Stomach Cancer
    "TCGA-TGCT",  # Testicular Cancer
    "TCGA-THYM",  # Thymoma
    "TCGA-THCA",  # Thyroid Cancer
    "TCGA-UCS",   # Uterine Carcinosarcoma
]

DATA_TYPES = [
    "gene-level_ascat3",       # Copy Number (Gene Level)
    "allele_cnv_ascat3",       # Allele-specific Copy Number Segment
    #"methylation450",          # DNA Methylation - Illumina Human Methylation
    "methylation27",           # DNA Methylation - Illumina Human Methylation
    "somaticmutation_wxs",     # Somatic Mutation
    "mirna",                   # miRNA Expression
    "star_counts",             # Gene Expression (STAR - counts)
    "star_fpkm",               # Gene Expression (STAR - FPKM)
    "star_tpm",                # Gene Expression (STAR - TPM)
    "protein",                 # Protein Expression
    "clinical",                # Clinical Data
    "survival",                # Survival Data
]

########
STUDY = "TCGA-BRCA" # CHANGE HERE TO SELECT STUDY
########
ROOT_PATH = f"../data/omics/{STUDY}/"
ROOT_MAPS = f"../data/mappings/"

PROBEMAP_PATHS = [
    # gene code annotation files
    f"{ROOT_MAPS}gencode.v36.annotation.gtf.gene.probemap",
    # methylation probe map files
    f"{ROOT_MAPS}HM450.hg38.manifest.gencode.v36.probeMap",
    f"{ROOT_MAPS}HM27.hg38.manifest.gencode.v36.probeMap",
]

OMICS_PATHS = {
    data_type: f"{ROOT_PATH}{STUDY}.{data_type}.tsv.gz"
    for data_type in DATA_TYPES
}
METADATA_PATHS = {
    data_type: f"{ROOT_PATH}{STUDY}.{data_type}.tsv.json"
    for data_type in DATA_TYPES
}
OMICS_PATHS

PANCAN_PATHS = {
    "phenotype": "../data/omics/PANCAN/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
    "subtypes": "../data/omics/PANCAN/TCGASubtype.20170308.tsv.gz",
}


def print_df_heads():
    """Print the heads of all omics dataframes for the selected study.
    print jys forst 5 columns and first 3 rows of each dataframe.
    """
    for data_type, file_path in OMICS_PATHS.items():
        if os.path.exists(file_path):
            print(f"\nData Type: {data_type}")
            df = pd.read_csv(file_path, sep='\t', compression='gzip', low_memory=False, nrows=3)
            # print(f"Stats: {df.shape[0]} rows x {df.shape[1]} columns")
            print(df.iloc[:3, :5])  # print first 5 columns and first 3 rows

        else:
            print(f"\nData Type: {data_type} - File not found: {file_path}")

def read_data_type(data_type):
    """Read a specific omics data type dataframe for the selected study."""
    file_path = OMICS_PATHS.get(data_type)
    if file_path and os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t', compression='gzip', low_memory=False)
        return df
    else:
        print(f"File not found for data type {data_type}: {file_path}")
        return None


#%%  Load PKT NODELABELS ===============================
pkt_nodelabels_path = "../data/pkt/builds/v3.0.2/PKT_NodeLabels_with_metadata_v3.0.2.csv"
pkt_nodelabels = pd.read_csv(pkt_nodelabels_path, dtype=str)
pkt_nodelabels
#%%
pkt_nodelabels['bioentity_type'].value_counts()
"""
in teoria sono mappabili su PKT
bioentity_type:
- gene --> entrez_id
- protein --> uniprot_id
- variant --> rsid
- disease --> mondo_id
- phenotype --> hp_id
"""

# serch in db_connection for collection "nodes"in nodes containing "tumor or cancer ion labels"
get_nodes_by_pattern(db_connection, "nodes", "label", "%tumor%")
#%%
nodes = get_nodes_by_pattern(db_connection, "nodes", "class_code", "%HP%")
print(f"Found {len(nodes)} nodes with class_code containing 'HP'")
# example of nodes properties
for node in nodes[:5]:
    for key, value in node.items():
        print(f"{key}: {value}")
    print("-----")

#%% serch "vital" in HP labels
from tqdm import tqdm
nodes_query = get_nodes_by_pattern(db_connection, "nodes", "label", "% vital%")
for n in tqdm(nodes_query):
    print(f"{n['class_code']}: {n['label']}")
#%%  Clinical_df SUB-DATAFRAMES
clinical_df = read_data_type("clinical")
def make_clinical_sub_df(clinical_df, suffix):
    # la prima colonna di sub_df deve essere sempre case_id
    """Create a sub-dataframe from clinical_df based on column suffix."""
    sub_df = clinical_df[[col for col in clinical_df.columns if col.endswith(suffix)]]
    sub_df.insert(0, 'sample', clinical_df['sample'])
    # remove the suffix from column names
    sub_df.columns = [col.replace(f".{suffix}", "") for col in sub_df.columns]
    return sub_df

df_tissue_source_site = make_clinical_sub_df(clinical_df, 'tissue_source_site')
df_diagnoses = make_clinical_sub_df(clinical_df, 'diagnoses')
df_demographic = make_clinical_sub_df(clinical_df, 'demographic')
df_annotations = make_clinical_sub_df(clinical_df, 'annotations')
df_project = make_clinical_sub_df(clinical_df, 'project')
df_samples = make_clinical_sub_df(clinical_df, 'samples')
df_project
#%%  MONDO ID MAPPING to Disease Types in PANCAN CLINICAL DATA

# load PANCAN data
pancan_phenotype_df = pd.read_csv(PANCAN_PATHS["phenotype"], sep='\t', compression='gzip', low_memory=False)
pancan_phenotype_df.describe(include='object').T
#%%
pancan_subtypes_df = pd.read_csv(PANCAN_PATHS["subtypes"], sep='\t', compression='gzip', low_memory=False)
pancan_subtypes_df.describe(include='object').T

# load mondo mapping
mondo_map_path = f"{ROOT_MAPS}mondo_to_primarydisease.tsv"
mondo_map_df = pd.read_csv(mondo_map_path, sep='\t', dtype=str)

# find mondo nodes in nodelabels
mondo_nodes = pkt_nodelabels[pkt_nodelabels['bioentity_type'] == "disease"]
mondo_nodes['class_code'].value_counts()
mondo_nodes["entity"] 

# merge mondo_map_df with mondo_nodes on entity, anche check missing ones
merged_mondo = pd.merge(mondo_map_df, mondo_nodes, on="entity", how="left", indicator=True)
merged_mondo._merge.value_counts()
merged_mondo[merged_mondo['_merge'] == 'left_only']
DISEASE_MONDO_MAP = merged_mondo[['_primary_disease', 'entity', 'class_code']]
DISEASE_MONDO_MAP
#%%

# serch 'Squamous Cell Neoplasms' in pancan_subtypes_df disease_type column
cols = ['sampleID', 
         'Subtype_mRNA', 'Subtype_DNAmeth', 'Subtype_protein',
       'Subtype_miRNA', 'Subtype_CNA', 'Subtype_Integrative', 'Subtype_other',
       'Subtype_Selected']
for col in cols:
    print(f"Value counts for column: {col}")
    display(pancan_subtypes_df[col].value_counts())

#%%
# for pancan_phenotype_df disease_type column
cols = pancan_phenotype_df.columns.to_list()
for col in cols:
    print(f"Value counts for column: {col}")
    display(pancan_phenotype_df[col].value_counts())
#%%  BRCA Example mapping disease types to MONDO IDs

disease_types = {'Squamous Cell Neoplasms', 'Adnexal and Skin Appendage Neoplasms', 'Epithelial Neoplasms, NOS', 'Complex Epithelial Neoplasms', 'Fibroepithelial Neoplasms', 'Cystic, Mucinous and Serous Neoplasms', 'Basal Cell Neoplasms', 'Adenomas and Adenocarcinomas', 'Ductal and Lobular Neoplasms'}
# il tuo compito è fare una ricerca per ogni termine e darmi il mapping con l'id dell'ontologia MONDO per ognuno di essi, che sia corretto, lo volgio in un json

nodes_query = get_nodes_by_pattern(db_connection, "nodes", "label", "%Fibroepithelial%")
for n in nodes_query:
    print(f"{n['class_code']}: {n['label']}, _id: {n['_key']}")
    
# nodes_query = get_nodes_by_ids(db_connection, "nodes", ["MONDO_0002297"])
# for n in nodes_query:
#     print(f"{n['class_code']}: {n['label']}")
#%%

mondo_mapping = {
  "Squamous Cell Neoplasms": "MONDO:0002532",
  "Adnexal and Skin Appendage Neoplasms": "MONDO:0002297",
  "Epithelial Neoplasms, NOS": "MONDO:0005626",
  "Complex Epithelial Neoplasms": "MONDO:0005626",
  "Fibroepithelial Neoplasms": "MONDO:0021046",
  "Cystic, Mucinous and Serous Neoplasms": "MONDO:0006720",
  "Basal Cell Neoplasms": "MONDO:0020799",
  "Adenomas and Adenocarcinomas":  "MONDO:0004972|MONDO:0004970",
  "Ductal and Lobular Neoplasms": 
  "MONDO:0004953|MONDO:0002486"
#   }
}
mondo_ids = mondo_mapping.values()

for id in mondo_ids:
    id = id.replace(":", "_")
    id = id.split('|')[0]  # take first id only
    nodes_query = get_nodes_by_pattern(db_connection, "nodes", "_key", f"{id}")
    for n in nodes_query:
        if n['class_code'] == "MONDO":
            print(f"{n['class_code']}: {n['label']}")

# for dt in disease_types:
#     dt = dt.split(',')[0]  # take first word only
#     dt = dt.rstrip('s')  
#     nodes_query = get_nodes_by_pattern(db_connection, "nodes", "label", f"%{dt}%")
#     for n in nodes_query:
#         print(f"{n['class_code']}: {n['label']}")

#%%
# from omics_utils import make_mirna_hgcn_map
# make_mirna_hgcn_map()
pkt_nodelabels[pkt_nodelabels['bioentity_type'] == "sequence"].class_code.value_counts()

