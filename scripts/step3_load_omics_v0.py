#%%
from arangodb_utils import *
import pandas as pd
# --- FUNZIONI DI BASE DI ARANGODB ---
db_connection = setup_arangodb_connection()

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

STUDY = "TCGA-BRCA"
ROOT_PATH = f"../data/omics/{STUDY}/"
ROOT_MAPS = f"../data/omics/maps/"

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
#%%
# from omics_connector import make_mirna_hgcn_map
# make_mirna_hgcn_map()

# %%
# print working directory
import os
print(f"Current working directory: {os.getcwd()}")

# print_df_heads()
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
df= read_data_type(DATA_TYPES[5])
print(df.columns)
df

#%%
df.Ensembl_ID.nunique()
df["Ensembl_ID_nover"] = df["Ensembl_ID"].str.split('.').str[0]
df["Ensembl_ID_nover"].nunique()

#%%
gene_mappings = pd.read_csv(PROBEMAP_PATHS[0], sep='\t')
gene_mappings.gene.nunique()

# count only gene that ends with "P"
number_of_pseudo = gene_mappings[gene_mappings['gene'].str.endswith('P')].gene.nunique()   
print(f"Number of pseudo genes: {number_of_pseudo}")


# count only gene that starts with "MIR"
number_of_mirna = gene_mappings[gene_mappings['gene'].str.startswith('MIR')].gene.nunique()   
print(f"Number of miRNA genes: {number_of_mirna}")

# count only gene that contains with "."
number_of_genes_with_dot = gene_mappings[gene_mappings['gene'].str.contains('\.')].gene.nunique()
print(f"Number of genes containing '.': {number_of_genes_with_dot}")


gene_mappings["gene_nodot"] = gene_mappings["gene"].str.split('.').str[0]
genenodot_num = gene_mappings["gene_nodot"].nunique()
print(f"Number of unique genes without dot: {genenodot_num}")
#%%
# get all annotation from biomart of all gene_mappings["gene"]
from pybiomart import Dataset

dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')
# get all attributes
attributes = dataset.attributes
# remove attributes containig "homolog"
for key in list(attributes.keys()):
    if 'homolog' in key:
        del attributes[key]
    if 'structure' in key:
        del attributes[key]
print(attributes.keys())
#%%
gene_biotypes = dataset.query(attributes=['ensembl_gene_id', 'ensembl_gene_id_version', 'gene_biotype'],
              
            #   filters={'chromosome_name': ['1']}
              )

# show distribution of gene_biotype
gene_biotype_counts = gene_biotypes['Gene type'].value_counts()
print(gene_biotype_counts)
#%%
gene_biotypes["Gene stable ID"].drop_duplicates().shape