#%%
from sqlalchemy import null
from arangodb_utils import *
import pandas as pd
# --- FUNZIONI DI BASE DI ARANGODB ---
db_connection = setup_arangodb_connection("PKT_test10000")

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

#%%
# serch "vital" in HP labels
nodes_query = get_nodes_by_pattern(db_connection, "nodes", "label", "% female%")
for n in nodes_query:
    print(f"{n['class_code']}: {n['label']}")
#%%
# volgio codice python per creare diversi datasets salvati in csv che risolvono la struttura nested di clinica_df 
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
#%%

disease_types = {'Squamous Cell Neoplasms', 'Adnexal and Skin Appendage Neoplasms', 'Epithelial Neoplasms, NOS', 'Complex Epithelial Neoplasms', 'Fibroepithelial Neoplasms', 'Cystic, Mucinous and Serous Neoplasms', 'Basal Cell Neoplasms', 'Adenomas and Adenocarcinomas', 'Ductal and Lobular Neoplasms'}
# il tuo compito è fare una ricerca per ogni termine e darmi il mapping con l'id dell'ontologia MONDO per ognuno di essi, che sia corretto, lo volgio in un json

mondo_mapping = {
  "Squamous Cell Neoplasms": "MONDO:0002532",
  "Adnexal and Skin Appendage Neoplasms": "MONDO:0002297",
  "Epithelial Neoplasms, NOS": "MONDO:0005626",
  "Complex Epithelial Neoplasms": "MONDO:0005626",
  "Fibroepithelial Neoplasms": "MONDO:0021045",
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
# from omics_connector import make_mirna_hgcn_map
# make_mirna_hgcn_map()
pkt_nodelabels[pkt_nodelabels['bioentity_type'] == "sequence"].class_code.value_counts()




# %%   ================== LOAD OMICS DATAFRAMES ==============================
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
# df= read_data_type(DATA_TYPES[5])
# print(df.columns)
# df

#%%  Load omics metadata json files
METADATA_PATHS
def load_omics_metadata(data_type):
    """Load metadata JSON for a specific omics data type."""
    file_path = METADATA_PATHS.get(data_type)
    if file_path and os.path.exists(file_path):
        with open(file_path, 'r', encoding='utf-8') as f:
            metadata = json.load(f)
        return metadata
    else:
        print(f"Metadata file not found for data type {data_type}: {file_path}")
        return None

load_omics_metadata("star_tpm")

#%% ------ LOAD all MAPPINGS
ROOT_MAPS

maps_file ={ 
  "bimart_genes_map":   "biomart_gene_mappings.tsv",
  "hm27_map":  "HM27.hg38.manifest.gencode.v36.probeMap",
  "hm450_map":  "HM450.hg38.manifest.gencode.v36.probeMap",
  "mirna_map":  "mirna_hgcn_map.tsv",
  "rppa_map":  "TCPA_RPPA500_metadata_mapping.tsv",
  "rsid_map":  "rsid_checkpoint.tsv", 
  #"rsid_hgcn":  "rsid_hgcn_mapping.tsv",
  }


# load all mapping files into a dictionary of dataframes with keys as map names
mapping_dfs = {}
for map_name, file_name in maps_file.items():
    file_path = f"{ROOT_MAPS}{file_name}"
    if os.path.exists(file_path):
        df_map = pd.read_csv(file_path, sep='\t', dtype=str)
        mapping_dfs[map_name] = df_map
        print(f"Loaded mapping: {map_name}, shape: {df_map.shape}")
    else:
        print(f"Mapping file not found: {file_path}")

mapping_dfs.keys()

#%% Load PANCAN data
pancan_phenotype_path = "../data/omics/PANCAN/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
pancan_phenotype_df = pd.read_csv(pancan_phenotype_path, sep='\t', compression='gzip', low_memory=False)
pancan_phenotype_df.columns

pancan_subtypes_path = "../data/omics/PANCAN/TCGASubtype.20170308.tsv.gz"
pancan_subtypes_df = pd.read_csv(pancan_subtypes_path, sep='\t', compression='gzip', low_memory=False)
pancan_subtypes_df.head()

#%%
 
mapping_dfs['bimart_genes_map'].head()
#%% --- Load to ARANGODB collections ---

# Ad ognuno anda aggiunta la chiave epr il mapping con il KNOWLEDGE GRAPH !

# 1. Collezione PROJECTS (una per progetto TCGA)
# data_source : "clinical"

clinical_df = read_data_type("clinical")
clinical_df.disease_type.unique()

#%% 1. Collezione PROJECT

PROJECT_example = {
  "_key": "TCGA-BRCA",
    "name": "Breast Invasive Carcinoma",
    "program": "TCGA",
    "primary_site": "Breast",
    "disease_type": [
      "Squamous Cell Neoplasms",
      "Adnexal and Skin Appendage Neoplasms",
      "Epithelial Neoplasms, NOS",
      "Complex Epithelial Neoplasms",
      "Fibroepithelial Neoplasms",
      "Cystic, Mucinous and Serous Neoplasms",
      "Basal Cell Neoplasms",
      "Adenomas and Adenocarcinomas",
      "Ductal and Lobular Neoplasms"
    ]
  }

PROJECT_keys = PROJECT_example.keys()

for key in PROJECT_keys:
    if key in clinical_df.columns:
        print(f"{key}: {clinical_df[key].dtype}, unique values: {clinical_df[key].nunique()}")
    else:
        print(f"{key}: NOT FOUND in clinical_df columns")

import ast
PROJECT_load = {
    "_key": clinical_df.loc[0, 'project_id.project'],
    "name": clinical_df.loc[0, 'name.project'],
    "program": clinical_df.loc[0, 'name.program.project'],
    "primary_site": clinical_df.loc[0, 'primary_site.project'],
    "disease_types": ast.literal_eval(clinical_df.loc[0, 'disease_type.project'])
}
PROJECT_load
#%% 2. Collezione SAMPLES

SAMPLES_example = {
  "_key": "TCGA-D8-A1XU-01A",
  "case_id": "TCGA-D8-A1XU",
  "sample_id": "9ad80473-9148-45a7-ad6d-9936a899e022",
  "sample_type": "Primary Tumor",
  "sample_type_id": 1,
  "tumor_descriptor": "Primary",
  "specimen_type": "Solid Tissue",
  "tissue_type": "Tumor",
  "composition": "Not Reported",
  "preservation_method": "OCT",
  "oct_embedded": True,
  "is_ffpe": False,
  "initial_weight": 120.0,
  "days_to_collection": 85,
  "days_to_sample_procurement": None,
  "pathology_report_uuid": "801A4E2F-E26E-424F-BF42-CD0D9CD62BCE",
  "annotations": None
}
SAMPLES_keys = SAMPLES_example.keys()
for key in SAMPLES_keys:
    """
    attenzione! clinical_df_columns sono innestate con . 
    esempio sample_type_id.samples
    """
    col_name = f"{key}.samples"
    if col_name in clinical_df.columns:
        print(f"{key}: {clinical_df[col_name].dtype}, unique values: {clinical_df[col_name].nunique()}")
    else:
        print(f"{key}: NOT FOUND in clinical_df columns")


SAMPLES_load = {
    "_key": clinical_df.loc[0, 'sample'],
    "submitter_id": clinical_df.loc[0, 'submitter_id'],
    "tcga_project": clinical_df.loc[0, 'project_id.project'],
    # "case_id": clinical_df.loc[0, 'case_id'],
    # "sample_id": clinical_df.loc[0, 'sample_id.samples'],
    "sample_type": clinical_df.loc[0, 'sample_type.samples'],
    "sample_type_id": clinical_df.loc[0, 'sample_type_id.samples'],
    # "disease_type": clinical_df.loc[0, 'disease_type'],
    "tumor_descriptor": clinical_df.loc[0, 'tumor_descriptor.samples'],
    "specimen_type": clinical_df.loc[0, 'specimen_type.samples'],
    "tissue_type": clinical_df.loc[0, 'tissue_type.samples'],
    "composition": clinical_df.loc[0, 'composition.samples'],
    "preservation_method": clinical_df.loc[0, 'preservation_method.samples'],
    "oct_embedded": clinical_df.loc[0, 'oct_embedded.samples'],
    "is_ffpe": clinical_df.loc[0, 'is_ffpe.samples'],
    "initial_weight": clinical_df.loc[0, 'initial_weight.samples'],
    "days_to_collection": clinical_df.loc[0, 'days_to_collection.samples'],
    # "days_to_sample_procurement": clinical_df.loc[0, 'days_to_sample_procurement.samples'],
    "pathology_report_uuid": clinical_df.loc[0, 'pathology_report_uuid.samples'],
    # "annotations": clinical_df.loc[0, 'annotations.samples'],
}
SAMPLES_load
#%% 3. Collezione CASES (pazienti/casi con dati clinici)
# data_source : "clinical"
clinical_df = read_data_type("clinical")
survival_df = read_data_type("survival")

CASES_example = {
  "_key": "TCGA-BH-A0W3",
  "id": "3c612e12-6de8-44fa-a095-805c45474821",
  "submitter_id": "TCGA-BH-A0W3",
  "tcga_project": "TCGA-BRCA",
  "primary_site": "Breast",
  "disease_type": "Ductal and Lobular Neoplasms",
  
  "demographic": {
    "gender": "female",
    "race": "white",
    "ethnicity": "not hispanic or latino",
    "vital_status": "Alive",
    "age_at_index": 58,
    "days_to_birth": -21369,
    "year_of_birth": 1952,
    "year_of_death": None,
    "days_to_death": None
  },
  
  "exposures": {
    "alcohol_history": "Not Reported"
  },
  
  "diagnoses": [
    {
      "primary_diagnosis": "Infiltrating duct carcinoma, NOS",
      "tissue_or_organ_of_origin": "Breast, NOS",
      "site_of_resection_or_biopsy": "Breast, NOS",
      "morphology": "8500/3",
      "icd_10_code": "C50.9",
      "age_at_diagnosis": 21369,
      "age_at_diagnosis_years": 58.55,
      "days_to_diagnosis": 0,
      "year_of_diagnosis": 2010,
      "tumor_grade": "Not Reported",
      "classification_of_tumor": "not reported",
      "prior_malignancy": "no",
      "synchronous_malignancy": "No",
      "prior_treatment": "No",
      "progression_or_recurrence": "not reported",
      "last_known_disease_status": "not reported",
      "days_to_last_follow_up": 728,
      "ajcc_pathologic_stage": "Stage IIA",
      "ajcc_staging_system_edition": "7th",
      "ajcc_pathologic_t": "T1c",
      "ajcc_pathologic_n": "N1a",
      "ajcc_pathologic_m": "M0"
    }
  ],
  
  "treatments": [
    {
      "treatment_id": "1c5f5df0-9317-51f8-95fb-88357c15289d",
      "submitter_id": "TCGA-BH-A0W3_treatment_1",
      "treatment_type": "Pharmaceutical Therapy, NOS",
      "treatment_or_therapy": "no",
      "created_datetime": "2019-04-28T13:48:22.214306-05:00",
      "updated_datetime": "2019-07-31T21:37:34.195388-05:00",
      "state": "released"
    },
    {
      "treatment_id": "9c5a28c6-258e-5f41-a6bf-4893cb797ac3",
      "submitter_id": "TCGA-BH-A0W3_treatment",
      "treatment_type": "Radiation Therapy, NOS",
      "treatment_or_therapy": "yes",
      "created_datetime": None,
      "updated_datetime": "2019-07-31T21:37:34.195388-05:00",
      "state": "released"
    }
  ],
  
  "tissue_source_site": {
    "id": "ad5db77f-ce9a-53c8-b7ff-7944acf5c0c6",
    "code": "BH",
    "name": "University of Pittsburgh",
    "project": "Breast invasive carcinoma",
    "bcr_id": "NCH"
  },
    
  
}
CASES_keys = CASES_example.keys()
# crea una versione delle keys nested con i punti come3 nel dataset clinical
CASES_keys_nestes_with_dot = []
for key in CASES_keys:
    if isinstance(CASES_example[key], dict):
        for subkey in CASES_example[key].keys():
            CASES_keys_nestes_with_dot.append(f"{subkey}.{key}")
    elif isinstance(CASES_example[key], list):
        # skip lists for now
        continue
    else:
        CASES_keys_nestes_with_dot.append(key)
CASES_keys_nestes_with_dot

for key in CASES_keys_nestes_with_dot:
    if key in clinical_df.columns:
        print(f"{key}: {clinical_df[key].dtype}, unique values: {clinical_df[key].nunique()}")
    else:
        print(f"{key}: NOT FOUND in clinical_df columns")

# Add to this collection : CLINICAL_SURVIVAL
# data_source : "survival"

survival_df.head().to_clipboard()
#%%

row = 0 # first row
sample_id = clinical_df.loc[row, 'sample']

CASES_load = {
    # "_key": clinical_df.loc[0, 'case_id'],
    "_key": clinical_df.loc[row, 'submitter_id'],
    # "id": clinical_df.loc[0, 'id'],
    # "submitter_id": clinical_df.loc[0, 'submitter_id'],
    "tcga_project": clinical_df.loc[row, 'project_id.project'],
    "primary_site": clinical_df.loc[row, 'primary_site'],
    "disease_type": clinical_df.loc[row, 'disease_type'],
    
    "demographic": {
        "gender": clinical_df.loc[row, 'gender.demographic'],
        "race": clinical_df.loc[row, 'race.demographic'],
        "ethnicity": clinical_df.loc[row, 'ethnicity.demographic'],
        "vital_status": clinical_df.loc[row, 'vital_status.demographic'],
        "age_at_index": clinical_df.loc[row, 'age_at_index.demographic'],
        "days_to_birth": clinical_df.loc[row, 'days_to_birth.demographic'],
        "year_of_birth": clinical_df.loc[row, 'year_of_birth.demographic'],
        "year_of_death": clinical_df.loc[row, 'year_of_death.demographic'],
        "days_to_death": clinical_df.loc[row, 'days_to_death.demographic']
    },

    "diagnoses": [{
        "primary_diagnosis": clinical_df.loc[row, 'primary_diagnosis.diagnoses'],
        "tissue_or_organ_of_origin": clinical_df.loc[row, 'tissue_or_organ_of_origin.diagnoses'],
        "site_of_resection_or_biopsy": clinical_df.loc[row, 'site_of_resection_or_biopsy.diagnoses'],
        "morphology": clinical_df.loc[row, 'morphology.diagnoses'],
        "icd_10_code": clinical_df.loc[row, 'icd_10_code.diagnoses'],
        "age_at_diagnosis": clinical_df.loc[row, 'age_at_diagnosis.diagnoses'],
        "age_at_diagnosis_years": clinical_df.loc[row, 'age_at_earliest_diagnosis_in_years.diagnoses.xena_derived'],
        "days_to_diagnosis": clinical_df.loc[row, 'days_to_diagnosis.diagnoses'],
        "year_of_diagnosis": clinical_df.loc[row, 'year_of_diagnosis.diagnoses'],
        "tumor_grade": clinical_df.loc[row, 'tumor_grade.diagnoses'],
        "classification_of_tumor": clinical_df.loc[row, 'classification_of_tumor.diagnoses'],
        "prior_malignancy": clinical_df.loc[row, 'prior_malignancy.diagnoses'],
        "synchronous_malignancy": clinical_df.loc[row, 'synchronous_malignancy.diagnoses'],
        "prior_treatment": clinical_df.loc[row, 'prior_treatment.diagnoses'],
        "progression_or_recurrence": clinical_df.loc[row, 'progression_or_recurrence.diagnoses'],
        "last_known_disease_status": clinical_df.loc[row, 'last_known_disease_status.diagnoses'],
        "days_to_last_follow_up": clinical_df.loc[row, 'days_to_last_follow_up.diagnoses'],
        "ajcc_pathologic_stage": clinical_df.loc[row, 'ajcc_pathologic_stage.diagnoses'],
        "ajcc_staging_system_edition": clinical_df.loc[row, 'ajcc_staging_system_edition.diagnoses'],
        "ajcc_pathologic_t": clinical_df.loc[row, 'ajcc_pathologic_t.diagnoses'],
        "ajcc_pathologic_n": clinical_df.loc[row, 'ajcc_pathologic_n.diagnoses'],
        "ajcc_pathologic_m": clinical_df.loc[row, 'ajcc_pathologic_m.diagnoses']
    }],
    
    "treatments": [{
        "treatment_id": clinical_df.loc[row, 'treatment_id.treatments.diagnoses'],
        "submitter_id": clinical_df.loc[row, 'submitter_id.treatments.diagnoses'],
        "treatment_type": clinical_df.loc[row, 'treatment_type.treatments.diagnoses'],
        "treatment_or_therapy": clinical_df.loc[row, 'treatment_or_therapy.treatments.diagnoses'],
        "created_datetime": clinical_df.loc[row, 'created_datetime.treatments.diagnoses'],
        "updated_datetime": clinical_df.loc[row, 'updated_datetime.treatments.diagnoses'],
        "state": clinical_df.loc[row, 'state.treatments.diagnoses']
    }],
    
    "exposures": {
        "alcohol_history": clinical_df.loc[row, 'alcohol_history.exposures']
    },
    
    "survival": {
        # use sample_id to link survival data
        "sample_id": sample_id,
        "OS_time": survival_df.loc[survival_df['sample'] == sample_id, 'OS.time'].values[0],
        "OS": survival_df.loc[survival_df['sample'] == sample_id, 'OS'].values[0]
    },

    # "tissue_source_site": {
    #     "id": clinical_df.loc[0, 'tissue_source_site_id.tissue_source_site'],
    #     "code": clinical_df.loc[0, 'code.tissue_source_site'],
    #     "name": clinical_df.loc[0, 'name.tissue_source_site'],
    #     "project": clinical_df.loc[0, 'project.tissue_source_site'],
    #     "bcr_id": clinical_df.loc[0, 'bcr_id.tissue_source_site']
    # }
}
CASES_load
#%%
survival_df.head()
#%% --------------------------------------------------------------------------
clinical_df.head()
#%% 4. Collezione GENE_EXPRESSION (normalizzato)
# data_source : ["star_counts", "star_fpkm", "star_tpm"]
GENE_EXPRESSION_example ={
  "_key": "expr_TCGA-D8-A1XU-01A_ENSG00000223972.5",
  "sample_id": "TCGA-D8-A1XU-01A",
  "gene_id": "ENSG00000223972.5",
  "star_counts": 11.737669863177452,
  "star_fpkm": 5.123456,
  "star_tpm": 6.543210,
  "platform": "STAR"
}
GENE_EXPRESSION_keys = GENE_EXPRESSION_example.keys()

star_counts_df = read_data_type("star_counts")
# star_counts_df.head()
star_fpkm_df = read_data_type("star_fpkm")
star_tpm_df = read_data_type("star_tpm")
#%% --------------------------------------------------------------------------
""" data structure
	Ensembl_ID	TCGA-D8-A146-01A	TCGA-AQ-A0Y5-01A	TCGA-C8-A274-01A	TCGA-BH-A0BD-01A
0	ENSG00000000003.15	11.737669863177452	9.78135971352466	13.122504484313904	11.016808287686557
1	ENSG00000000005.6	7.721099188707185	3.321928094887362	0.0	6.6865005271832185
2	ENSG00000000419.13	11.042343379793692	11.357552004618084	11.506307581100309	10.801708358916462

"""
#%%
GENE = star_tpm_df.loc[0, 'Ensembl_ID']
print(f"Processing gene: {GENE}")
# use mapping_dfs to find the node_id in PKT for this GENE
mapping_dfs["bimart_genes_map"].head()
# columns are gene_stable_id_version and entrez_id
# create mappings lookup dictionary
gene_to_entrez = pd.Series(mapping_dfs["bimart_genes_map"]['entrez_id'].values,index=mapping_dfs["bimart_genes_map"]['gene_stable_id']).to_dict()
node_id = str(gene_to_entrez.get(GENE.split(".")[0], "Unknown").split('.')[0])  # remove version if present 
print(f"Mapped to node_id (entrez_id): {node_id}")


GENE_EXPRESSION_load = {
        "_key": f"expr_{star_tpm_df.columns[1]}_{GENE}",
        "sample_id": star_tpm_df.columns[1],
        "gene_id": GENE,
        "node_id": node_id,  # KG mapping
        "star_counts": star_counts_df.loc[star_tpm_df['Ensembl_ID'] == GENE, star_tpm_df.columns[1]].values[0],
        "star_fpkm": star_fpkm_df.loc[star_tpm_df['Ensembl_ID'] == GENE, star_tpm_df.columns[1]].values[0],
        "star_tpm": star_tpm_df.loc[star_tpm_df['Ensembl_ID'] == GENE, star_tpm_df.columns[1]].values[0],
        "platform": "STAR"
    }
GENE_EXPRESSION_load
#%%
# query PKT to verify node_id exists
nodes_query = get_nodes_by_pattern(db_connection, "nodes", "_key", f"{node_id}")
for n in nodes_query:
    print(f"{n['class_code']}: {n['label']}")


#%%
GENES_EXP_LOAD = []

for GENE in star_tpm_df['Ensembl_ID'].head():
    node_id = str(gene_to_entrez.get(GENE.split(".")[0], "Unknown").split('.')[0])  # remove version if present 
    GENE_EXPRESSION_load = {
        "_key": f"expr_{star_tpm_df.columns[1]}_{GENE}",
        "sample_id": star_tpm_df.columns[1],
        "gene_id": GENE,
        "node_id": node_id,  # KG mapping
        "star_counts": star_counts_df.loc[star_tpm_df['Ensembl_ID'] == GENE, star_tpm_df.columns[1]].values[0],
        "star_fpkm": star_fpkm_df.loc[star_tpm_df['Ensembl_ID'] == GENE, star_tpm_df.columns[1]].values[0],
        "star_tpm": star_tpm_df.loc[star_tpm_df['Ensembl_ID'] == GENE, star_tpm_df.columns[1]].values[0],
        "platform": "STAR"
    }
    GENES_EXP_LOAD.append(GENE_EXPRESSION_load)

GENES_EXP_LOAD

#%% 5. Collezione CNV (Copy Number Variation - normalizzato)
# data_source : ["gene-level_ascat3", "allele_cnv_ascat3"]
{
  "_key": "cnv_TCGA-D8-A1XU-01A_chr1_17369_17436",
  "sample_id": "TCGA-D8-A1XU-01A",
  "chromosome": "chr1",
  "start": 17369,
  "end": 17436,
  "value": 2,
  "analysis_method": "allele_cnv_ascat3",
}
{"_key": "cnv_TCGA-D8-A1XU-01A_ENSG00000227232.5",
  "sample_id": "TCGA-D8-A1XU-01A",
  "value": 3,
  "analysis_method": "gene-level_ascat3",
}

cnv_gene_level_ascat3 = read_data_type("gene-level_ascat3")
cnv_gene_level_ascat3
#%%
for i in range(50):
    GENE = cnv_gene_level_ascat3.loc[i, 'Ensembl_ID']
    print(f"Processing gene: {GENE}")
    # use mapping_dfs to find the node_id in PKT for this GENE
    mapping_dfs["bimart_genes_map"].head()
    # columns are gene_stable_id_version and entrez_id
    # create mappings lookup dictionary
    gene_to_entrez = pd.Series(mapping_dfs["bimart_genes_map"]['entrez_id'].values,index=mapping_dfs["bimart_genes_map"]['gene_stable_id']).to_dict()
    node_id = str(str(gene_to_entrez.get(GENE.split(".")[0], "Unknown")).split('.')[0])  # remove version if present 
    print(f"Mapped to node_id (entrez_id): {node_id}")
#%%
mapping_dfs["bimart_genes_map"][["entrez_id", "hgnc_symbol"]].drop_duplicates()
#%% 6. Collezione MUTATIONS (normalizzato)
# data_source : "somaticmutation_wxs"


mapping_dfs['rsid_map']["gene_map"] = mapping_dfs['rsid_map']["38.gene"] +"_" + mapping_dfs['rsid_map']["38.pos"].astype(str)
mapping_dfs['rsid_map']["chr_map"] = mapping_dfs['rsid_map']["38.chr"] +"_" + mapping_dfs['rsid_map']["38.pos"].astype(str)
mapping_dfs['rsid_map']
#%%
somatic_df = read_data_type("somaticmutation_wxs")
somatic_df["gene_map"] = somatic_df["gene"]+"_" + somatic_df['start'].astype(str)
somatic_df["chr_map"] = somatic_df["chrom"].replace("chr", "") +"_" + somatic_df['start'].astype(str)
somatic_df[["gene", "start"]].drop_duplicates()
#%%
# merge somatic_df with mapping_dfs on chr_map
somatic_merged =  pd.merge(somatic_df, mapping_dfs['rsid_map'], on="gene_map", how="left")
somatic_merged[['gene', 'start', 'snp_id']]
somatic_merged["snp_id"].nunique()

# repeat the mapping using chromisome instes of gene ==============
# somatic_merged =  pd.merge(somatic_df, mapping_dfs['rsid_map'], on="chr_map", how="left")
# somatic_merged[['chrom', 'start', 'snp_id']]
# somatic_merged["snp_id"].nunique()



#%% 7. Collezione METHYLATION (normalizzato)
# data_source : ["methylation450", "methylation27"]


#%% 8. Collezione MIRNA (normalizzato)
# data_source : "mirna_seq"

#%% 9. Collezione PROTEIN (normalizzato)
# data_source : "protein"



#%% ----------------  [non farei queste collezioni ma mappere dentro le omics]  --------------------
# 11. Collezione GENE_MAP (reference biologico)
# data_source : "biomart_metadata"

#12 Collezione VARIANT_MAP (reference biologico)
# data_source : "rsid_metadata"