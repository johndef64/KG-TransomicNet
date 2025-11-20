#%%
import os 
import pandas as pd

OTHER_STUDIES = [
    "CGCI-BLGSP", "CGCI-HTMCP-CC", "CGCI-HTMCP-DLBCL", "CGCI-HTMCP-LC",
    "CMI-ASC", "CMI-MBC", "CMI-MPC",     
    "CPTAC-2", "CPTAC-3",
    "GDC_PANCAN",
    "APOLLO-LUAD",
    "BEATAML1.0-COHORT",
    "MMRF-COMMPASS",
    "HCMI-CMDC",
    "ORGANOID-PANCREATIC",
    "MP2PRT-ALL", "MP2PRT-WT",
    "REBC-THYR",
    "CTSP-DLBCL1",
    "NCICCR-DLBCL",
    "OHSU-CNL",
    "WCDT-MCRPC",
    "EXCEPTIONAL_RESPONDERS-ER",
    "CDDP_EAGLE-1",
]

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

GDC_ROOT_URL = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"

PROBEMAP_URLS = [
    # gene code annotation files
    f"{GDC_ROOT_URL}gencode.v36.annotation.gtf.gene.probemap",
    # methylation probe map files
    f"{GDC_ROOT_URL}HM450.hg38.manifest.gencode.v36.probeMap",
    f"{GDC_ROOT_URL}HM27.hg38.manifest.gencode.v36.probeMap",
]

# study = "TCGA-BRCA"
# data_type = "gene-level_ascat3"
# STUDY_URL = f"{GDC_ROOT_URL}{study}.{data_type}.tsv.gz"
# METADATA_URL = f"{GDC_ROOT_URL}{study}.{data_type}.tsv.json"


# get all data_types for the TCGA-n study 
def get_tcga_data_types(study_name):
    """
    Download all available data types for a TCGA study.
    
    Args:
        study_name: Nome dello studio (es. "TCGA-BRCA")
    
    Returns:
        tuple: (downloaded_files, missing_files)
    """
    import requests
    from pathlib import Path
    
    # Crea il percorso della cartella
    output_dir = Path("../data/omics") / study_name
    print(f"Output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Downloading data for {study_name}")
    print(f"{'='*60}\n")
    
    downloaded_files = []
    missing_files = []
    
    for data_type in DATA_TYPES:
        # Costruisci gli URL
        data_url = f"{GDC_ROOT_URL}{study_name}.{data_type}.tsv.gz"
        metadata_url = f"{GDC_ROOT_URL}{study_name}.{data_type}.tsv.json"
        
        # Nome dei file di output
        data_file = output_dir / f"{study_name}.{data_type}.tsv.gz"
        metadata_file = output_dir / f"{study_name}.{data_type}.tsv.json"
        
        # Verifica se il file esiste sul server
        try:
            response = requests.head(data_url, timeout=10)
            
            if response.status_code == 200:
                print(f"✓ Downloading {data_type}...")
                
                # Download del file dati
                data_response = requests.get(data_url, stream=True)
                with open(data_file, 'wb') as f:
                    for chunk in data_response.iter_content(chunk_size=8192):
                        f.write(chunk)
                
                # Download del metadata
                try:
                    metadata_response = requests.get(metadata_url)
                    if metadata_response.status_code == 200:
                        with open(metadata_file, 'wb') as f:
                            f.write(metadata_response.content)
                except:
                    pass
                
                downloaded_files.append(data_type)
                print(f"  → Saved to {data_file.name}")
                
            else:
                missing_files.append(data_type)
                print(f"✗ {data_type} - Not available (HTTP {response.status_code})")
                
        except requests.exceptions.RequestException as e:
            missing_files.append(data_type)
            print(f"✗ {data_type} - Error: {str(e)}")
    
    # Stampa il riepilogo
    print(f"\n{'='*60}")
    print(f"Download Summary for {study_name}")
    print(f"{'='*60}")
    print(f"Downloaded files: {len(downloaded_files)}/{len(DATA_TYPES)}")
    
    if downloaded_files:
        print(f"\n✓ Successfully downloaded:")
        for dt in downloaded_files:
            print(f"  - {dt}")
    
    if missing_files:
        print(f"\n✗ Missing/Failed:")
        for dt in missing_files:
            print(f"  - {dt}")

def get_all_probe_map(output_dir):
    """
    Download all probe map files.
    
    Args:
        output_dir: Directory di output per salvare i file
    """
    import requests
    from pathlib import Path
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Downloading Probe Maps")
    print(f"{'='*60}\n")
    
    for url in PROBEMAP_URLS:
        # Estrai il nome del file dall'URL
        file_name = url.split("/")[-1]
        probe_map_file = output_dir / file_name
        
        try:
            print(f"✓ Downloading {file_name}...")
            response = requests.get(url, stream=True)
            
            if response.status_code == 200:
                with open(probe_map_file, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"  → Saved to {probe_map_file.name}")
            else:
                print(f"✗ Error: HTTP {response.status_code}")
                
        except requests.exceptions.RequestException as e:
            print(f"✗ Error downloading {file_name}: {str(e)}")
    
    print(f"\n{'='*60}\n")

#### OMICS DATA Manipulation Functions

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
def make_mirna_hgcn_map():
    mirna_df = read_data_type("mirna")
    """Make maps for mirnaid to HGCNc gene symbol"""
    mirna_df["hgcn_id"] = mirna_df["miRNA_ID"].str.replace("hsa-mir-", "MIR").str.replace("hsa-let-", "MIRLET")
    # make al letter uppercase
    mirna_df["hgcn_id"] = mirna_df["hgcn_id"].str.upper()
    maps = mirna_df[["hgcn_id", "miRNA_ID"]]
    maps.to_csv(f"{ROOT_MAPS}mirna_hgcn_map.tsv", sep='\t', index=False)
    return maps

#%%

# Esempio di utilizzo
if __name__ == "__main__":
    # Download per un singolo studio
    get_tcga_data_types("TCGA-BRCA")
    
    # Oppure per tutti gli studi TCGA
    # for study in TCGA_STUDIES:
    #     get_tcga_data_types(study)

    get_all_probe_map("../data/omics/probe_maps")
#%%



#%%
"""
downloadhttps://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.somaticmutation_wxs.tsv.gz; Full metadata
samples986
version08-01-2024
type of datasomatic mutation (SNPs and small INDELs)
assemblyhg38
platformIllumina
authorGenomic Data Commons
wranglingHugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, HGVSp_Short and Consequence data are renamed accordingly and presented; dna_vaf data is added and is calculated by t_alt_count / t_depth.
input data formatVariant by Position (i.e. mutationVector)

alt	amino-acid	chrom	chromstart	chromend	effect	ref	sampleID
TCGA-AC-A2FK-01A				-1	-1			TCGA-AC-A2FK-01A
TCGA-EW-A2FV-01A	-	p.E306Kfs*36	chr1	953261	953261	frameshift_variant	C	TCGA-EW-A2FV-01A
TCGA-GM-A2DM-01A	G	p.W157R	chr1	956911	956911	missense_variant	A	TCGA-GM-A2DM-01A
TCGA-GM-A2DM-01A	A	p.R156S	chr1	956912	956912	missense_variant	T	TCGA-GM-A2DM-01A
TCGA-D8-A1XM-01A	A	p.V133M	chr1	961658	961658	missense_variant	G	TCGA-D8-A1XM-01A



miRNA Expression Quantification

TCGA-D8-A146-01A	TCGA-AQ-A0Y5-01A	TCGA-C8-A1HI-01A	TCGA-A7-A0CD-01A	TCGA-5L-AAT1-01A	TCGA-BH-A0C0-01A	TCGA-B6-A1KC-01B	TCGA-PL-A8LV-01A	TCGA-GM-A2DO-01A	TCGA-AN-A0XW-01A
hsa-let-7a-1	12.82	12.97	13.26	12.87	13.87	12.55	13.85	13.06	14.39	12.10
hsa-let-7a-2	12.82	12.97	13.25	12.87	13.86	12.52	13.86	13.05	14.39	12.13
hsa-let-7a-3	12.86	12.97	13.27	12.90	13.87	12.55	13.87	13.07	14.40	12.14
hsa-let-7b	14.73	14.50	13.84	14.94	14.61	14.48	15.29	13.07	14.00	14.58
hsa-let-7c	11.17	11.87	
"""
