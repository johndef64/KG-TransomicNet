
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


#%% --------------------------------------------------------------------------
from omics_connector import read_data_type, PROBEMAP_PATHS

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
from pybiomart import Server

# Connessione al server Ensembl
server = Server(host='http://www.ensembl.org')
dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
print(f"Dataset: {dataset.display_name}")
print(f"Name: {dataset.name}")
# alt access
dataset = Dataset(name='hsapiens_gene_ensembl',
                  host='http://www.ensembl.org')

# Accedi al VARIATION MART (non al gene mart)
variation_mart = server.marts['ENSEMBL_MART_SNP']
somatic_dataset = variation_mart.datasets['hsapiens_snp_som']
somatic_attributes = somatic_dataset.list_attributes()

#%%

# get all attributes
attributes = dataset.attributes
dataset.list_attributes()
#%%
# remove attributes containig "homolog"
for key in list(attributes.keys()):
    if 'homolog' in key:
        del attributes[key]
    if 'paralog' in key:
        del attributes[key]
    if 'structure' in key:
        del attributes[key]
    if 'sequence' in key:
        del attributes[key]
    if 'affy' in key or 'array' in key or 'illumina' in key:
        del attributes[key]
print(attributes.keys())


ensembl_attributes = [
    # === IDENTIFICATORI GENE ===
    # Identificatori principali del gene Ensembl
    'ensembl_gene_id',
    'ensembl_gene_id_version',
    'external_gene_name',
    'external_gene_source',
    'gene_biotype',
    'version',
    'source',
    
    # === IDENTIFICATORI TRASCRITTO ===
    # Identificatori del trascritto e annotazioni di qualità
    'ensembl_transcript_id',
    'ensembl_transcript_id_version',
    'external_transcript_name',
    'external_transcript_source_name',
    'transcript_biotype',
    'transcript_source',
    'transcript_version',
    'transcript_count',
    'transcript_tsl',  # Transcript Support Level
    'transcript_gencode_basic',
    'transcript_primary_basic',
    'transcript_appris',  # Annotation of principal isoforms
    'transcript_is_canonical',
    'transcript_mane_select',  # MANE Select transcript
    'transcript_mane_plus_clinical',  # MANE Plus Clinical
    
    # === IDENTIFICATORI PEPTIDE/PROTEINA ===
    # Identificatori delle sequenze proteiche
    'ensembl_peptide_id',
    'ensembl_peptide_id_version',
    'peptide_version',
    'protein_id',
    
    # === IDENTIFICATORI ESONE ===
    # Identificatori degli esoni
    'ensembl_exon_id',
    'exon_id',
    'end_exon_id',
    'start_exon_id',
    
    # === DESCRIZIONE E ANNOTAZIONE ===
    # Descrizioni testuali e annotazioni funzionali
    'description',
    'external_synonym',
    
    # === COORDINATE GENOMICHE GENE ===
    # Posizione genomica del gene
    'chromosome_name',
    'start_position',
    'end_position',
    'strand',
    'band',  # Cytogenetic band
    
    # === COORDINATE TRASCRITTO ===
    # Posizione e lunghezza del trascritto
    'transcript_start',
    'transcript_end',
    'transcription_start_site',
    'transcript_length',
    
    # === COORDINATE ESONE ===
    # Coordinate degli esoni e informazioni strutturali
    'exon_chrom_start',
    'exon_chrom_end',
    'is_constitutive',
    'rank',
    'phase',
    'end_phase',
    
    # === COORDINATE REGIONI CODIFICANTI (CDS) ===
    # Coordinate delle regioni codificanti
    'cdna_coding_start',
    'cdna_coding_end',
    'genomic_coding_start',
    'genomic_coding_end',
    'cds_start',
    'cds_end',
    'cds_length',
    'cds_start_2076',
    'cds_end_2076',
    'somatic_cds_start_2076',
    'somatic_cds_end_2076',
    'coding_start_offset',
    'coding_end_offset',
    
    # === COORDINATE UTR ===
    # Coordinate delle regioni non tradotte (5' e 3' UTR)
    '5_utr_start',
    '5_utr_end',
    '3_utr_start',
    '3_utr_end',
    '5utr',
    '3utr',
    
    # === CARATTERISTICHE GENOMICHE ===
    # Contenuto GC e caratteristiche della sequenza
    'percentage_gene_gc_content',
    'struct_transcript_count',
    'codon_table_id',
    
    # === FENOTIPI E STUDI ===
    # Associazioni fenotipiche e studi
    'phenotype_description',
    'Source_name',
    'study_external_id',
    'strain_name',
    'strain_gender',
    'p_value',
    
    # === GENE ONTOLOGY (GO) ===
    # Annotazioni Gene Ontology
    'go_id',
    'name_1006',
    'definition_1006',
    'go_linkage_type',
    'namespace_1003',
    'goslim_goa_accession',
    'goslim_goa_description',
    'go',
    
    # === DATABASE ESTERNI - GENE ===
    # Cross-references a database esterni per geni
    'biogrid',  # Protein-protein interactions
    'ccds',  # Consensus CDS
    'chembl',  # Bioactive molecules
    'genecards',
    'hgnc_symbol',  # HUGO Gene Nomenclature
    'hgnc_id',
    'hgnc_trans_name',
    'entrezgene_id',
    'entrezgene_accession',
    'entrezgene_description',
    'entrezgene_trans_name',
    'reactome',  # Pathway database
    'reactome_gene',
    'reactome_transcript',
    'wikigene_id',
    'wikigene_name',
    'wikigene_description',
    
    # === DATABASE ESTERNI - SEQUENZE ===
    # Cross-references per sequenze nucleotidiche e proteiche
    'embl',  # European nucleotide archive
    'refseq_mrna',
    'refseq_mrna_predicted',
    'refseq_ncrna',
    'refseq_ncrna_predicted',
    'refseq_peptide',
    'refseq_peptide_predicted',
    'rnacentral',  # Non-coding RNA
    'rfam',  # RNA families
    'rfam_trans_name',
    'ens_lrg_gene',  # Locus Reference Genomic
    'ens_lrg_transcript',
    
    # === DATABASE ESTERNI - PROTEINE ===
    # Cross-references per proteine e domini
    'uniparc',  # UniProt Archive
    'uniprot_gn_id',
    'uniprot_gn_symbol',
    'uniprot_isoform',
    'uniprotswissprot',
    'uniprotsptrembl',
    'pdb',  # Protein structure
    'string',  # Protein-protein interaction network
    'hpa_id',  # Human Protein Atlas
    'hpa_accession',
    'merops',  # Peptidase database
    
    # === DATABASE ESTERNI - MALATTIE ===
    # OMIM (Online Mendelian Inheritance in Man)
    'mim_gene_description',
    'mim_gene_accession',
    'mim_morbid_description',
    'mim_morbid_accession',
    
    # === DATABASE ESTERNI - miRNA ===
    # microRNA databases
    'mirbase_id',
    'mirbase_accession',
    'mirbase_trans_name',
    
    # === DATABASE ESTERNI - SPLICING ===
    # Alternative splicing databases
    'dbass3_id',  # Database of aberrant 3' splice sites
    'dbass3_name',
    'dbass5_id',  # Database of aberrant 5' splice sites
    'dbass5_name',
    
    # === DATABASE ESTERNI - ALTRI ===
    # Altri cross-references
    'arrayexpress',
    'ucsc',  # UCSC Genome Browser
    
    # === ARRAY E MICROARRAY ===
    # Identificatori di sonde su piattaforme microarray
    # Affymetrix arrays
    # Agilent arrays
    # Altri arrays
    
    # === DOMINI PROTEICI - DATABASE SPECIFICI ===
    # Database di domini e famiglie proteiche con coordinate
    # CDD - Conserved Domain Database
    'cdd',
    'cdd_start',
    'cdd_end',
    # Gene3D
    'gene3d',
    'gene3d_start',
    'gene3d_end',
    # HAMAP
    'hamap',
    'hamap_start',
    'hamap_end',
    # NCBI Protein Family Models
    'ncbifam',
    'ncbifam_start',
    'ncbifam_end',
    # PANTHER
    'hmmpanther',
    'hmmpanther_start',
    'hmmpanther_end',
    # Pfam
    'pfam',
    'pfam_start',
    'pfam_end',
    # PIRSF
    'pirsf',
    'pirsf_start',
    'pirsf_end',
    # PRINTS
    'prints',
    'prints_start',
    'prints_end',
    # PROSITE
    'scanprosite',
    'scanprosite_start',
    'scanprosite_end',
    'pfscan',
    'pfscan_start',
    'pfscan_end',
    # SFLD
    'sfld',
    'sfld_start',
    'sfld_end',
    # SMART
    'smart',
    'smart_start',
    'smart_end',
    # SUPERFAMILY
    'superfamily',
    'superfamily_start',
    'superfamily_end',
    # TIGRFAMs
    'tigrfam',
    'tigrfam_start',
    'tigrfam_end',
    
    # === DOMINI PROTEICI - INTERPRO ===
    # InterPro unified protein domain database
    'interpro',
    'interpro_short_description',
    'interpro_description',
    'interpro_start',
    'interpro_end',
    
    # === PREDIZIONI STRUTTURALI PROTEICHE ===
    # Predizioni di struttura e caratteristiche proteiche
    # AlphaFold
    'alphafold',
    'alphafold_start',
    'alphafold_end',
    # MobiDB-lite (disordini intrinseci)
    'mobidblite',
    'mobidblite_start',
    'mobidblite_end',
    # Ncoils (coiled-coils)
    'ncoils',
    'ncoils_start',
    'ncoils_end',
    # Phobius (transmembrane topology)
    'phobius',
    'phobius_start',
    'phobius_end',
    # SEG (low-complexity regions)
    'seg',
    'seg_start',
    'seg_end',
    # SIFTS (Structure Integration with Function, Taxonomy and Sequences)
    'sifts_import',
    'sifts_import_start',
    'sifts_import_end',
    # SignalP (signal peptides)
    'signalp',
    'signalp_start',
    'signalp_end',
    # Gram-negative signal peptides
    'signal_gn',
    'signal_gn_start',
    'signal_gn_end',
    # Gram-positive signal peptides
    'signal_gp',
    'signal_gp_start',
    'signal_gp_end',
    # TMHMM (transmembrane helices)
    'tmhmm',
    'tmhmm_start',
    'tmhmm_end',

    # === VARIANTI GERMINALI (SNP) - IDENTIFICATORI GENE/TRASCRITTO ===
    # Identificatori associati a SNP germinali
    'snp_ensembl_gene_id',
    'snp_gene_stable_id_version',
    'snp_gene_version',
    'snp_ensembl_transcript_id',
    'snp_transcript_stable_id_version',
    'snp_transcript_version',
    'snp_ensembl_peptide_id',
    'snp_translation_stable_id_version',
    'snp_peptide_version',
    
    # === VARIANTI GERMINALI (SNP) - COORDINATE GENOMICHE ===
    # Coordinate genomiche per SNP germinali
    'snp_chromosome_name',
    'snp_start_position',
    'snp_end_position',
    'snp_strand',
    'snp_band',
    'snp_external_gene_name',
    'snp_external_gene_source',
    
    # === VARIANTI GERMINALI (SNP) - CARATTERISTICHE ===
    # Lunghezze e caratteristiche delle sequenze con SNP
    'snp_ensembl_CDS_length',
    'snp_ensembl_cDNA_length',
    'snp_ensembl_peptide_length',
    'snp_transcript_count',
    'snp_percentage_gc_content',
    'snp_description',
    
    # === VARIANTI GERMINALI (SNP) - INFORMAZIONI VARIANTE ===
    # Dettagli della variante germinale
    'variation_name',
    'germ_line_variation_source',
    'source_description',
    'allele',
    'validated',
    'mapweight',
    'minor_allele',
    'minor_allele_freq',
    'minor_allele_count',
    'clinical_significance',
    
    # === VARIANTI GERMINALI (SNP) - POSIZIONI E CONSEGUENZE ===
    # Posizioni della variante e predizioni funzionali
    'transcript_location',
    'snp_chromosome_strand',
    'peptide_location',
    'chromosome_start',
    'chromosome_end',
    'polyphen_prediction_2076',  # PolyPhen-2 prediction
    'polyphen_score_2076',
    'sift_prediction_2076',  # SIFT prediction
    'sift_score_2076',
    'distance_to_transcript_2076',
    'peptide_shift',
    'synonymous_status',
    'allele_string_2076',
    
    # === VARIANTI SOMATICHE - IDENTIFICATORI GENE/TRASCRITTO ===
    # Identificatori associati a varianti somatiche
    'snp_som_ensembl_gene_id',
    'snp_som_gene_stable_id_version',
    'snp_som_gene_version',
    'snp_som_ensembl_transcript_id',
    'snp_som_transcript_stable_id_version',
    'snp_som_transcript_version',
    'snp_som_ensembl_peptide_id',
    'snp_som_translation_stable_id_version',
    'snp_som_peptide_version',
    
    # === VARIANTI SOMATICHE - COORDINATE GENOMICHE ===
    # Coordinate genomiche per varianti somatiche
    'snp_som_chromosome_name',
    'snp_som_start_position',
    'snp_som_end_position',
    'snp_som_strand',
    'snp_som_band',
    'snp_som_external_gene_name',
    'snp_som_external_gene_source',
    
    # === VARIANTI SOMATICHE - CARATTERISTICHE ===
    # Lunghezze e caratteristiche per varianti somatiche
    'snp_som_ensembl_CDS_length',
    'snp_som_ensembl_cDNA_length',
    'snp_som_ensembl_peptide_length',
    'snp_som_transcript_count',
    'snp_som_percentage_gc_content',
    'snp_som_description',
    
    # === VARIANTI SOMATICHE - INFORMAZIONI VARIANTE ===
    # Dettagli della variante somatica
    'somatic_variation_name',
    'somatic_source_name',
    'somatic_source_description',
    'somatic_allele',
    'somatic_validated',
    'somatic_mapweight',
    
    # === VARIANTI SOMATICHE - POSIZIONI E CONSEGUENZE ===
    # Posizioni e conseguenze delle varianti somatiche
    'somatic_transcript_location',
    'somatic_snp_chromosome_strand',
    'somatic_peptide_location',
    'somatic_chromosome_start',
    'somatic_chromosome_end',
    'mart_transcript_variation_som__dm_distance_to_transcript_2076',
    'somatic_synonymous_status',
    'mart_transcript_variation_som__dm_allele_string_2076',
    
    # === MODIFICHE POST-TRASCRIZIONALI ===
    # RNA editing e modifiche delle sequenze
    'seq_edits',
    'rna_seq_edits',
    
    # === SEQUENZE E REGIONI ===
    # Chiavi per recuperare vari tipi di sequenze
    'transcript_id_key',
    'transcript_exon_intron',  # Sequenza esoni+introni del trascritto
    'gene_exon_intron',  # Sequenza esoni+introni del gene
    'transcript_flank',  # Sequenze fiancheggianti il trascritto
    'gene_flank',  # Sequenze fiancheggianti il gene
    'coding_transcript_flank',  # Sequenze fiancheggianti la regione codificante del trascritto
    'coding_gene_flank',  # Sequenze fiancheggianti la regione codificante del gene
    'gene_exon',  # Sequenza esonicadel gene
    'cdna',  # Sequenza cDNA
    'coding',  # Sequenza codificante (CDS)
    'peptide',  # Sequenza peptidica
    'upstream_flank',  # Sequenza upstream
    'downstream_flank',  # Sequenza downstream
]


#%%
somatic_dataset.attributes.keys()   
somatic_dataset.filters.keys()

get_somatic= somatic_dataset.query(attributes=[
    'refsnp_id', 'refsnp_source', 'refsnp_source_description', 'chr_name'
    ],
    filters={'chr_name': ['21']}
    )
get_somatic["Variant source"].value_counts()
#%%


gene_mappings = dataset.query(attributes=[
    'ensembl_gene_id', 
    'ensembl_gene_id_version', 

    'entrezgene_id', 
    'hgnc_symbol',
    # 'ensembl_transcript_id',
    'gene_biotype',

    # 'external_gene_name',
    # 'ensembl_peptide_id',
    'description',
    'external_synonym',

    # 'reactome',  
    # 'reactome_gene',
	# 'mirbase_id',
    # 'mirbase_accession',
	
	'uniprot_gn_id',
    # 'uniprot_gn_symbol',
    # 'uniprot_isoform',
    # 'uniprotswissprot',
    # 'uniprotsptrembl',

    # # Posizione genomica del gene
    'chromosome_name',
    'start_position',
    'end_position',
    'strand',

    ],
              
    filters={'chromosome_name': ['21']}
    )


gene_mappings


#%%

# il tuo compito è comnporre il dataset di di gene mappings con le annotazioni di biomart
# ricorda che si possono usare solo 3 annotazioni esternae alla volta per query
# quindi devi fare piu query e poi unire i risultati, le annotazioni esterne sono quelle commentate sopra piu 
# 'entrezgene_id', 
# 'hgnc_symbol', 
# 'uniprot_gn_id',
# 'reactome',
# 'reactome_gene',
# 'mirbase_id',
# 'mirbase_accession', etc


#%%
# Compose comprehensive gene mappings with BioMart annotations
# Split external annotations into batches of 3 and merge results

# Base attributes (always included)
base_attributes = [
    'ensembl_gene_id', 
    'ensembl_gene_id_version',
    'gene_biotype',
    'description',
    # 'external_synonym',
    'chromosome_name',
    'start_position',
    'end_position',
    'strand',
]

# External annotations split into batches of 3
external_batches = [
    ['entrezgene_id', 'hgnc_symbol', 'mirbase_id'],
    # ['mirbase_accession']
    # ['reactome', 'reactome_gene'],
    # [ 'uniprot_gn_id'], # 'uniprot_gn_symbol'
    ['uniprotswissprot', 'uniprot_isoform' ] # "uniprotsptrembl"
]
def get_biomart_annotations(chromosome):
    # Initialize with base query
    print("Querying base attributes...")
    gene_mappings_full = dataset.query(
        attributes=base_attributes,
        filters={'chromosome_name': [chromosome]}
    )

    # Query each batch and merge
    for i, batch in enumerate(external_batches, 1):
        print(f"Querying batch {i}/{len(external_batches)}: {batch}")
        
        # Query with base ID + current batch
        temp_df = dataset.query(
            attributes=['ensembl_gene_id'] + batch,
            filters={'chromosome_name': [chromosome]}
        )
        
        # Merge with main dataframe
        gene_mappings_full = gene_mappings_full.merge(
            temp_df,
            on='Gene stable ID',
            how='left'
        )
        
        print(f"  Current shape: {gene_mappings_full.shape}")

    print(f"\nFinal dataset shape: {gene_mappings_full.shape}")
    print(f"Columns: {list(gene_mappings_full.columns)}")
    return gene_mappings_full
#%%
# Display result
chromosomes =["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"]
import time
for chrom in chromosomes:
    print(f"\nProcessing chromosome: {chrom}")
    gene_mappings = get_biomart_annotations(chrom)
    #save to file in temp folder
    gene_mappings.to_csv(f"../data/omics/maps/temp/biomart_gene_mappings_chr{chrom}.tsv", sep='\t', index=False)
    # merge or append to a master dataframe if needed
    if chrom == chromosomes[0]:
        master_gene_mappings = gene_mappings
    else:
        master_gene_mappings = pd.concat([master_gene_mappings, gene_mappings], ignore_index=True)

    time.sleep(2)  # Pause to avoid overwhelming the server

columns_dict = {
    'Gene stable ID': 'gene_stable_id',
    'Gene stable ID version': 'gene_stable_id_version',
    'Gene type': 'gene_type',
    'Gene description': 'gene_description',
    'Chromosome/scaffold name': 'chromosome',
    'Gene start (bp)': 'gene_start_bp',
    'Gene end (bp)': 'gene_end_bp',
    'Strand': 'strand',
    'NCBI gene (formerly Entrezgene) ID': 'entrez_id',
    'HGNC symbol': 'hgnc_symbol',
    'miRBase ID': 'mirbase_id',
    'UniProtKB/Swiss-Prot ID': 'uniprot_swissprot_id',
    'UniProtKB isoform ID': 'uniprot_isoform_id',

    'Reactome ID': 'reactome_id', 
    'Reactome gene ID': 'reactome_gene_id'
}
# replace column name and save as new file
master_gene_mappings.rename(columns=columns_dict, inplace=True)
master_gene_mappings.to_csv("../data/omics/maps/biomart_gene_mappings.tsv", sep='\t', index=False)
#%%
chr2 = get_biomart_annotations("2")
chr2.to_csv("../data/omics/maps/biomart_gene_mappings_chr2_alluniprot.tsv", sep='\t', index=False)
#%%
# master_gene_mappings.to_csv("../data/omics/maps/biomart_gene_mappings.tsv", sep='\t', index=False)    

#%%
# show a table  showing counts of uniques values on proprtey columns ["Gene type"] master_gene_mappings 
for col in ["gene_type", "hgnc_symbol", "uniprot_swissprot_id"]:
    print(f"\n=== Unique value counts for: {col} ===")
    print(f"Total unique values for: {col} = {master_gene_mappings[col].nunique()}")
    print(master_gene_mappings[col].value_counts())

#%%

# base_attributes = [
#     'ensembl_gene_id', 
#     'ensembl_gene_id_version',
#     'gene_biotype',
#     'description',
#     # 'external_synonym',
#     'chromosome_name',
#     'start_position',
#     'end_position',
#     'strand',
# ]

# gene_mappings_3 = dataset.query(
#     attributes=base_attributes,
#     filters={'chromosome_name': ["3"]}
# )
# gene_mappings_3

#%%
#%%
# Summary statistics
print("\n=== SUMMARY STATISTICS ===")
print(f"Total genes: {master_gene_mappings['gene_stable_id'].nunique()}")
print(f"\nGene biotype distribution:")
print(master_gene_mappings['gene_type'].value_counts())

print(f"\nExternal IDs coverage:")
for col in ['entrez_id', 'hgnc_symbol', 'uniprot_swissprot_id', 'reactome', 'mirbase_id']:
    if col in master_gene_mappings.columns:
        non_null = master_gene_mappings[col].notna().sum()
        print(f"  {col}: {non_null}/{len(master_gene_mappings)} ({non_null/len(gene_mappings_full)*100:.1f}%)")



#%%
