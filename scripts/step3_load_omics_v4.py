#!/usr/bin/env python3
"""
step3_load_omics_v4.py - Sample-Centric Vector-Based TCGA Data Loading for ArangoDB

Implementazione basata su architettura IRR (Information Representation and Reuse) v4:
- SAMPLE-CENTRIC STRATEGY: Un documento per campione (non per gene×sample)
- VECTOR-BASED STORAGE: Vettori di espressione paralleli con indice gene condiviso
- DRASTICA RIDUZIONE DOCUMENTI: da ~20M a ~1K documenti per coorte
- OTTIMIZZATO PER: estrazione matrici, ML/GNN, analisi cohort-level

Architettura:
1. Semantic Layer: GENES (nodi stabili), PROJECTS, CASES, SAMPLES
2. Quantitative Index Layer: EXPRESSION_INDEX, CNV_INDEX (mapping gene→posizione)
3. Quantitative Vector Layer: GENE_EXPRESSION_SAMPLES (vettori per sample)

Versione: 4.0
Data: 2025-12-10
Autore: KG-Transomics Framework
Riferimento: Strategie Sample-Centric per dati trans-omici scalabili
"""

from sqlalchemy import null
from arangodb_utils import *
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
import os
from typing import Dict, List, Optional, Tuple
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ========================================================================
# CONFIGURATION
# ========================================================================

TCGA_STUDIES = [
    "TCGA-LAML", "TCGA-ACC", "TCGA-CHOL", "TCGA-BLCA", "TCGA-BRCA",
    "TCGA-CESC", "TCGA-COAD", "TCGA-UCEC", "TCGA-ESCA", "TCGA-GBM",
    "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-DLBC",
    "TCGA-LIHC", "TCGA-LGG", "TCGA-LUAD", "TCGA-LUSC", "TCGA-SKCM",
    "TCGA-MESO", "TCGA-UVM", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG",
    "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-STAD", "TCGA-TGCT",
    "TCGA-THYM", "TCGA-THCA", "TCGA-UCS"
]

DATA_TYPES = [
    "gene-level_ascat3",    # Copy Number (Gene Level)
    "allele_cnv_ascat3",    # Allele-specific Copy Number Segment
    "methylation27",         # DNA Methylation - Illumina Human Methylation
    "somaticmutation_wxs",  # Somatic Mutation
    "mirna",                 # miRNA Expression
    "star_counts",           # Gene Expression (STAR - counts)
    "star_fpkm",             # Gene Expression (STAR - FPKM)
    "star_tpm",              # Gene Expression (STAR - TPM)
    "protein",               # Protein Expression (RPPA)
    "clinical",              # Clinical Data
    "survival"               # Survival Data
]

# Paths configuration
STUDY = "TCGA-BRCA"  # Change as needed
ROOT_PATH = f"../data/omics/{STUDY}/"
ROOT_MAPS = f"../data/mappings/"
OUTPUT_PATH = f"../data/arangodb_collections/{STUDY}/"

# Ensure output directory exists
os.makedirs(OUTPUT_PATH, exist_ok=True)

# ========================================================================
# DATA LOADING UTILITIES
# ========================================================================

def load_omics_data(data_type: str) -> Optional[pd.DataFrame]:
    """Load TCGA omics data file for a specific data type."""
    file_path = f"{ROOT_PATH}{STUDY}.{data_type}.tsv.gz"
    if os.path.exists(file_path):
        logger.info(f"Loading {data_type} from {file_path}")
        return pd.read_csv(file_path, sep='\t', compression='gzip', low_memory=False)
    else:
        logger.warning(f"File not found: {file_path}")
        return None

def load_mappings() -> Dict[str, pd.DataFrame]:
    """Load all ID mapping files."""
    maps_files = {
        "biomart_genes": "biomart_gene_mappings_filled_extended.tsv",
        "hm27": "HM27.hg38.manifest.gencode.v36.probeMap",
        "hm450": "HM450.hg38.manifest.gencode.v36.probeMap",
        "mirna": "mirna_hgcn_map.tsv",
        "rppa": "TCPA_RPPA500_metadata_mapping.tsv",
        "rsid": "rsid_checkpoint.tsv",
        "mondo": "mondo_to_primarydisease.tsv",
        "ensg_enst": "biomart_ensg_enst_mapping.tsv"
    }

    mapping_dfs = {}
    for map_name, file_name in maps_files.items():
        file_path = f"{ROOT_MAPS}{file_name}"
        if os.path.exists(file_path):
            mapping_dfs[map_name] = pd.read_csv(file_path, sep='\t', dtype=str)
            logger.info(f"Loaded mapping: {map_name}, shape: {mapping_dfs[map_name].shape}")
        else:
            logger.warning(f"Mapping file not found: {file_path}")

    return mapping_dfs

def create_lookup_dicts(mapping_dfs: Dict[str, pd.DataFrame]) -> Dict[str, Dict]:
    """Create lookup dictionaries for fast ID mapping."""
    logger.info("Creating lookup dictionaries...")

    biomart = mapping_dfs["biomart_genes"]

    lookups = {
        "ensembl_to_entrez": pd.Series(biomart['entrez_id'].values, 
                                       index=biomart['gene_stable_id']).to_dict(),
        "ensembl_to_hgnc": pd.Series(biomart['hgnc_symbol'].values, 
                                     index=biomart['gene_stable_id']).to_dict(),
        "ensembl_to_uniprot": pd.Series(biomart['uniprot_swissprot_id'].values, 
                                        index=biomart['gene_stable_id']).to_dict(),
        "ensembl_to_mirbase": pd.Series(biomart['mirbase_id'].values, 
                                        index=biomart['gene_stable_id']).to_dict(),
    }

    # ENSG to ENST mapping (one-to-many)
    if "ensg_enst" in mapping_dfs:
        lookups["ensg_to_ensts"] = mapping_dfs["ensg_enst"].groupby('gene_id')['transcript_id'].apply(list).to_dict()

    logger.info(f"Created {len(lookups)} lookup dictionaries")
    return lookups

# ========================================================================
# COLLECTION BUILDERS - IRR v4 ARCHITECTURE
# ========================================================================

class CollectionBuilder:
    """Base class for building ArangoDB collections following IRR v4 principles."""

    def __init__(self, study: str, mapping_dfs: Dict, lookups: Dict):
        self.study = study
        self.mapping_dfs = mapping_dfs
        self.lookups = lookups

    def save_collection(self, documents: List[Dict], collection_name: str):
        """Save documents to JSON file for ArangoDB import."""
        output_file = f"{OUTPUT_PATH}{collection_name}.json"
        with open(output_file, 'w', encoding='utf-8') as f:
            for doc in documents:
                json.dump(doc, f, ensure_ascii=False)
                f.write('\n')
        logger.info(f"✓ Saved {len(documents)} documents to {output_file}")

# ========================================================================
# SEMANTIC LAYER: Stable Biological Entities (unchanged from v3)
# ========================================================================

class SemanticLayerBuilder(CollectionBuilder):
    """Build semantic layer collections (stable biological entities)."""

    def build_gene_nodes(self) -> List[Dict]:
        """
        Build GENES collection - semantic nodes with ontological annotations.
        One node per gene (not replicated per sample).
        """
        logger.info("Building GENES semantic collection...")

        biomart = self.mapping_dfs["biomart_genes"]
        genes = []

        for idx, row in tqdm(biomart.iterrows(), total=len(biomart), desc="Processing genes"):
            gene_id = row['gene_stable_id']

            # Parse transcript list if available
            transcript_ids = []
            if pd.notna(row.get('transcript_id')):
                try:
                    import ast
                    transcript_ids = ast.literal_eval(row['transcript_id']) if isinstance(row['transcript_id'], str) else [row['transcript_id']]
                except:
                    transcript_ids = []

            gene_doc = {
                "_key": gene_id,
                "entity_type": "gene",
                "bioentity_type": "gene",
                "gene_stable_id": gene_id,
                "gene_stable_id_version": row.get('gene_stable_id_version'),
                "hgnc_symbol": row.get('hgnc_symbol') if pd.notna(row.get('hgnc_symbol')) else None,
                "entrez_id": str(int(float(row['entrez_id']))) if pd.notna(row.get('entrez_id')) else None,
                "uniprot_id": row.get('uniprot_swissprot_id') if pd.notna(row.get('uniprot_swissprot_id')) else None,
                "mirbase_id": row.get('mirbase_id') if pd.notna(row.get('mirbase_id')) else None,
                "gene_type": row.get('gene_type'),
                "gene_description": row.get('gene_description'),
                "chromosome": row.get('chromosome'),
                "gene_start_bp": row.get('gene_start_bp'),
                "gene_end_bp": row.get('gene_end_bp'),
                "strand": row.get('strand'),
                "transcript_ids": transcript_ids,
                "source": "Ensembl_BioMart",
                "source_version": "GRCh38.v36"
            }

            genes.append(gene_doc)

        return genes

    def build_project_node(self, clinical_df: pd.DataFrame) -> List[Dict]:
        """Build PROJECT collection - one document per TCGA study."""
        logger.info("Building PROJECT collection...")

        import ast
        project_doc = {
            "_key": clinical_df.loc[0, 'project_id.project'],
            "name": clinical_df.loc[0, 'name.project'],
            "program": clinical_df.loc[0, 'name.program.project'],
            "primary_site": clinical_df.loc[0, 'primary_site.project'],
            "disease_types": ast.literal_eval(clinical_df.loc[0, 'disease_type.project']),
            "entity_type": "project",
            "n_cases": int(clinical_df['submitter_id'].nunique()),
            "n_samples": int(clinical_df['sample'].nunique())
        }

        return [project_doc]

# ========================================================================
# QUANTITATIVE INDEX LAYER: Gene/Feature Position Mappings
# ========================================================================

class QuantitativeIndexBuilder(CollectionBuilder):
    """Build index collections for interpreting quantitative vectors."""

    def build_expression_index(self, star_tpm_df: pd.DataFrame) -> List[Dict]:
        """
        Build EXPRESSION_INDEX collection - gene position mapping for expression vectors.

        This collection defines the ordered list of genes for the expression matrix.
        One document per cohort that maps gene_position → gene_id.
        """
        logger.info("Building EXPRESSION_INDEX collection...")

        # Extract gene IDs in order from the dataframe
        gene_ids_full = star_tpm_df['Ensembl_ID'].values
        gene_ids_base = [g.split('.')[0] for g in gene_ids_full]

        # Create mappings to semantic identifiers
        gene_mappings = []
        for i, (gene_full, gene_base) in enumerate(zip(gene_ids_full, gene_ids_base)):
            entrez_id = self.lookups["ensembl_to_entrez"].get(gene_base)
            hgnc_symbol = self.lookups["ensembl_to_hgnc"].get(gene_base)

            gene_mappings.append({
                "position": i,
                "gene_id_ensembl": gene_full,
                "gene_id_base": gene_base,
                "gene_ref": f"genes/{gene_base}",
                "entrez_id": str(int(float(entrez_id))) if pd.notna(entrez_id) else None,
                "hgnc_symbol": hgnc_symbol if pd.notna(hgnc_symbol) else None
            })

        # Create index document
        index_doc = {
            "_key": f"expr_index_{self.study}",
            "cohort": self.study,
            "data_type": "gene_expression_index",
            "n_genes": len(gene_ids_base),
            "platform": "STAR",
            "genome_version": "GRCh38.v36",
            "gene_mappings": gene_mappings,
            "description": "Gene position mapping for expression vectors in sample documents"
        }

        logger.info(f"  Created expression index with {len(gene_mappings)} genes")
        return [index_doc]

    def build_cnv_index(self, cnv_df: pd.DataFrame) -> List[Dict]:
        """Build CNV_INDEX collection - gene position mapping for CNV vectors."""
        logger.info("Building CNV_INDEX collection...")

        gene_ids_full = cnv_df['Ensembl_ID'].values
        gene_ids_base = [g.split('.')[0] for g in gene_ids_full]

        gene_mappings = []
        for i, (gene_full, gene_base) in enumerate(zip(gene_ids_full, gene_ids_base)):
            entrez_id = self.lookups["ensembl_to_entrez"].get(gene_base)
            enst_ids = self.lookups.get("ensg_to_ensts", {}).get(gene_base, [])

            gene_mappings.append({
                "position": i,
                "gene_id_ensembl": gene_full,
                "gene_id_base": gene_base,
                "gene_ref": f"genes/{gene_base}",
                "entrez_id": str(int(float(entrez_id))) if pd.notna(entrez_id) else None,
                "enst_ids": enst_ids if enst_ids else None
            })

        index_doc = {
            "_key": f"cnv_index_{self.study}",
            "cohort": self.study,
            "data_type": "cnv_index",
            "n_genes": len(gene_ids_base),
            "analysis_method": "ASCAT3",
            "gene_mappings": gene_mappings,
            "description": "Gene position mapping for CNV vectors in sample documents"
        }

        logger.info(f"  Created CNV index with {len(gene_mappings)} genes")
        return [index_doc]

    def build_mirna_index(self, mirna_df: pd.DataFrame) -> List[Dict]:
        """Build MIRNA_INDEX collection - miRNA position mapping."""
        logger.info("Building MIRNA_INDEX collection...")

        mirna_ids = mirna_df['miRNA_ID'].values
        mirna_hgcn = mirna_df['hgcn_id'].values if 'hgcn_id' in mirna_df.columns else [None]*len(mirna_ids)

        # Map to gene IDs if available
        mirna_mappings = []
        for i, mirna_id in enumerate(mirna_ids):
            mirna_mappings.append({
                "position": i,
                "mirna_id": mirna_id,
                "mirbase_id": mirna_id,
                "hgnc_symbol": mirna_hgcn[i] if pd.notna(mirna_hgcn[i]) else None,
                "description": f"miRNA {mirna_id}"
            })

        index_doc = {
            "_key": f"mirna_index_{self.study}",
            "cohort": self.study,
            "data_type": "mirna_index",
            "n_mirnas": len(mirna_ids),
            "platform": "Illumina",
            "mirna_mappings": mirna_mappings,
            "description": "miRNA position mapping for expression vectors in sample documents"
        }

        logger.info(f"  Created miRNA index with {len(mirna_ids)} miRNAs")
        return [index_doc]

    def build_protein_index(self, protein_df: pd.DataFrame) -> List[Dict]:
        """Build PROTEIN_INDEX collection - protein/peptide position mapping (RPPA)."""
        logger.info("Building PROTEIN_INDEX collection...")

        peptide_ids = protein_df['peptide_target'].values

        # Load RPPA mapping
        rppa_map = self.mapping_dfs.get("rppa")

        protein_mappings = []
        for i, peptide_id in enumerate(peptide_ids):
            # Try to map to gene
            if rppa_map is not None:
                match = rppa_map[rppa_map['RPPA_Protein_ID'] == peptide_id]
                if not match.empty:
                    entrez_id = match.iloc[0].get('Entrez_Gene_ID')
                    gene_symbol = match.iloc[0].get('Gene_Symbol')
                    protein_type = match.iloc[0].get('Protein_Type')
                else:
                    entrez_id = None
                    gene_symbol = None
                    protein_type = None
            else:
                entrez_id = None
                gene_symbol = None
                protein_type = None

            protein_mappings.append({
                "position": i,
                "peptide_target": peptide_id,
                "entrez_id": entrez_id,
                "gene_symbol": gene_symbol,
                "protein_type": protein_type
            })

        index_doc = {
            "_key": f"protein_index_{self.study}",
            "cohort": self.study,
            "data_type": "protein_index",
            "n_proteins": len(peptide_ids),
            "platform": "RPPA",
            "protein_mappings": protein_mappings,
            "description": "Protein/peptide position mapping for RPPA vectors in sample documents"
        }

        logger.info(f"  Created protein index with {len(peptide_ids)} peptides")
        return [index_doc]

# ========================================================================
# QUANTITATIVE VECTOR LAYER: Sample-Centric Storage
# ========================================================================

class SampleCentricVectorBuilder(CollectionBuilder):
    """Build sample-centric quantitative collections with vector storage."""

    def build_gene_expression_samples(self, 
                                      star_counts_df: pd.DataFrame = None,
                                      star_fpkm_df: pd.DataFrame = None,
                                      star_tpm_df: pd.DataFrame = None) -> List[Dict]:
        """
        Build GENE_EXPRESSION_SAMPLES collection - one document per sample.

        Strategy: Sample-centric with parallel vectors for counts, FPKM, TPM.
        Each sample document contains:
        - sample_id (key)
        - expression_index_ref (reference to index collection)
        - values_counts, values_fpkm, values_tpm (parallel vectors)

        Advantages:
        - ~1,000 documents instead of ~20,000,000
        - Direct matrix extraction: sample → vector
        - Efficient for ML/GNN feature engineering
        - Cohort-level queries with simple filters
        """
        logger.info("Building GENE_EXPRESSION_SAMPLES collection (sample-centric vectors)...")

        # Get sample IDs (columns, excluding first column which is Ensembl_ID)
        sample_ids = star_tpm_df.columns[1:].tolist()
        logger.info(f"  Processing {len(sample_ids)} samples")

        documents = []

        for sample_id in tqdm(sample_ids, desc="Building sample expression vectors"):
            # Extract vectors for this sample
            # gene_vectors = star_tpm_df['Ensembl_ID'].values
            if star_counts_df is not None:
                counts_vector = pd.to_numeric(star_counts_df[sample_id], errors='coerce').tolist()
            if star_fpkm_df is not None:
                fpkm_vector = pd.to_numeric(star_fpkm_df[sample_id], errors='coerce').tolist()

            tpm_vector = pd.to_numeric(star_tpm_df[sample_id], errors='coerce').tolist()

            

            # Create sample document
            sample_doc = {
                "_key": sample_id,
                "sample_id": sample_id,
                "cohort": self.study,
                "data_type": "gene_expression_vector",
                "expression_index_ref": f"expression_index/expr_index_{self.study}",
                "platform": "STAR",
                "n_genes": len(tpm_vector),
                # "gene_ids": gene_vectors.tolist(),  
                # Parallel vectors (same order as expression_index gene_mappings)
                "values_counts": counts_vector if star_counts_df is not None else None,
                "values_fpkm": fpkm_vector if star_fpkm_df is not None else None,
                "values_tpm": tpm_vector,
                "normalization": "TPM/FPKM/counts"
            }

            documents.append(sample_doc)

        logger.info(f"  Created {len(documents)} sample expression documents")
        return documents

    def build_cnv_samples(self, cnv_df: pd.DataFrame) -> List[Dict]:
        """Build CNV_SAMPLES collection - one document per sample with CNV vector."""
        logger.info("Building CNV_SAMPLES collection (sample-centric vectors)...")

        sample_ids = cnv_df.columns[1:].tolist()
        logger.info(f"  Processing {len(sample_ids)} samples")

        documents = []

        for sample_id in tqdm(sample_ids, desc="Building sample CNV vectors"):
            cnv_vector = pd.to_numeric(cnv_df[sample_id], errors='coerce').tolist()

            sample_doc = {
                "_key": sample_id,
                "sample_id": sample_id,
                "cohort": self.study,
                "data_type": "cnv_vector",
                "cnv_index_ref": f"cnv_index/cnv_index_{self.study}",
                "analysis_method": "ASCAT3",
                "n_genes": len(cnv_vector),
                "values_copy_number": cnv_vector
            }

            documents.append(sample_doc)

        logger.info(f"  Created {len(documents)} sample CNV documents")
        return documents

    def build_mirna_samples(self, mirna_df: pd.DataFrame) -> List[Dict]:
        """Build MIRNA_SAMPLES collection - one document per sample with miRNA vector."""
        logger.info("Building MIRNA_SAMPLES collection (sample-centric vectors)...")

        sample_ids = mirna_df.columns[1:].tolist()
        logger.info(f"  Processing {len(sample_ids)} samples")

        documents = []

        for sample_id in tqdm(sample_ids, desc="Building sample miRNA vectors"):
            mirna_vector = pd.to_numeric(mirna_df[sample_id], errors='coerce').tolist()

            sample_doc = {
                "_key": sample_id,
                "sample_id": sample_id,
                "cohort": self.study,
                "data_type": "mirna_vector",
                "mirna_index_ref": f"mirna_index/mirna_index_{self.study}",
                "platform": "Illumina",
                "n_mirnas": len(mirna_vector),
                "values_expression": mirna_vector
            }

            documents.append(sample_doc)

        logger.info(f"  Created {len(documents)} sample miRNA documents")
        return documents

    def build_protein_samples(self, protein_df: pd.DataFrame) -> List[Dict]:
        """Build PROTEIN_SAMPLES collection - one document per sample with RPPA vector."""
        logger.info("Building PROTEIN_SAMPLES collection (sample-centric vectors)...")

        sample_ids = protein_df.columns[1:].tolist()
        logger.info(f"  Processing {len(sample_ids)} samples")

        documents = []

        for sample_id in tqdm(sample_ids, desc="Building sample protein vectors"):
            protein_vector = pd.to_numeric(protein_df[sample_id], errors='coerce').tolist()

            sample_doc = {
                "_key": sample_id,
                "sample_id": sample_id,
                "cohort": self.study,
                "data_type": "protein_vector",
                "protein_index_ref": f"protein_index/protein_index_{self.study}",
                "platform": "RPPA",
                "n_proteins": len(protein_vector),
                "values_abundance": protein_vector
            }

            documents.append(sample_doc)

        logger.info(f"  Created {len(documents)} sample protein documents")
        return documents

# ========================================================================
# CLINICAL/SAMPLE METADATA LAYER (unchanged from v3)
# ========================================================================

class MetadataLayerBuilder(CollectionBuilder):
    """Build clinical and sample metadata collections."""

    def build_samples(self, clinical_df: pd.DataFrame) -> List[Dict]:
        """Build SAMPLES collection."""
        logger.info("Building SAMPLES collection...")

        samples = []
        for idx, row in tqdm(clinical_df.iterrows(), total=len(clinical_df), desc="Processing samples"):
            sample_doc = {
                "_key": row['sample'],
                "submitter_id": row['submitter_id'],
                "case_ref": f"cases/{row['submitter_id']}",
                "project_ref": f"projects/{row['project_id.project']}",
                "sample_type": row.get('sample_type.samples'),
                "sample_type_id": row.get('sample_type_id.samples'),
                "tissue_type": row.get('tissue_type.samples'),
                "tumor_descriptor": row.get('tumor_descriptor.samples'),
                "specimen_type": row.get('specimen_type.samples'),
                "composition": row.get('composition.samples'),
                "preservation_method": row.get('preservation_method.samples'),
                "days_to_collection": row.get('days_to_collection.samples'),
                "entity_type": "sample"
            }
            samples.append(sample_doc)

        return samples

    def build_cases(self, clinical_df: pd.DataFrame, survival_df: pd.DataFrame) -> List[Dict]:
        """Build CASES collection with clinical and survival data."""
        logger.info("Building CASES collection...")

        # Group by patient (submitter_id) to avoid duplicates
        patient_groups = clinical_df.groupby('submitter_id')

        cases = []
        for patient_id, group in tqdm(patient_groups, desc="Processing cases"):
            row = group.iloc[0]  # Take first row for patient-level data

            # Find survival data
            sample_id = row['sample']
            survival_match = survival_df[survival_df['sample'] == sample_id]

            case_doc = {
                "_key": patient_id,
                "project_ref": f"projects/{row['project_id.project']}",
                "primary_site": row.get('primary_site'),
                "disease_type": row.get('disease_type'),
                # "mondo_id": row.get('mondo_id'),
                "demographic": {
                    "gender": row.get('gender.demographic'),
                    "race": row.get('race.demographic'),
                    "ethnicity": row.get('ethnicity.demographic'),
                    "vital_status": row.get('vital_status.demographic'),
                    "age_at_index": row.get('age_at_index.demographic'),
                    "days_to_birth": row.get('days_to_birth.demographic'),
                    "year_of_birth": row.get('year_of_birth.demographic'),
                    "year_of_death": row.get('year_of_death.demographic'),
                    "days_to_death": row.get('days_to_death.demographic')
                },
                "diagnoses": {
                    "primary_diagnosis": row.get('primary_diagnosis.diagnoses'),
                    "age_at_diagnosis": row.get('age_at_diagnosis.diagnoses'),
                    "tumor_grade": row.get('tumor_grade.diagnoses'),
                    "ajcc_pathologic_stage": row.get('ajcc_pathologic_stage.diagnoses'),
                    "ajcc_pathologic_t": row.get('ajcc_pathologic_t.diagnoses'),
                    "ajcc_pathologic_n": row.get('ajcc_pathologic_n.diagnoses'),
                    "ajcc_pathologic_m": row.get('ajcc_pathologic_m.diagnoses')
                },
                "entity_type": "case"
            }

            # Add survival data if available
            if not survival_match.empty:
                case_doc["survival"] = {
                    "overall_survival_time": float(survival_match.iloc[0]['OS.time']) if pd.notna(survival_match.iloc[0]['OS.time']) else None,
                    "overall_survival_status": int(survival_match.iloc[0]['OS']) if pd.notna(survival_match.iloc[0]['OS']) else None
                }

            cases.append(case_doc)

        return cases

# ========================================================================
# MAIN EXECUTION PIPELINE v4
# ========================================================================

def main():
    """Main execution pipeline for loading TCGA data into ArangoDB collections (v4)."""

    logger.info("=" * 70)
    logger.info("TCGA Data Loading Pipeline v4.0 - Sample-Centric Vector Strategy")
    logger.info("=" * 70)
    logger.info(f"Study: {STUDY}")
    logger.info(f"Output: {OUTPUT_PATH}")
    logger.info("")

    # STEP tester
    STEP1 = True  # Load Mappings
    STEP2 = False  # Semantic Layer - GENES
    STEP3 = False  # Semantic Layer - PROJECT
    STEP4 = False  # Metadata Layer - SAMPLES & CASES
    STEP5 = False  # Quantitative Index Layer - EXPRESSION_INDEX
    STEP6 = True  # Quantitative Vector Layer - GENE_EXPRESSION_SAMPLES
    STEP7 = False  # CNV Index & Vectors
    STEP8 = False  # miRNA Index & Vectors
    STEP9 = False  # Protein Index & Vectors
    ALL_STEPS = False  # If True, run all steps



    # ========== STEP 1: Load Mappings ==========
    logger.info("[Step 1/9] Loading ID mappings...")
    mapping_dfs = load_mappings()
    lookups = create_lookup_dicts(mapping_dfs)
    logger.info("")

    # initialize builders
    index_builder = QuantitativeIndexBuilder(STUDY, mapping_dfs, lookups)
    vector_builder = SampleCentricVectorBuilder(STUDY, mapping_dfs, lookups)

    # ========== STEP 2: Semantic Layer - GENES ==========
    if STEP2 or ALL_STEPS:
        logger.info("[Step 2/9] Building Semantic Layer - GENES...")
        semantic_builder = SemanticLayerBuilder(STUDY, mapping_dfs, lookups)
        genes_collection = semantic_builder.build_gene_nodes()
        semantic_builder.save_collection(genes_collection, "genes")
        logger.info("")

    # ========== STEP 3: Semantic Layer - PROJECT ==========
    if STEP3 or ALL_STEPS:
        logger.info("[Step 3/9] Building Semantic Layer - PROJECT...")
        clinical_df = load_omics_data("clinical")
        # add to clinical_df , mondo_id column from mapping file
        if clinical_df is not None and "mondo" in mapping_dfs:
            mondo_map = mapping_dfs["mondo"]
            clinical_df = clinical_df.merge(mondo_map[['_disease_type', 'mondo_id']], 
                                            left_on='disease_type.project', 
                                            right_on='_disease_type', 
                                            how='left')
            clinical_df.drop(columns=['_disease_type'], inplace=True)
            print(f"  Merged mondo_id into clinical_df, new shape: {clinical_df.shape}")
        if clinical_df is not None:
            project_collection = semantic_builder.build_project_node(clinical_df)
            semantic_builder.save_collection(project_collection, "projects")
        logger.info("")

    # ========== STEP 4: Metadata Layer - SAMPLES & CASES ==========
    if STEP4 or ALL_STEPS:
        logger.info("[Step 4/9] Building Metadata Layer - SAMPLES & CASES...")
        if clinical_df is not None:
            survival_df = load_omics_data("survival")
            metadata_builder = MetadataLayerBuilder(STUDY, mapping_dfs, lookups)

            samples_collection = metadata_builder.build_samples(clinical_df)
            metadata_builder.save_collection(samples_collection, "samples")

            if survival_df is not None:
                cases_collection = metadata_builder.build_cases(clinical_df, survival_df)
                metadata_builder.save_collection(cases_collection, "cases")
        logger.info("")

    # ========== STEP 5: Quantitative Index Layer - EXPRESSION_INDEX ==========
    
    if STEP5 or ALL_STEPS:
        logger.info("[Step 5/9] Building Quantitative Index Layer - EXPRESSION_INDEX...")
        star_tpm = load_omics_data("star_tpm")
        if star_tpm is not None:
            # index_builder = QuantitativeIndexBuilder(STUDY, mapping_dfs, lookups)
            expr_index = index_builder.build_expression_index(star_tpm)
            index_builder.save_collection(expr_index, "expression_index")
        logger.info("")

    # ========== STEP 6: Quantitative Vector Layer - GENE_EXPRESSION_SAMPLES ==========
    if STEP6 or ALL_STEPS:
        logger.info("[Step 6/9] Building Quantitative Vector Layer - GENE_EXPRESSION_SAMPLES...")
        star_tpm = load_omics_data("star_tpm")
        LITE_VERSION = True  # If True, only load TPM data
        if LITE_VERSION:
            star_counts = pd.DataFrame()
            star_fpkm = pd
        else:
            star_counts = load_omics_data("star_counts")
            star_fpkm = load_omics_data("star_fpkm")

        if all([star_counts is not None, star_fpkm is not None, star_tpm is not None]):
            # vector_builder = SampleCentricVectorBuilder(STUDY, mapping_dfs, lookups)
            if LITE_VERSION:
                expr_samples = vector_builder.build_gene_expression_samples(None, None, star_tpm)
            else:
                expr_samples = vector_builder.build_gene_expression_samples(star_counts, star_fpkm, star_tpm)
            print(f"  Built gene expression samples, total documents: {len(expr_samples)}")
            print("Saving gene expression samples collection...")
            vector_builder.save_collection(expr_samples, f"gene_expression_samples_{STUDY}")



        logger.info("")

    # ========== STEP 7: CNV Index & Vectors ==========
    if STEP7 or ALL_STEPS:
        logger.info("[Step 7/9] Building CNV Index & Vectors...")
        cnv_df = load_omics_data("gene-level_ascat3")
        if cnv_df is not None:
            cnv_index = index_builder.build_cnv_index(cnv_df)
            index_builder.save_collection(cnv_index, "cnv_index")

            cnv_samples = vector_builder.build_cnv_samples(cnv_df)
            vector_builder.save_collection(cnv_samples, f"cnv_samples_{STUDY}")
        logger.info("")

    # ========== STEP 8: miRNA Index & Vectors ==========
    if STEP8 or ALL_STEPS:
        logger.info("[Step 8/9] Building miRNA Index & Vectors...")
        mirna_df = load_omics_data("mirna")
        # add to mirna_df , hgcn column from mapping file
        if mirna_df is not None and "mirna" in mapping_dfs:
            mirna_map = mapping_dfs["mirna"]
            #hgcn_id	miRNA_ID
            mirna_df = mirna_df.merge(mirna_map[['miRNA_ID', 'hgcn_id']], 
                                    on='miRNA_ID', 
                                    how='left')
            print(f"  Merged HGNC_Symbol into mirna_df, new shape: {mirna_df.shape}")
        "$$$"
        if mirna_df is not None:
            mirna_index = index_builder.build_mirna_index(mirna_df)
            index_builder.save_collection(mirna_index, "mirna_index")

            mirna_samples = vector_builder.build_mirna_samples(mirna_df)
            vector_builder.save_collection(mirna_samples, f"mirna_samples_{STUDY}")
        logger.info("")

    # ========== STEP 9: Protein Index & Vectors ==========
    if STEP9 or ALL_STEPS:
        logger.info("[Step 9/9] Building Protein Index & Vectors...")
        protein_df = load_omics_data("protein")
        if protein_df is not None:
            protein_index = index_builder.build_protein_index(protein_df)
            index_builder.save_collection(protein_index, "protein_index")

            protein_samples = vector_builder.build_protein_samples(protein_df)
            vector_builder.save_collection(protein_samples, f"protein_samples_{STUDY}")
        logger.info("")

    # STEP 10 load all to ArangoDB usig python-arango or arangoimport
    logger.info("[Step 10/10] Loading JSON files into ArangoDB...")
    logger.info("  (This step is not implemented in this script. Use arangoimport or python-arango to load the JSON files.)")
    logger.info("")


    # ========== Summary ==========
    logger.info("=" * 70)
    logger.info("✓ Pipeline v4.0 completed successfully!")
    logger.info("=" * 70)
    logger.info(f"Output directory: {OUTPUT_PATH}")
    logger.info("")
    logger.info("Collection Summary:")
    logger.info("  Semantic Layer:")
    logger.info("    - genes.json (~20K documents)")
    logger.info("    - projects.json (1 document)")
    logger.info("  Metadata Layer:")
    logger.info("    - samples.json (~1K documents)")
    logger.info("    - cases.json (~1K documents)")
    logger.info("  Quantitative Index Layer:")
    logger.info("    - expression_index.json (1 document with gene mappings)")
    logger.info("    - cnv_index.json (1 document)")
    logger.info("    - mirna_index.json (1 document)")
    logger.info("    - protein_index.json (1 document)")
    logger.info("  Quantitative Vector Layer:")
    logger.info(f"    - gene_expression_samples_{STUDY}.json (~1K documents with vectors)")
    logger.info(f"    - cnv_samples_{STUDY}.json (~1K documents with vectors)")
    logger.info(f"    - mirna_samples_{STUDY}.json (~1K documents with vectors)")
    logger.info(f"    - protein_samples_{STUDY}.json (~1K documents with vectors)")
    logger.info("")
    logger.info("Next Steps:")
    logger.info("  1. Import JSON files into ArangoDB:")
    logger.info("     arangoimport --file genes.json --collection genes --type jsonl")
    logger.info("  2. Create indexes:")
    logger.info("     - Index on sample_id for *_samples collections")
    logger.info("     - Index on cohort for multi-cohort queries")
    logger.info("  3. Query examples:")
    logger.info("     - Get sample vector: FOR s IN gene_expression_samples_TCGA_BRCA")
    logger.info("                          FILTER s.sample_id == 'TCGA-A1-A0SB-01A'")
    logger.info("                          RETURN s.values_tpm")
    logger.info("     - Get gene across cohort: see index gene_mappings + extract from vectors")
    logger.info("")
    logger.info("Advantages of v4 vs v3:")
    logger.info("  - Documents: ~1,000 vs ~20,000,000 (20,000x reduction!)")
    logger.info("  - Storage: ~100MB vs ~10GB per cohort")
    logger.info("  - Query speed: O(1) sample lookup vs O(n_genes) scan")
    logger.info("  - ML-ready: direct vector extraction for feature matrices")
    logger.info("=" * 70)

if __name__ == "__main__":
    main()
