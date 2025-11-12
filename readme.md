# KG-TransomicNet

A semantic approach for multi-omics data representation and reuse in a knowledge graph for systems biology and precision medicine applications.

## Overview

KG-TransomicNet integrates multi-omics data into a knowledge graph infrastructure based on ArangoDB, enabling exploratory analysis and descriptive statistics for cancer research.

## Project Goals

- Build a multi-omics knowledge graph using PheKnowLator (PKT) as foundation
- Integrate TCGA multi-omics datasets from Xena/cBioPortal
- Enable semantic representation and reuse of heterogeneous biological data
- Support exploratory analysis and visualization of integrated multi-omics data

## Architecture

### 1. Knowledge Graph Construction
- **Base KG**: PheKnowLator (PKT) converted from RDF triples to ArangoDB property graph
- **Build Configuration**: OWLNETS, Instance-based (Abox), Inverse relationships
- **Storage**: ArangoDB graph database

### 2. Multi-Omics Data Integration
- **Data Sources**: TCGA data from Xena browser and cBioPortal
- **Data Types**:
  - Gene expression
  - Somatic mutations
  - DNA methylation
  - Copy number variations
  - Proteomics
  - Metabolomics
- **ID Normalization**: Ensembl, Entrez, HGNC mappings
- **Clinical Data**: Survival, treatment response

### 3. Data Representation
Multi-omics measurements are represented as properties of knowledge graph nodes:
- Gene expression → Gene/Transcript nodes
- Mutations → Gene/Protein nodes
- Methylation → Gene/Promoter nodes
- CNV → Gene/Chromosome region nodes
- Proteomics → Protein nodes
- Metabolomics → Metabolite nodes

## Key Scripts

- `kg_to_arangodb.py`: Converts PKT TSV tables to ArangoDB-compatible JSON format

## Methodology

The approach follows the **IRR (Information Representation and Reuse)** framework for:
1. Semantic integration of heterogeneous biological data
2. Ontology-based data normalization
3. Graph-based exploratory analysis

## Requirements

- ArangoDB
- Python 3.x
- PheKnowLator dataset
- TCGA data from Xena/cBioPortal

## License

[Add your license here]

## Citation

[Add citation information when published]

## Contact

[Add contact information]