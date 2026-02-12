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
  - miRNA expression
  - Proteomics

- **ID Normalization**: Ensembl, Entrez, HGNC mappings
- **Clinical Data**: Survival, treatment response

### 3. Data Representation
Multi-omics measurements are represented as properties of knowledge graph nodes:
- Gene expression → Gene/Transcript nodes
- Mutations → Gene/Protein nodes
- Methylation → Gene/Promoter nodes
- CNV → Gene/Chromosome region nodes
- miRNA expression → miRNA nodes
- Proteomics → Protein nodes

Detailed data structure are described in the [Database Structure Documentation](docs/readme_db_structure.md).
<!-- add a link to docs/readme_db_structure.md -->


## Key Scripts

1. `scripts\build_property_graph.py`: Converts PKT TSV tables to ArangoDB-compatible JSON format
2. `scripts\load_graph_to_arangodb.py`: Loads JSON graph data into ArangoDB
3. `scripts\build_omics_collections.py`: Builds ArangoDB collections for multi-omics data
4. `scripts\load_omics_collections_to_arangodb.py`: Integrates TCGA multi-omics data into the knowledge graph
5. `scripts\build_kg_transomic.py`: Builds sample-specific trans-omic networks for a given TCGA sample
6. `scripts\analyze_kg_transomics.py`: Analyzes and summarizes properties of the trans-omic subgraphs


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

## Versioning
- Current Version: v0.0.1 (alpha)

## License

[Add your license here]

## Citation

[Add citation information when published]

## Contact

[Add contact information]