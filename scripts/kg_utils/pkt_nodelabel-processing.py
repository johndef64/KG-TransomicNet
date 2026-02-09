#%%
"""
Simple functional script to read and process the PKT NodeLabels file.
Reads the NodeLabels tar.gz (tab-separated), converts integer_id to string,
extracts an `entity` column from `entity_uri`, computes `class_code` and
applies the same replacements as the original `read-edit.py`, and writes a CSV.

Usage:
    python read_node_labels_simple.py --input PATH_TO/PKT_NodeLabels.txt.tar.gz --output out.csv

Defaults assume the file at ../data/pkt/builds/v3.0.2/PKT_NodeLabels.txt.tar.gz
relative to this script.
"""

import argparse
import os
import tarfile
import pandas as pd
import re
from  pkt_utils import *
set_working_directory()


def read_node_labels_from_tar(tar_path):
    if not os.path.exists(tar_path):
        raise FileNotFoundError(f"Input tar.gz not found: {tar_path}")

    with tarfile.open(tar_path, 'r:gz') as tar:
        members = tar.getmembers()
        if not members:
            raise ValueError("Tar file is empty")
        # assume the first member is the NodeLabels text file
        extracted_file = tar.extractfile(members[0])
        if extracted_file is None:
            raise ValueError("Could not extract member from tar")
        # read tab-separated file
        df = pd.read_csv(extracted_file, sep='\t', dtype=str)
        return df


def extract_entity_from_uri(uri):
    if pd.isna(uri):
        return ''
    s = uri.rstrip('>')
    # take last segment after '/' then after '=' if present
    part = s.split('/')[-1]
    part = part.split('=')[-1]
    return part


def extract_entity_class(entity):
    # remove trailing numbers (but not if they are part of a larger digit sequence preceding digits?)
    s = re.sub(r'(?<!\d)\d+$', '', entity).rstrip('_-')
    for prefix in ('PR_', 'CHR_', 'GNO_'):
        if s.startswith(prefix):
            return prefix.rstrip('_')
    # remove anything after a '#' if present
    return s.split('#')[0]


def process_node_labels(df):
    # ensure expected columns exist
    # common columns seen in original script: entity_type, integer_id, entity_uri, label, description/definition, synonym
    # convert integer_id to string (keep NaNs as empty strings)
    if 'integer_id' in df.columns:
        df['integer_id'] = df['integer_id'].astype(str)

    # create entity from entity_uri
    if 'entity_uri' in df.columns:
        df['entity'] = df['entity_uri'].apply(extract_entity_from_uri)
    else:
        df['entity'] = ''

    # compute class_code
    df['class_code'] = df['entity'].apply(extract_entity_class)

    # replacements from original script
    df['class_code'] = df['class_code'].replace('UMLS_C', 'UMLS')
    df['class_code'] = df['class_code'].replace('22-rdf-syntax-ns', 'RDF')
    # drop entity_class if present
    if 'entity_class' in df.columns:
        df = df.drop(columns=['entity_class'])

    # replace empty strings with EntrezID (only exact empty strings)
    df['class_code'] = df['class_code'].replace('', 'EntrezID')
    # replace exact 'rs' entries with 'dbSNP'
    df['class_code'] = df['class_code'].replace('rs', 'dbSNP')

    return df

#%%

# Add class code to NodeLabels and save to CSV

input_path = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'builds', 'v3.0.2', 'PKT_NodeLabels.txt.tar.gz'))
out_path = 'PKT_NodeLabels_with_class_code_v3.0.2.csv'

if not os.path.exists(out_path):
    NodeLabels = read_node_labels_from_tar(input_path)
    NodeLabels = process_node_labels(NodeLabels)
    NodeLabels.to_csv(out_path, index=False)
else:
    NodeLabels = pd.read_csv(out_path)
    print(NodeLabels.head())

#%%

# load simple the output datqaset
import pandas as pd
NodeLabels = pd.read_csv("PKT_NodeLabels_with_class_code_v3.0.2.csv")
NodeLabels.head()
#%%
# describe NodeLabels.class_code entries groupby
NodeLabels.groupby('class_code').size().sort_values(ascending=False)

#%%
edge_list_dict = [
    # ── 1. PROTEIN-CENTRIC INTERACTIONS ─────────────────────────────
    'protein-protein',      # direct physical or functional links
    'protein-catalyst',     # enzyme ↔ substrate/product
    'protein-cofactor',     # enzyme ↔ helper molecule
    'protein-gobp',         # protein ↔ biological process
    'protein-gomf',         # protein ↔ molecular function
    'protein-gocc',         # protein ↔ cellular component
    'protein-anatomy',      # protein ↔ tissue / organ
    'protein-cell',         # protein ↔ cell type
    'protein-pathway',      # protein ↔ pathway membership

    # ── 2. GENE-CENTRIC INTERACTIONS ────────────────────────────────
    'gene-gene',            # co-expression / epistasis
    'gene-protein',         # gene ↔ its product
    'gene-rna',             # gene ↔ transcript
    'gene-phenotype',       # gene ↔ observable trait
    'gene-disease',         # gene ↔ disorder
    'gene-pathway',         # gene ↔ pathway

    # ── 3. CHEMICAL-CENTRIC INTERACTIONS ─────────────────────────────
    'chemical-protein',     # drug ↔ target
    'chemical-gene',      # compound ↔ regulator gene
    'chemical-disease',     # drug / metabolite ↔ disease
    'chemical-phenotype',   # chemical ↔ phenotypic effect
    'chemical-pathway',   # compound ↔ pathway
    'chemical-gobp',        # chemical ↔ biological process
    'chemical-gomf',        # chemical ↔ molecular function
    'chemical-gocc',        # chemical ↔ cellular component

    # ── 4. VARIANT-CENTRIC INTERACTIONS ──────────────────────────────
    'variant-gene',         # mutation ↔ affected gene
    'variant-disease',        # mutation ↔ disease
    'variant-phenotype',      # mutation ↔ phenotypic outcome

    # ── 5. RNA-CENTRIC INTERACTIONS ─────────────────────────────────
    'rna-protein',          # RNA-binding protein interactions
    'rna-anatomy',          # RNA ↔ tissue localization
    'rna-cell',             # RNA ↔ cell-type expression

    # ── 6. ONTOLOGICAL / HIGHER-ORDER RELATIONS ──────────────────────
    'pathway-gocc',         # pathway ↔ cellular location
    'pathway-gomf',        # pathway ↔ molecular function
    'gobp-pathway',         # biological process ↔ pathway
    'disease-phenotype',    # disease ↔ clinical trait
]


#%%
# I want a mini dataset where for each code_class there are 3 examples of entity

# NodeLabels[['entity',
#     'class_code']].drop_duplicates()

# ...existing code...

# volgio un mini dataset dove per ogni code_class ci sono 3 esempi di enity
# mini_dataset = (NodeLabels[['entity', 'class_code']]
#                 .drop_duplicates()
#                 .groupby('class_code')
#                 .head(3)
#                 .reset_index(drop=True))

# print(f"Mini dataset shape: {mini_dataset.shape}")
# print(f"Number of unique class_codes: {mini_dataset['class_code'].nunique()}")
# print(f"Examples per class_code:")
# # print(mini_dataset.groupby('class_code').size().head(10))

# # Preview the mini dataset
# print("\n--- Mini dataset preview ---")
# print(mini_dataset.head(15))

# # Save mini dataset if needed
# mini_dataset.to_csv("LabelClassDataset.csv", index=False)

# enrich manually the metadata for each class_code
# %%

metadata_label = pd.read_csv("LabelClassDataset_metadata.csv")
# meta= metadata_label[['class_code', 'entity_type', 'source', 'source_type']].drop_duplicates()
# meta.to_csv("LabelClassDataset_metadata.csv", index=False)

# %%
# add metadata on NodeLabels merging on class_code
NodeLabels_meta = NodeLabels.merge(metadata_label, on='class_code', how='left')
NodeLabels_meta.head()
#%%
NodeLabels_meta.to_csv("PKT_NodeLabels_with_metadata_v3.0.2.csv", index=False)
#%%
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

def plot_entity_type_counts(df = NodeLabels_meta, column = 'bioentity_type', save_=False):
    entity_type_counts = df[column].value_counts()
    # adapt height to number of unique types
    plt_height = max(4, len(entity_type_counts) * 0.3)

    plt.figure(figsize=(10, plt_height))
    ax = entity_type_counts.plot(kind='barh', color='steelblue')
    plt.title('Bioentity Types')
    plt.xlabel('Count')
    plt.gca().invert_yaxis()
    sns.despine(left=True, bottom=True)

    os.makedirs("../plots", exist_ok=True)
    plt.tight_layout()
    if save_:
        plt.savefig("../plots/bioentity_type_counts.png", dpi=300)
    plt.show()

plot_entity_type_counts(save_=True)
# %%.

# split NodeLabels_meta where entity_type is NODE or else
NodeLabels_nodes = NodeLabels_meta[NodeLabels_meta['entity_type'] == 'NODES']
NodeLabels_relations = NodeLabels_meta[NodeLabels_meta['entity_type'] == 'RELATIONS']
NodeLabels_others = NodeLabels_meta[~NodeLabels_meta['entity_type'].isin(['NODES', 'RELATIONS'])]
print(f"NodeLabels_nodes shape: {NodeLabels_nodes.shape}")
print(f"NodeLabels_relations shape: {NodeLabels_relations.shape}")

NodeLabels_relations
#%%
print(f"Relations total counts is:{NodeLabels_relations.shape[0]}")
plot_entity_type_counts(df=NodeLabels_relations, save_=True)
print(f"Nodes total counts is:{NodeLabels_nodes.shape[0]}")
plot_entity_type_counts(df=NodeLabels_nodes, save_=True)
print(f"Others total counts is:{NodeLabels_others.shape[0]}")
plot_entity_type_counts(df=NodeLabels_others, save_=False)
#%%
plot_entity_type_counts(df=NodeLabels_relations, column='class_code', save_=True)
plot_entity_type_counts(df=NodeLabels_nodes, column='class_code', save_=True)
plot_entity_type_counts(df=NodeLabels_others, column='class_code', save_=False)