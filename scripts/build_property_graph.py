#%%

"""
Docstring for scripts.build_property_graph
Version 0.2 
"""

import pandas as pd
import json
import re
from collections import defaultdict
from urllib.parse import urlparse
import sys, os
import xml.etree.ElementTree as ET
from tqdm import tqdm

from  kg_utils.pkt_utils import *
# print_file_contents(os.getcwd())
#%%

# --- CONFIGURAZIONE ---
pkt_build_dir = "../data/pkt/builds/v3.0.2/"
pkt_file = f"{pkt_build_dir}PKT.nt.tar.gz"
node_label_file = f"{pkt_build_dir}PKT_NodeLabels_with_metadata_v3.0.2.csv"
output_dir = "../temp_dir/" # Assicurati che questa directory esista!
os.makedirs(output_dir, exist_ok=True)

USE_SAMPLE = True # Usa un campione di 10.000 triple per il test rapido
if USE_SAMPLE:
    print("USING SAMPLE OF 10,000 TRIPLES FOR TESTING PURPOSES")
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")    
    output_dir = f"../temp_dir/sample{timestamp}/"


""" 
PKT RDF triples preview:
subject	predicate	object
0	<http://purl.obolibrary.org/obo/PR_000003305>	<http://purl.obolibrary.org/obo/RO_0002436>	<http://purl.obolibrary.org/obo/PR_000003304>
1	<http://purl.obolibrary.org/obo/PR_000003305>	<http://purl.obolibrary.org/obo/RO_0002436>	<http://purl.obolibrary.org/obo/PR_000003306>

NodeLabels contents preview:
	entity_type	integer_id	entity_uri	label	description/definition	synonym	entity	class_code	bioentity_type	source	source_type
0	NODES	536450	<https://www.ncbi.nlm.nih.gov/snp/rs876660433>	NM_000535.7(PMS2):c.567T>C (p.His189=)	This variant is a germline single nucleotide variant located on chromosome 7 (NC_000007.14, start:5999246/stop:5999246 positions, cytogenetic location:7p22.1) and has clinical significance 'Likely benign'. This entry is for the GRCh38 and was last reviewed on May 03, 2020 with review status 'criteria provided, multiple submitters, no conflicts'.		rs876660433	dbSNP	variant	Database of Single Nucleotide Polymorphisms	Database
1	NODES	338642	<https://www.ncbi.nlm.nih.gov/snp/rs35349917>	NM_001105206.3(LAMA4):c.280G>A (p.Gly94Ser)	This variant is a germline/unknown single nucleotide variant located on chromosome 6 (NC_000006.12, start:112216385/stop:112216385 positions, cytogenetic location:6q21) and has clinical significance 'Benign'. This entry is for the GRCh38 and was last reviewed on Dec 07, 2020 with review status 'criteria provided, multiple submitters, no conflicts'.		rs35349917	dbSNP	variant	Database of Single Nucleotide Polymorphisms	Database
2	NODES	253428	<http://purl.obolibrary.org/obo/CHEBI_177646>	3,4,15-Trihydroxy-12,13-epoxytrichothec-9-en-8-yl 3-methylbutanoate		[(1S,2R,4S,7R,9R,10R,11S)-10,11-dihydroxy-2-(hydroxymethyl)-1,5-dimethylspiro[8-oxatricyclo[7.2.1.02,7]dodec-5-ene-12,2'-oxirane]-4-yl] 3-methylbutanoate	CHEBI_177646	CHEBI	chemical	Chemical Entities of Biological Interest	Ontology
2923	RELATIONS	247	<http://purl.obolibrary.org/obo/RO_0002513>	ribosomally translates to	inverse of ribosomal translation of		RO_0002513	RO	phenotype	Relation Ontology	Ontology
4033	RELATIONS	309922	<http://purl.obolibrary.org/obo/BFO_0000067>	contains process	[copied from inverse property 'occurs in'] b occurs_in c =def b is a process and c is a material entity or immaterial entity& there exists a spatiotemporal region r and b occupies_spatiotemporal_region r.& forall(t) if b exists_at t then c exists_at t & there exist spatial regions s and s’ where & b spatially_projects_onto s at t& c is occupies_spatial_region s’ at t& s is a proper_continuant_part_of s’ at t		BFO_0000067	BFO	go	Basic Formal Ontology	Ontology
10123	RELATIONS	167760	<http://purl.obolibrary.org/obo/RO_0002568>	has muscle antagonist	m1 has_muscle_antagonist m2 iff m1 acts in opposition to m2, and m2 is responsible for returning the structure to its initial position.		RO_0002568	RO	phenotype	Relation Ontology	Ontology
"""

# --- METADATA LOOKUP ---

def create_metadata_lookup(nodelabels_df):
    """
    Create a lookup dictionary from NodeLabels metadata.
    Maps entity_uri (without angle brackets) to all metadata fields.
    """
    print("Creating metadata lookup from NodeLabels...")
    
    # Remove angle brackets from entity_uri if present
    nodelabels_df['clean_uri'] = nodelabels_df['entity_uri'].str.strip('<>')
    
    # Create lookup dictionary: URI -> metadata
    metadata_lookup = {}
    for _, row in nodelabels_df.iterrows():
        uri = row['clean_uri']
        metadata_lookup[uri] = {
            'entity_type_cat': row['entity_type'],  # NODES or RELATIONS
            'integer_id': row['integer_id'],
            'label': row['label'] if pd.notna(row['label']) else '',
            'description': row['description/definition'] if pd.notna(row['description/definition']) else '',
            'synonym': row['synonym'] if pd.notna(row['synonym']) else '',
            'entity': row['entity'] if pd.notna(row['entity']) else '',
            'class_code': row['class_code'] if pd.notna(row['class_code']) else '',
            'bioentity_type': row['bioentity_type'] if pd.notna(row['bioentity_type']) else '',
            'source': row['source'] if pd.notna(row['source']) else '',
            'source_type': row['source_type'] if pd.notna(row['source_type']) else ''
        }
    
    print(f"Created metadata lookup with {len(metadata_lookup)} entries")
    return metadata_lookup


def get_entity_metadata(uri, metadata_lookup):
    """
    Get metadata for an entity from the metadata lookup.
    If not found, extract basic info from URI.
    """
    # Remove angle brackets if present
    clean_uri = uri.strip('<>')
    
    if clean_uri in metadata_lookup:
        meta = metadata_lookup[clean_uri]
        # Extract namespace from URI
        parsed = urlparse(clean_uri)
        namespace = parsed.netloc if parsed.netloc else "unknown"
        
        return {
            'uri': clean_uri,
            'namespace': namespace,
            'entity_id': meta['entity'],
            'class_code': meta['class_code'],
            'label': meta['label'],
            'bioentity_type': meta['bioentity_type'],
            'description': meta['description'],
            'synonym': meta['synonym'],
            'source': meta['source'],
            'source_type': meta['source_type'],
            'integer_id': meta['integer_id']
        }
    else:
        # Fallback: extract basic info from URI
        parsed = urlparse(clean_uri)
        namespace = parsed.netloc if parsed.netloc else "unknown"
        entity_id = clean_uri.split('/')[-1]
        
        return {
            'uri': clean_uri,
            'namespace': namespace,
            'entity_id': entity_id,
            'class_code': 'unknown',
            'label': entity_id,
            'bioentity_type': 'unknown',
            'description': '',
            'synonym': '',
            'source': 'unknown',
            'source_type': 'unknown',
            'integer_id': None
        }
  
# --- ESTRAZIONE E CONVERSIONE ---
#%%
# Caricamento dei dati
print(f"Loading data from {pkt_file}...")
pkt_kg = read_tar_rdf(pkt_file)
print(f"Loaded {len(pkt_kg)} RDF triples.")
# Nota: La riga di campionamento è commentata, se i dati sono grandi potresti volerla decommentare
# pkt_kg.sample(10000).to_csv(r'..\temp_dir\sample_10000.csv', index=False)
# pkt_kg = pd.read_csv(r'..\temp_dir\sample_10000.csv')
#%%


# if f"{output_dir}metadata_lookup.json" exists, then load it back
metadata_lookup_path = f"../temp_dir/metadata_lookup.json"
if os.path.exists(metadata_lookup_path):
    print(f"Loading existing metadata lookup from {metadata_lookup_path}...")
    # Load existing metadata lookup from the zip file wothout decompressing
    with zipfile.ZipFile(f"../temp_dir/metadata_lookup.zip", 'r') as zipf:
        with zipf.open('metadata_lookup.json') as f:
            metadata_lookup = json.load(f)
    print(f"Loaded existing metadata lookup from {metadata_lookup_path}")
    
else:
    print(f"Loading NodeLabels from {node_label_file}...")
    NodeLabels = pd.read_csv(node_label_file)
    print(f"Loaded {len(NodeLabels)} NodeLabels.")
    import zipfile

    # Create metadata lookup from NodeLabels
    metadata_lookup = create_metadata_lookup(NodeLabels)
    # save metadatalookup as json in temp dir 
    with open(f"{output_dir}metadata_lookup.json", 'w', encoding='utf-8') as f:
        json.dump(metadata_lookup, f, indent=2, ensure_ascii=False)
        # compess the json file in zip
        print(f"Saved metadata lookup to {output_dir}metadata_lookup.json")
        with zipfile.ZipFile(f"{output_dir}metadata_lookup.zip", 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(f"{output_dir}metadata_lookup.json", arcname="metadata_lookup.json")
        # delete the uncompressed json file
        os.remove(f"{output_dir}metadata_lookup.json")
        print(f"Compressed metadata lookup to {output_dir}metadata_lookup.zip")

#%%
# Funzione per convertire RDF in formato Grafo di Proprietà
def convert_rdf_to_property_graph(df, metadata_lookup):
    """Convert RDF dataframe to property graph with nodes and edges using metadata"""
    
    # Estrarre nodi unici (soggetti e oggetti)
    nodes = {}
    edges = []
    
    print("Converting RDF triples to property graph...")
    
    for idx, row in df.iterrows():
        subject = row['subject']
        predicate = row['predicate']
        obj = row['object']
        
        # 1. Processare il nodo soggetto
        if subject not in nodes:
            subject_meta = get_entity_metadata(subject, metadata_lookup)
            # ArangoDB richiede un _key per l'importazione
            safe_key = subject_meta['entity_id'].replace(':', '_').replace('/', '_').replace('?', '_').replace('=', '_').replace('#', '_').replace('@', '_').replace('.', '_')
            nodes[subject] = {
                '_key': safe_key,
                'uri': subject_meta['uri'],
                'namespace': subject_meta['namespace'],
                'entity_id': subject_meta['entity_id'],
                'class_code': subject_meta['class_code'],
                'label': subject_meta['label'],
                'bioentity_type': subject_meta['bioentity_type'],
                'description': subject_meta['description'],
                'synonym': subject_meta['synonym'],
                'source': subject_meta['source'],
                'source_type': subject_meta['source_type'],
                'integer_id': subject_meta['integer_id']
            }
        
        # 2. Processare il nodo oggetto
        if obj not in nodes:
            object_meta = get_entity_metadata(obj, metadata_lookup)
            # ArangoDB richiede un _key per l'importazione
            safe_key = object_meta['entity_id'].replace(':', '_').replace('/', '_').replace('?', '_').replace('=', '_').replace('#', '_').replace('@', '_').replace('.', '_')
            nodes[obj] = {
                '_key': safe_key,
                'uri': object_meta['uri'],
                'namespace': object_meta['namespace'],
                'entity_id': object_meta['entity_id'],
                'class_code': object_meta['class_code'],
                'label': object_meta['label'],
                'bioentity_type': object_meta['bioentity_type'],
                'description': object_meta['description'],
                'synonym': object_meta['synonym'],
                'source': object_meta['source'],
                'source_type': object_meta['source_type'],
                'integer_id': object_meta['integer_id']
            }
        
        # 3. Creare l'arco - get metadata for predicate too
        predicate_meta = get_entity_metadata(predicate, metadata_lookup)
        edge = {
            # ArangoDB crea automaticamente _key per gli archi se non fornito
            'edge_id': f"edge_{idx}",
            'source_uri': subject.strip('<>'),
            'target_uri': obj.strip('<>'),
            'predicate_uri': predicate_meta['uri'],
            'predicate_label': predicate_meta['label'],
            'predicate_class_code': predicate_meta['class_code'],
            'predicate_bioentity_type': predicate_meta['bioentity_type'],
            'predicate_source': predicate_meta['source'],
            # Questi campi verranno popolati in ArangoDB_Import.py per gli archi
        }
        edges.append(edge)
        
        if idx % 10000 == 0 and idx > 0:
            print(f"Processed {idx} triples...")
    
    return list(nodes.values()), edges

# Esecuzione della conversione
if USE_SAMPLE:
    print("Using sample of 10,000 triples for conversion...")
    sample_df = pkt_kg.sample(10000, random_state=42)
    nodes, edges = convert_rdf_to_property_graph(sample_df, metadata_lookup)
else:
    nodes, edges = convert_rdf_to_property_graph(pkt_kg, metadata_lookup)

# Calcolo delle statistiche
node_types = defaultdict(int)
edge_types = defaultdict(int)

for node in nodes:
    node_types[node['bioentity_type']] += 1

for edge in edges:
    edge_types[edge['predicate_label']] += 1

print(f"\nConversion completed!")
print(f"Total nodes: {len(nodes)}")
print(f"Total edges: {len(edges)}")

print("\nNode type distribution (by bioentity_type):")
for node_type, count in sorted(node_types.items(), key=lambda x: x[1], reverse=True)[:20]:
    print(f"  {node_type}: {count}")

print("\nEdge type distribution (by predicate_label):")
for edge_type, count in sorted(edge_types.items(), key=lambda x: x[1], reverse=True)[:20]:
    print(f"  {edge_type}: {count}")

# --- ESPORTAZIONE JSON PER ARANGODB ---

def export_to_json_for_arangodb(nodes, edges, output_dir):
    """Export property graph to JSON format compatible with ArangoDB's import utility"""
    
    os.makedirs(output_dir, exist_ok=True)

    # 1. Esportare nodi separatamente (già con _key)
    with open(f"{output_dir}nodes.json", 'w', encoding='utf-8') as f:
        # ArangoDB può importare un array di JSON. dump è sufficiente.
        json.dump(nodes, f, indent=2, ensure_ascii=False) 
    
    # 2. Esportare archi separatamente (senza _from e _to, saranno aggiunti al momento dell'importazione)
    # Nota: Stiamo esportando gli archi così come sono, i campi _from e _to verranno aggiunti nel passaggio 2.
    with open(f"{output_dir}edges.json", 'w', encoding='utf-8') as f:
        json.dump(edges, f, indent=2, ensure_ascii=False)
    
    print(f"\nJSON exports for ArangoDB completed in {output_dir}")
    print(f"Files: nodes.json (Nodes), edges.json (Edges)")

# Esecuzione dell'esportazione
export_to_json_for_arangodb(nodes, edges, output_dir)

print("\n" + "="*50)
print("PROPERTY GRAPH CONVERSION COMPLETED (STEP 1)")
print("="*50)