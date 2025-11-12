#%%
import pandas as pd
import json
import re
from collections import defaultdict
from urllib.parse import urlparse
import sys, os
import xml.etree.ElementTree as ET

from  pkt_utils import *
set_working_directory()
print_file_contents()
#%%

# --- CONFIGURAZIONE ---
pkt_file = "PKT.nt.tar.gz"
node_label_file = "PKT_NodeLabels_with_metadata_v3.0.2.csv"
output_dir = "../temp_dir/" # Assicurati che questa directory esista!
os.makedirs(output_dir, exist_ok=True)

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

# --- LOGICA DI MAPPATURA (mantenuta) ---

# non è più necessario creare i "mappings" perchè ora i metadati sono gia contentuto in NodeLabels_with_metadata_v3.0.2.csv (il preprocessing è stato gia fatto)
# sistrivi di conseguenza e accuratamente le funziono di per construire i JSON per il Property Graph di ArangoDB, i nomi dei nodi e degli archi (uri) e le loro proprietà prese dal file NodeLabels_with_metadata_v3.0.2.csv


def load_ro_mappings_from_owl(owl_file_path="../data/owl/ro.owl"):
    """Load RO relation mappings from ro.owl file (mantenuto)"""
    # [Contenuto della funzione load_ro_mappings_from_owl]
    # ... (Il corpo completo della funzione load_ro_mappings_from_owl del file originale)
    try:
        # Parsing del file OWL (logica di fallback se il parsing fallisce)
        tree = ET.parse(owl_file_path)
        root = tree.getroot()
        
        namespaces = {
            'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
            'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
            'owl': 'http://www.w3.org/2002/07/owl#',
            'obo': 'http://purl.obolibrary.org/obo/'
        }
        
        relation_mappings = {}
        for obj_prop in root.findall('.//owl:ObjectProperty', namespaces):
            about = obj_prop.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about')
            if about and 'RO_' in about:
                ro_code = about.split('/')[-1]
                label_elem = obj_prop.find('.//rdfs:label', namespaces)
                if label_elem is not None and label_elem.text:
                    clean_label = label_elem.text.lower().replace(' ', '_').replace('-', '_')
                    relation_mappings[ro_code] = clean_label
        
        print(f"Loaded {len(relation_mappings)} RO mappings from {owl_file_path}")
        return relation_mappings
        
    except Exception as e:
        print(f"Error loading RO mappings from OWL file: {e}. Using hardcoded fallback.")
        return {
            'RO_0002436': 'molecularly_interacts_with',
            'RO_0002434': 'interacts_with',
            'RO_0002511': 'inheres_in',
            'RO_0002212': 'negatively_regulates',
            'RO_0002213': 'positively_regulates',
            'RO_0002211': 'regulates',
            'RO_0002331': 'involved_in',
            'RO_0002202': 'develops_from',
            'RO_0002220': 'adjacent_to'
        }

def load_bfo_mappings_from_owl(owl_file_path="../data/owl/bfo_classes_only.owl"):
    """Load BFO class mappings from bfo_classes_only.owl file (mantenuto)"""
    # [Contenuto della funzione load_bfo_mappings_from_owl]
    # ... (Il corpo completo della funzione load_bfo_mappings_from_owl del file originale)
    try:
        # Parsing del file OWL (logica di fallback se il parsing fallisce)
        tree = ET.parse(owl_file_path)
        root = tree.getroot()
        
        namespaces = {
            'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
            'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
            'owl': 'http://www.w3.org/2002/07/owl#',
            'obo': 'http://purl.obolibrary.org/obo/'
        }
        
        class_mappings = {}
        for owl_class in root.findall('.//owl:Class', namespaces):
            about = owl_class.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about')
            if about and 'BFO_' in about:
                bfo_code = about.split('/')[-1]
                label_elem = owl_class.find('.//obo:BFO_0000179', namespaces) or owl_class.find('.//rdfs:label', namespaces)
                if label_elem is not None and label_elem.text:
                    clean_label = label_elem.text.lower().replace(' ', '_').replace('-', '_')
                    class_mappings[bfo_code] = clean_label
        
        print(f"Loaded {len(class_mappings)} BFO class mappings from {owl_file_path}")
        return class_mappings
        
    except Exception as e:
        print(f"Error loading BFO mappings from OWL file: {e}. Using hardcoded fallback.")
        return {
            'BFO_0000001': 'entity', 'BFO_0000002': 'continuant', 'BFO_0000003': 'occurrent',
            'BFO_0000004': 'independent_continuant', 'BFO_0000006': 'spatial_region',
            'BFO_0000008': 'temporal_region', 'BFO_0000015': 'process',
            'BFO_0000016': 'disposition', 'BFO_0000017': 'realizable_entity',
            'BFO_0000019': 'quality', 'BFO_0000020': 'specifically_dependent_continuant',
            'BFO_0000023': 'role', 'BFO_0000024': 'fiat_object_part',
            'BFO_0000027': 'object_aggregate', 'BFO_0000029': 'site',
            'BFO_0000030': 'object', 'BFO_0000031': 'generically_dependent_continuant',
            'BFO_0000034': 'function', 'BFO_0000035': 'process_boundary',
            'BFO_0000038': 'one_dimensional_temporal_region', 'BFO_0000040': 'material_entity',
            'BFO_0000141': 'immaterial_entity'
        }

def extract_entity_info(uri):
    """Extract namespace, entity type, and ID from URI (mantenuto)"""
    parsed = urlparse(uri)
    namespace = parsed.netloc if parsed.netloc else "unknown"
    path = parsed.path
    entity_id = uri.split('/')[-1]
    
    if not hasattr(extract_entity_info, '_bfo_mappings'):
        extract_entity_info._bfo_mappings = load_bfo_mappings_from_owl()
    
    bfo_mappings = extract_entity_info._bfo_mappings
    
    entity_type = "unknown"
    if "obo/PR_" in uri:
        entity_type = "protein"
    elif "obo/CHEBI_" in uri:
        entity_type = "chemical"
    elif "obo/RO_" in uri:
        entity_type = "relation"
    elif "obo/BFO_" in uri:
        bfo_code = entity_id
        entity_type = bfo_mappings.get(bfo_code, f"bfo_{bfo_code}")
    elif "gene/" in uri:
        entity_type = "gene"
    elif "ensembl.org" in uri:
        entity_type = "transcript"
    elif "obo/GO_" in uri:
        entity_type = "gene_ontology"
    elif "obo/HP_" in uri:
        entity_type = "phenotype"
    elif "obo/CL_" in uri:
        entity_type = "cell_type"
    elif "obo/UBERON_" in uri:
        entity_type = "anatomy"
    
    return {
        'namespace': namespace,
        'entity_type': entity_type,
        'entity_id': entity_id,
        'full_uri': uri
    }

def extract_predicate_info(predicate_uri):
    """Extract relation type from predicate URI using complete RO ontology (mantenuto)"""
    if not hasattr(extract_predicate_info, '_ro_mappings'):
        extract_predicate_info._ro_mappings = load_ro_mappings_from_owl()
    
    relation_mappings = extract_predicate_info._ro_mappings
    
    ro_match = re.search(r'RO_(\d+)', predicate_uri)
    if ro_match:
        ro_code = f"RO_{ro_match.group(1)}"
        return relation_mappings.get(ro_code, ro_code)
    
    return predicate_uri.split('/')[-1]

# --- ESTRAZIONE E CONVERSIONE ---

# Caricamento dei dati
print(f"Loading data from {pkt_file}...")
pkt_kg = read_tar_rdf_csv(pkt_file)
print(f"Loaded {len(pkt_kg)} RDF triples.")
# Nota: La riga di campionamento è commentata, se i dati sono grandi potresti volerla decommentare
# pkt_kg.sample(10000).to_csv(r'..\temp_dir\sample_10000.csv', index=False)
# pkt_kg = pd.read_csv(r'..\temp_dir\sample_10000.csv')
#%%
NodeLabels = pd.read_csv(node_label_file)
print(f"Loaded {len(NodeLabels)} NodeLabels.")
#%%

#%%
# Funzione per convertire RDF in formato Grafo di Proprietà
def convert_rdf_to_property_graph(df):
    """Convert RDF dataframe to property graph with nodes and edges"""
    
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
            subject_info = extract_entity_info(subject)
            # ArangoDB richiede un _key per l'importazione
            safe_key = subject_info['entity_id'].replace(':', '_').replace('/', '_').replace('?', '_').replace('=', '_').replace('#', '_').replace('@', '_')
            nodes[subject] = {
                '_key': safe_key,
                'node_id': subject,
                'label': subject_info['entity_type'],
                'namespace': subject_info['namespace'],
                'entity_id': subject_info['entity_id'],
                'full_uri': subject_info['full_uri']
            }
        
        # 2. Processare il nodo oggetto
        if obj not in nodes:
            object_info = extract_entity_info(obj)
            # ArangoDB richiede un _key per l'importazione
            safe_key = object_info['entity_id'].replace(':', '_').replace('/', '_').replace('?', '_').replace('=', '_').replace('#', '_').replace('@', '_')
            nodes[obj] = {
                '_key': safe_key,
                'node_id': obj,
                'label': object_info['entity_type'],
                'namespace': object_info['namespace'],
                'entity_id': object_info['entity_id'],
                'full_uri': object_info['full_uri']
            }
        
        # 3. Creare l'arco
        edge = {
            # ArangoDB crea automaticamente _key per gli archi se non fornito
            'edge_id': f"edge_{idx}",
            'source_uri': subject,
            'target_uri': obj,
            'relationship': extract_predicate_info(predicate),
            'predicate_uri': predicate,
            # Questi campi verranno popolati in ArangoDB_Import.py per gli archi
        }
        edges.append(edge)
        
        if idx % 10000 == 0 and idx > 0:
            print(f"Processed {idx} triples...")
    
    return list(nodes.values()), edges

# Esecuzione della conversione
nodes, edges = convert_rdf_to_property_graph(pkt_kg)

# Calcolo delle statistiche (mantenuto)
node_types = defaultdict(int)
edge_types = defaultdict(int)

for node in nodes:
    node_types[node['label']] += 1

for edge in edges:
    edge_types[edge['relationship']] += 1

print(f"\nConversion completed!")
print(f"Total nodes: {len(nodes)}")
print(f"Total edges: {len(edges)}")

print("\nNode type distribution:")
for node_type, count in sorted(node_types.items()):
    print(f"  {node_type}: {count}")

print("\nEdge type distribution:")
for edge_type, count in sorted(edge_types.items()):
    print(f"  {edge_type}: {count}")

# --- ESPORTAZIONE JSON PER ARANGODB ---

def export_to_json_for_arangodb(nodes, edges, output_dir):
    """Export property graph to JSON format compatible with ArangoDB's import utility"""
    
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