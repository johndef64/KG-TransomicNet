
#%%
import os
from  pkt_utils import *
set_working_directory()

file_root = "PheKnowLator_v3.0.2_full_instance_inverseRelations_OWLNETS_INSTANCE_purified"
file_root = "PKT"

kg_files = {
# other files
# 1:f"{file_root}_NetworkxMultiDiGraph.gpickle.tar.gz",
# 2:f"{file_root}_decoding_dict.pkl.tar.gz",

# ntriples
3:f"{file_root}.nt.tar.gz",

# tables
5:f"{file_root}_NodeLabels.txt.tar.gz",
6:f"{file_root}_Triples_Identifiers.txt.tar.gz",
7:f"{file_root}_Triples_Integers.txt.tar.gz",

}

import os
# from networkx import display
import pandas as pd

def read_kg_file(file_id, base_path=""):
    file_name = kg_files[file_id]
    full_path = os.path.join(base_path, file_name)
    
    print(f"Reading file ID {file_id}: {full_path}")
    
    if file_id == 3:  # NTriples
        import tarfile
        with tarfile.open(full_path, 'r:gz') as tar:
            members = tar.getmembers()
            extracted_file = tar.extractfile(members[0])
            df = pd.read_csv(extracted_file, sep=' ', header=None, names=['subject', 'predicate', 'object', '.'])
            return df[['subject', 'predicate', 'object']]
    
    elif file_id in [5, 6, 7]:  # Tables
        import tarfile
        with tarfile.open(full_path, 'r:gz') as tar:
            members = tar.getmembers()
            extracted_file = tar.extractfile(members[0])
            df = pd.read_csv(extracted_file, sep='\t')
            return df
    
    else:
        raise ValueError("Unsupported file ID")

file6 = read_kg_file(6)  # Example usage to read Triples_Identifiers
file6.head()
#%%
file6.describe().to_clipboard()
"""
	subject	predicate	object
count	11132839	11132839	11132839
unique	780752	296	532497
top	<http://purl.obolibrary.org/obo/UBERON_0000473>	<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>	<http://purl.obolibrary.org/obo/SO_0000673>
freq	21869	4550215	191121

"""
#%%
file3 = read_kg_file(3)  # Example usage to read NTriples
file3.head()
#%%
desc3 = file3.describe()
display(desc3)
desc3.to_clipboard()
"""
	subject	predicate	object
count	11132839	11132839	11132839
unique	780752	296	532497
top	<http://purl.obolibrary.org/obo/UBERON_0000473>	<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>	<http://purl.obolibrary.org/obo/SO_0000673>
freq	21869	4550215	191121

"""
#%%
NodeLables = read_kg_file(5)  # Example usage to read NodeLabels
NodeLables.head()
NodeLables.integer_id = NodeLables.integer_id.astype('str')  # convert to string for description
desc5 = NodeLables.describe(include='object') # describe only objects
display(desc5)
desc5.to_clipboard()

"""
	entity_type	integer_id	entity_uri	label	description/definition	synonym
count	780025	780119	780119	780025	606751	288251
unique	2	780119	780119	778810	605205	287910
top	NODES	763062	<https://www.ncbi.nlm.nih.gov/snp/rs368043574>	aroA from Shigella	A immortal human B cell line cell that has the characteristics: Human B lymphocytes transformed by Epstein-Barr Virus. This B cell lines derived from Mongoloid Minority Groups in South America.	E3
freq	779732	1	1	4	147	8
"""
#%%

NodeLables.head().T.to_clipboard()
"""
	entity_type	integer_id	entity_uri	label	description/definition	synonym
0	NODES	536450	<https://www.ncbi.nlm.nih.gov/snp/rs876660433>	NM_000535.7(PMS2):c.567T>C (p.His189=)	This variant is a germline single nucleotide variant located on chromosome 7 (NC_000007.14, start:5999246/stop:5999246 positions, cytogenetic location:7p22.1) and has clinical significance 'Likely benign'. This entry is for the GRCh38 and was last reviewed on May 03, 2020 with review status 'criteria provided, multiple submitters, no conflicts'.	
1	NODES	338642	<https://www.ncbi.nlm.nih.gov/snp/rs35349917>	NM_001105206.3(LAMA4):c.280G>A (p.Gly94Ser)	This variant is a germline/unknown single nucleotide variant located on chromosome 6 (NC_000006.12, start:112216385/stop:112216385 positions, cytogenetic location:6q21) and has clinical significance 'Benign'. This entry is for the GRCh38 and was last reviewed on Dec 07, 2020 with review status 'criteria provided, multiple submitters, no conflicts'.	
2	NODES	253428	<http://purl.obolibrary.org/obo/CHEBI_177646>	3,4,15-Trihydroxy-12,13-epoxytrichothec-9-en-8-yl 3-methylbutanoate		[(1S,2R,4S,7R,9R,10R,11S)-10,11-dihydroxy-2-(hydroxymethyl)-1,5-dimethylspiro[8-oxatricyclo[7.2.1.02,7]dodec-5-ene-12,2'-oxirane]-4-yl] 3-methylbutanoate
3	NODES	430342	<http://purl.obolibrary.org/obo/PR_000031649>	Ly6/PLAUR domain-containing protein 6	A protein that is a translation product of the human LYPD6 gene or a 1:1 ortholog thereof.	LYPD6|UNQ3023/PRO9821
4	NODES	506783	<https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000522820>	NECAB1-205	Transcript NECAB1-205 is classified as type 'protein_coding'.	


	0	1	2	3	4
entity_type	NODES	NODES	NODES	NODES	NODES
integer_id	536450	338642	253428	430342	506783
entity_uri	<https://www.ncbi.nlm.nih.gov/snp/rs876660433>	<https://www.ncbi.nlm.nih.gov/snp/rs35349917>	<http://purl.obolibrary.org/obo/CHEBI_177646>	<http://purl.obolibrary.org/obo/PR_000031649>	<https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000522820>
label	NM_000535.7(PMS2):c.567T>C (p.His189=)	NM_001105206.3(LAMA4):c.280G>A (p.Gly94Ser)	3,4,15-Trihydroxy-12,13-epoxytrichothec-9-en-8-yl 3-methylbutanoate	Ly6/PLAUR domain-containing protein 6	NECAB1-205
description/definition	This variant is a germline single nucleotide variant located on chromosome 7 (NC_000007.14, start:5999246/stop:5999246 positions, cytogenetic location:7p22.1) and has clinical significance 'Likely benign'. This entry is for the GRCh38 and was last reviewed on May 03, 2020 with review status 'criteria provided, multiple submitters, no conflicts'.	This variant is a germline/unknown single nucleotide variant located on chromosome 6 (NC_000006.12, start:112216385/stop:112216385 positions, cytogenetic location:6q21) and has clinical significance 'Benign'. This entry is for the GRCh38 and was last reviewed on Dec 07, 2020 with review status 'criteria provided, multiple submitters, no conflicts'.		A protein that is a translation product of the human LYPD6 gene or a 1:1 ortholog thereof.	Transcript NECAB1-205 is classified as type 'protein_coding'.
synonym			[(1S,2R,4S,7R,9R,10R,11S)-10,11-dihydroxy-2-(hydroxymethyl)-1,5-dimethylspiro[8-oxatricyclo[7.2.1.02,7]dodec-5-ene-12,2'-oxirane]-4-yl] 3-methylbutanoate	LYPD6|UNQ3023/PRO9821	

"""
# entityuri_label_map = dict(zip(NodeLables.entity_uri, NodeLables.label))

# create new column "entity" parsing "entity_uri"
# usle last string from last "/" or "=" if present
# remove trailing ">" if present
NodeLables['entity'] = NodeLables['entity_uri'].apply(lambda x: x.rstrip('>').split('/')[-1].split('=')[-1])
NodeLables[['entity_uri', 'entity', "label"]]#.head()#.to_clipboard()

#%%
# search PKT documentation for all entity types

# crea column per entity_class
# removing every numers from right side of "entity" column
import re
def extract_entity_class(entity):
    s = re.sub(r'(?<!\d)\d+$', '', entity).rstrip('_-')
    for prefix in ('PR_', 'CHR_', "GNO_"):
        if s.startswith(prefix):
            return prefix.rstrip('_')
    return s.split('#')[0]
NodeLables['class_code'] = NodeLables['entity'].apply(extract_entity_class)

NodeLables[['entity', 'class_code', 'label']]
#%%
NodeLables['class_code'].value_counts()

#%%
NodeLables[['entity', 'class_code', 'label']]. drop_duplicates(subset=['class_code']).sort_values('class_code').to_clipboard(index=False)



#%%
# NodeLables[NodeLables['entity_type']!='NODES']
NodeLables[NodeLables['entity_uri'].str.contains("#")]

#%% ==================================================
# node_metadata_dict is simply == NodeLabels 

with open("node_metadata_dict.pkl", "rb") as f:
    node_metadata_dict = pickle.load(f)
#%%
type(node_metadata_dict)
len(node_metadata_dict)
node_metadata_dict.keys()
node_metadata_dict["relations"]
# list(node_metadata_dict.items())[:5]

#%% ==================================================
#load and read Master_Edge_List_Dict.json
import json
with open("Master_Edge_List_Dict.json", "r") as f:
    master_edge_list_dict = json.load(f)
#%%
type(master_edge_list_dict)
len(master_edge_list_dict)
master_edge_list_dict.keys()

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


"""
protein
chemical
disease
gene
pathway
rna
variant
anatomy
phenotype
gomf
gocc
gobp
cell
catalyst
cofactor
"""

#%%
edge_type = "variant-gene"
edge_type = "variant-phenotype"
for edge_type in edge_list_dict:
    print(f"Number of edges in {edge_type}: {len(master_edge_list_dict[edge_type]['edge_list'])}")
    # print 5 sample edges
    print("Sample edges:")
    for edge in master_edge_list_dict[edge_type]['edge_list'][:5]:
        print(edge)
# master_edge_list_dict[edge_type]

#%%

#%%
# sostitisci "UMLS_C" in class_code con "UMLS"
NodeLables['class_code'] = NodeLables['class_code'].replace('UMLS_C', 'UMLS')
# sostituisci "22-rdf-syntax-ns" in class_code con "RDF"
NodeLables['class_code'] = NodeLables['class_code'].replace('22-rdf-syntax-ns', 'RDF')
# elimina la colonna "entity_class"
NodeLables = NodeLables.drop(columns=['entity_class'], errors='ignore')

#%%
NodeLables['class_code'] = NodeLables['class_code'].replace('', 'EntrezID')
NodeLables['class_code'] = NodeLables['class_code'].replace('rs', 'dbSNP')
NodeLables.to_csv("PKT_NodeLabels_with_class_code_v3.0.2.csv", index=False)
#%%
NodeLables

#%%
file7 = read_kg_file(7)  # Example usage to read Triples_Integers
file7.head()


#%%

NodeLables[NodeLables.entity_uri.str.contains("MF_")]