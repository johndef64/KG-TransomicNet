#%%
"""
Utility functions for querying omics data from ArangoDB.
Supports: Gene Expression, CNV, Methylation, miRNA, Protein data.

UNIFIED COLLECTION STRUCTURE (v3):
Each omic collection (GENE_EXPRESSION, CNV, METHYLATION, MIRNA, PROTEIN) contains
two types of documents distinguished by the `data_type` field:
- *_index: Index document mapping features to positions
- *_vector: Sample document with quantitative vectors
"""

from typing import Optional, List, Dict, Any, Union

# =============================================================================
# COLLECTION CONFIGURATION
# =============================================================================
# Mapping of omic types to their unified collection and field names
OMIC_CONFIG = {
    "gene_expression": {
        "collection": "GENE_EXPRESSION",
        "index_data_type": "gene_expression_index",
        "vector_data_type": "gene_expression_vector",
        "index_key_prefix": "expr_index",
        "mapping_field": "gene_mappings",
        "value_fields": ["values_tpm", "values_fpkm", "values_counts"]
    },
    "cnv": {
        "collection": "CNV",
        "index_data_type": "cnv_index",
        "vector_data_type": "cnv_vector",
        "index_key_prefix": "cnv_index",
        "mapping_field": "gene_mappings",
        "value_fields": ["values_copy_number"]
    },
    "methylation": {
        "collection": "METHYLATION",
        "index_data_type": "methylation_index",
        "vector_data_type": "methylation_vector",
        "index_key_prefix": "methylation_index",
        "mapping_field": "probe_mappings",
        "value_fields": ["values_beta"]
    },
    "mirna": {
        "collection": "MIRNA",
        "index_data_type": "mirna_index",
        "vector_data_type": "mirna_vector",
        "index_key_prefix": "mirna_index",
        "mapping_field": "mirna_mappings",
        "value_fields": ["values_expression"]
    },
    "protein": {
        "collection": "PROTEIN",
        "index_data_type": "protein_index",
        "vector_data_type": "protein_vector",
        "index_key_prefix": "protein_index",
        "mapping_field": "protein_mappings",
        "value_fields": ["values_abundance"]
    }
}


KG_CONFING = {
 "nodes": {
    "key": "_key",
    "uri": "uri",
    "label": "label",
    "bioentity_type": "bioentity_type",
    "class_code": "class_code"
 },
 "edges": {
    "key": "_key",
    "uri": "uri",
    "label": "label",
    "predicate_label": "predicate_label",
    "source": "_from",
    "target": "_to"
 }
}


# =============================================================================
# KNOWLEDGE GRAPH QUERIES
# =============================================================================

# basic KG functions

def get_key_value_counts(db_connection, collection_name: str, key: str, limit: int = None, verbose = True) -> Dict[str, int]:
    """
    Get all unique values for a key in a collection with occurrence counts.
    
    Args:
        db_connection: ArangoDB database connection
        collection_name: Name of the collection to query (e.g., "nodes", "edges")
        key: The document key/field to aggregate (e.g., "bioentity_type", "predicate_label")
        limit: Optional limit on number of distinct values to return
    
    Returns:
        Dict mapping each unique value to its count, sorted by count descending
    """
    limit_clause = f"LIMIT {limit}" if limit else ""
    
    aql = f"""
    FOR doc IN {collection_name}
        FILTER doc.@key != null
        COLLECT value = doc.@key WITH COUNT INTO count
        SORT count DESC
        {limit_clause}
        RETURN {{value: value, count: count}}
    """
    
    bind_vars = {"key": key}
    results = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    if verbose:
        print(f"Unique values for '{key}' in '{collection_name}':")
        for r in results:
            print(f"  {r['value']}: {r['count']}")

    return {r["value"]: r["count"] for r in results}


def get_node_by_key(db_connection, node_key: str) -> Optional[Dict]:
    """
    Get a node document by its _key.
    
    Args:
        db_connection: ArangoDB database connection
        node_key: The _key of the node (e.g., "80724")
    
    Returns:
        Node document or None if not found
    """
    aql = """
    FOR n IN nodes
        FILTER n._key == @nodeKey
        RETURN n
    """
    result = list(db_connection.aql.execute(aql, bind_vars={"nodeKey": node_key}))
    return result[0] if result else None


def get_node_by_uri(db_connection, uri: str) -> Optional[Dict]:
    """
    Get a node document by its URI.
    
    Args:
        db_connection: ArangoDB database connection
        uri: The URI of the node
    
    Returns:
        Node document or None if not found
    """
    aql = """
    FOR n IN nodes
        FILTER n.uri == @uri
        RETURN n
    """
    result = list(db_connection.aql.execute(aql, bind_vars={"uri": uri}))
    return result[0] if result else None

def get_node_by_chosen_field(db_connection, field_name: str, field_value: str, 
                              limit: int = 100, verbose: bool = True) -> Union[Optional[Dict], List[Dict]]:
    """
    Get a node document by a chosen field and value.
    
    Args:
        db_connection: ArangoDB database connection
        field_name: The name of the field to filter by (e.g., "label", "class_code", "bioentity_type")
        field_value: The value to match in the specified field (supports '*' wildcard for regex)
        limit: Maximum number of results to return (for regex searches)
        verbose: Whether to print debug information
    
    Returns:
        Node document or None if not found (exact match)
        List of matching nodes if wildcard pattern is used
    
    Examples:
        get_node_by_chosen_field(db, "label", "TP53 (human)")     # Exact match
        get_node_by_chosen_field(db, "label", "BRCA*")            # All labels starting with BRCA
        get_node_by_chosen_field(db, "class_code", "EntrezID")    # All EntrezID nodes
        get_node_by_chosen_field(db, "bioentity_type", "*RNA*")   # Any bioentity containing RNA
    """
    # Check if field_value contains wildcard for regex search
    use_regex = '*' in field_value
    
    if use_regex:
        # Convert wildcard pattern to regex: '*' -> '.*'
        regex_pattern = '^' + field_value.replace('*', '.*') + '$'
        
        aql = """
        FOR n IN nodes
            FILTER REGEX_TEST(n.@fieldName, @regexPattern, true)
            LIMIT @limit
            RETURN n
        """
        bind_vars = {
            "fieldName": field_name,
            "regexPattern": regex_pattern,
            "limit": limit
        }
        result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
        
        if verbose:
            print(f"Found {len(result)} nodes matching '{field_name}' ~ '{field_value}'")
            if result and len(result) <= 10:
                for r in result:
                    print(f"  - {r.get('_key')}: {r.get('label', 'N/A')}")
            elif result:
                print(f"  (showing first 10 of {len(result)})")
                for r in result[:10]:
                    print(f"  - {r.get('_key')}: {r.get('label', 'N/A')}")
        
        return result  # Return list for regex matches
    else:
        # Exact match query
        aql = """
        FOR n IN nodes
            FILTER n.@fieldName == @fieldValue
            LIMIT 1
            RETURN n
        """
        bind_vars = {
            "fieldName": field_name,
            "fieldValue": field_value
        }
        result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
        
        if verbose and result:
            print(f"Found node: {result[0].get('_key')} - {result[0].get('label', 'N/A')}")
        elif verbose:
            print(f"No node found with {field_name} = '{field_value}'")
        
        return result[0] if result else None

def get_edges_by_predicate(db_connection, predicate_label: str, limit: int = 100) -> List[Dict]:
    """
    Get edges filtered by predicate label.
    
    Args:
        db_connection: ArangoDB database connection
        predicate_label: The predicate label to filter by (e.g., "has gene product")
        limit: Maximum number of edges to return
    
    Returns:
        List of edge documents
    """
    aql = """
    FOR e IN edges
        FILTER e.predicate_label == @predicateLabel
        LIMIT @limit
        RETURN e
    """
    bind_vars = {"predicateLabel": predicate_label, "limit": limit}
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))

def get_node_neighbors(db_connection, node_key: str, direction: str = "any", limit: int = 100) -> List[Dict]:
    """
    Get neighboring nodes connected to a given node.
    
    Args:
        db_connection: ArangoDB database connection
        node_key: The _key of the source node
        direction: "outbound", "inbound", or "any"
        limit: Maximum number of neighbors to return
    
    Returns:
        List of dicts with edge and neighbor node info
    """
    direction_map = {"outbound": "OUTBOUND", "inbound": "INBOUND", "any": "ANY"}
    dir_clause = direction_map.get(direction, "ANY")
    
    aql = f"""
    FOR v, e IN 1..1 {dir_clause} CONCAT("nodes/", @nodeKey) edges
        LIMIT @limit
        RETURN {{
            edge: e,
            neighbor: v
        }}
    """
    bind_vars = {"nodeKey": node_key, "limit": limit}
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))

def search_nodes(db_connection, bioentity_type: str = None, 
                 class_code: str = None, 
                 label_contains: str = None,
                 limit: int = 100) -> List[Dict]:
    """
    Search nodes by various criteria.
    
    Args:
        db_connection: ArangoDB database connection
        bioentity_type: Filter by bioentity type (e.g., "gene", "protein")
        class_code: Filter by class code (e.g., "EntrezID", "PR")
        label_contains: Filter nodes whose label contains this substring
        limit: Maximum number of nodes to return
    
    Returns:
        List of matching node documents
    """
    filters = []
    bind_vars = {"limit": limit}
    
    if bioentity_type:
        filters.append("n.bioentity_type == @bioentityType")
        bind_vars["bioentityType"] = bioentity_type
    if class_code:
        filters.append("n.class_code == @classCode")
        bind_vars["classCode"] = class_code
    if label_contains:
        filters.append("CONTAINS(LOWER(n.label), LOWER(@labelContains))")
        bind_vars["labelContains"] = label_contains
    
    filter_clause = "FILTER " + " AND ".join(filters) if filters else ""
    
    aql = f"""
    FOR n IN nodes
        {filter_clause}
        LIMIT @limit
        RETURN n
    """
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


""" [KG COLLECTION SAMPLES]

"edge" collections:
_id:edges/edge_0
_rev:_kneF6ZO---
_key:edge_0
_from:nodes/80724
_to:nodes/PR_Q6JQN1
{
  "source_uri": "http://www.ncbi.nlm.nih.gov/gene/80724",
  "target_uri": "http://purl.obolibrary.org/obo/PR_Q6JQN1",
  "predicate_uri": "http://purl.obolibrary.org/obo/RO_0002205",
  "predicate_label": "has gene product",
  "predicate_class_code": "RO",
  "predicate_bioentity_type": "phenotype",
  "predicate_source": "Relation Ontology"
}

Unique values for 'predicate_class_code' in 'edges':
  RO: 6362632
  RDF: 4521549
  BFO: 119618
  chebi: 41661
  pr: 19848
  CLO: 11629
  VO: 2628
  so: 628
  BSPO: 522
  core: 473
  mondo: 369
  MONDO: 163
  pw: 136
  pato: 125
  cl: 38
  SIO: 38
  EFO: 33
  [...]


-------------------------------
"nodes" collections:
_id:nodes/80724
_rev:_kneFK2G---
_key:80724
{
  "uri": "http://www.ncbi.nlm.nih.gov/gene/80724",
  "namespace": "www.ncbi.nlm.nih.gov",
  "entity_id": "80724",
  "class_code": "EntrezID",
  "label": "ACAD10 (human)",
  "bioentity_type": "gene",
  "description": "A protein coding gene ACAD10 in human.",
  "synonym": "",
  "source": "NCBI Entrez Gene",
  "source_type": "Database",
  "integer_id": 1
}

Unique values for 'class_code' in 'nodes':
  ENST: 190613  
  CHEBI: 149902
  PR: 96085  
  GO: 43821
  CLO: 41738
  EntrezID: 26502 
  MONDO: 22305   
  HP: 16271
  UBERON: 14168
  R-HSA: 13504
  VO: 6240
  PW: 2598
  CL: 2365
  SO: 2362
  NCBITaxon: 2067
  unknown: 929
  DOID: 927
  GNO: 796
  PATO: 572
  ENVO: 397
  R-NUL: 278
  CHR: 249
  ECTO: 171
  NBO: 164
  FOODON: 162
  EFO: 144
  MOD: 91
  MPATH: 75
  MAXO: 72
  PO: 47
  OGG: 30
  NCIT_C: 27
  HsapDv: 16
  [...]
"""

if __name__ == "__main__":
    from arangodb_utils import *
    db_connection = setup_arangodb_connection("PKT_test10000")

    node_class_code = get_key_value_counts(db_connection, "nodes", "class_code")

    edge_class_code = get_key_value_counts(db_connection, "edges", "predicate_class_code")
#%%
edge_class_code


#%%
# =============================================================================
# GENE DATA QUERIES
# "GENE" collection is the pivot to map omics data to knowledge graph nodes. It contains comprehensive gene information and mappings to various identifiers.

# the mapping key is in the index of every omic collection, for example in GENE_EXPRESSION collection there is a document with data_type "gene_expression_index" that contains a field "gene_mappings" which is a list of dicts with keys like "entrez_id", "hgnc_symbol", "gene_id_base", "gene_id_ensembl" and "position". The position is the index of the gene in the expression vectors of the samples in the same collection with data_type "gene_expression_vector".

# to use HGCN symbol to retrive kg nodes remembver to remove "-" from the symbol, for example "MIRLET7-1" in hgnc_symbol becomes "MIRLET71" in the kg node label and in the GENES collection hgnc_symbol field.
# =============================================================================

def get_gene_info(db_connection, entrez_id: str = None, 
                  hgnc_symbol: str = None,
                  ensembl_id: str = None) -> Union[Optional[Dict], List[Dict]]:
    """
    Retrieve gene information by various identifiers.
    
    Args:
        db_connection: ArangoDB database connection
        entrez_id: Entrez gene ID
        hgnc_symbol: HGNC gene symbol (supports '*' wildcard for regex search)
        ensembl_id: Ensembl gene ID base
    
    Returns:
        Gene information dict or None if not found (exact match)
        List of matching genes if wildcard pattern is used
    
    Examples:
        get_gene_info(db, hgnc_symbol="TP53")        # Exact match
        get_gene_info(db, hgnc_symbol="MIRLET7*")    # All MIRLET7 family genes
        get_gene_info(db, hgnc_symbol="*MIR*")       # Any gene containing MIR
    """
    # Check if hgnc_symbol contains wildcard for regex search
    use_regex = hgnc_symbol is not None and '*' in hgnc_symbol
    
    if use_regex:
        # Convert wildcard pattern to regex: '*' -> '.*'
        # Anchor pattern for full match behavior
        regex_pattern = '^' + hgnc_symbol.replace('*', '.*') + '$'
        
        aql = """
        FOR g IN GENES
            FILTER REGEX_TEST(g.hgnc_symbol, @regexPattern, true)
            RETURN g
        """
        bind_vars = {"regexPattern": regex_pattern}
        result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
        return result  # Return list for regex matches
    else:
        # Exact match query
        aql = """
        FOR g IN GENES
            FILTER (@entrezId != null AND g.entrez_id == @entrezId)
                OR (@hgncSymbol != null AND g.hgnc_symbol == @hgncSymbol)
                OR (@ensemblId != null AND g.ensembl_id_base == @ensemblId)
            RETURN g
        """
        bind_vars = {
            "entrezId": entrez_id,
            "hgncSymbol": hgnc_symbol,
            "ensemblId": ensembl_id
        }
        result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
        return result[0] if result else None


def get_data_structure_info(db_connection, collection_name):
    """
    get dict atrigutes kyes and values  from on entry in collection_name
    """
    aql_query = f"""
    FOR doc IN {collection_name}
      LIMIT 1
      RETURN doc
    """
    result = db_connection.aql.execute(aql_query)
    for doc in result:
        return doc

def test_query(db_connection, query_aql):
    result = db_connection.aql.execute(query_aql)
    for doc in result:
        return doc

info = get_gene_info(db_connection, hgnc_symbol="MIRLET7*")
for i in info:
    print(i["hgnc_symbol"])
#%%
# =============================================================================
# GENE EXPRESSION QUERIES
# =============================================================================

def resolve_gene_position(db_connection, cohort: str, 
                          entrez_id: str = None, 
                          hgnc_symbol: str = None,
                          ensembl_id: str = None) -> Optional[int]:
    """
    Resolve gene identifier to position in expression vector.
    
    Args:
        db_connection: ArangoDB database connection
        cohort: TCGA cohort (e.g., "TCGA-BRCA")
        entrez_id: Entrez gene ID
        hgnc_symbol: HGNC gene symbol (e.g., "TP53")
        ensembl_id: Ensembl gene ID base (e.g., "ENSG00000141510")
    
    Returns:
        Position index in the expression vector, or None if not found
    """
    config = OMIC_CONFIG["gene_expression"]
    index_key = f"{config['index_key_prefix']}_{cohort}"
    
    aql = f"""
    LET idx = FIRST(
        FOR doc IN {config['collection']}
            FILTER doc.data_type == @indexDataType
            FILTER doc._key == @indexKey
            RETURN doc
    )
    LET pos = FIRST(
        FOR m IN idx.{config['mapping_field']}
            FILTER (@entrezId != null AND m.entrez_id == @entrezId)
                OR (@hgncSymbol != null AND m.hgnc_symbol == @hgncSymbol)
                OR (@ensemblId != null AND m.gene_id_base == @ensemblId)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
        "indexDataType": config['index_data_type'],
        "entrezId": entrez_id,
        "hgncSymbol": hgnc_symbol,
        "ensemblId": ensembl_id
    }
    
    result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    return result[0] if result and result[0] is not None else None


def get_gene_expression_by_gene(db_connection, cohort: str,
                                entrez_id: str = None,
                                hgnc_symbol: str = None,
                                value_type: str = "tpm",
                                sample_ids: List[str] = None) -> List[Dict]:
    """
    Get expression values for a specific gene across samples.
    
    Args:
        db_connection: ArangoDB database connection
        cohort: TCGA cohort
        entrez_id: Entrez gene ID
        hgnc_symbol: HGNC gene symbol
        value_type: "tpm", "fpkm", or "counts"
        sample_ids: Optional list of sample IDs to filter
    
    Returns:
        List of dicts with sample_id and expression value
    """
    pos = resolve_gene_position(db_connection, cohort, entrez_id=entrez_id, hgnc_symbol=hgnc_symbol)
    
    if pos is None:
        return []
    
    config = OMIC_CONFIG["gene_expression"]
    value_field = f"values_{value_type}"
    
    aql = f"""
    FOR s IN {config['collection']}
        FILTER s.data_type == @vectorDataType
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s._key IN @sampleIds
        FILTER s.{value_field} != null
        RETURN {{
            sample_id: s._key,
            value: s.{value_field}[@pos]
        }}
    """
    
    bind_vars = {
        "vectorDataType": config['vector_data_type'],
        "cohort": cohort,
        "pos": pos,
        "sampleIds": sample_ids
    }
    
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


def get_gene_expression_matrix(db_connection, cohort: str,
                               gene_ids: List[str],
                               id_type: str = "entrez_id",
                               value_type: str = "tpm",
                               sample_ids: List[str] = None) -> Dict:
    """
    Get expression matrix for multiple genes.
    
    Args:
        db_connection: ArangoDB database connection
        cohort: TCGA cohort
        gene_ids: List of gene identifiers
        id_type: "entrez_id", "hgnc_symbol", or "ensembl_id"
        value_type: "tpm", "fpkm", or "counts"
        sample_ids: Optional list of sample IDs
    
    Returns:
        Dict with 'samples', 'genes', and 'matrix' (samples x genes)
    """
    # Resolve all gene positions
    positions = {}
    for gid in gene_ids:
        kwargs = {id_type: gid} if id_type != "ensembl_id" else {"ensembl_id": gid}
        pos = resolve_gene_position(db_connection, cohort, **kwargs)
        if pos is not None:
            positions[gid] = pos
    
    if not positions:
        return {"samples": [], "genes": [], "matrix": []}
    
    pos_list = list(positions.values())
    gene_list = list(positions.keys())
    config = OMIC_CONFIG["gene_expression"]
    value_field = f"values_{value_type}"
    
    aql = f"""
    FOR s IN {config['collection']}
        FILTER s.data_type == @vectorDataType
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s._key IN @sampleIds
        FILTER s.{value_field} != null
        RETURN {{
            sample_id: s._key,
            values: (FOR p IN @positions RETURN s.{value_field}[p])
        }}
    """
    
    bind_vars = {
        "vectorDataType": config['vector_data_type'],
        "cohort": cohort,
        "positions": pos_list,
        "sampleIds": sample_ids
    }
    
    results = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    
    return {
        "samples": [r["sample_id"] for r in results],
        "genes": gene_list,
        "matrix": [r["values"] for r in results]
    }


def get_genes_vector_json(index_data, id_type = None):
    """
    Get genes vector by position from index data. for the requested id type :
    'gene_id_ensembl'
    "gene_id_base'
    "entrez_id"
    """
    # get the index_data.keys() that contains '_mappings'
    mapping_keys = [key for key in index_data.keys() if '_mappings' in key]
    mapping_data = index_data[mapping_keys[0]]
    if id_type in mapping_data[0].keys():
      genes_vector = []
      for mapping in mapping_data:
        genes_vector.append(mapping[id_type])
      return genes_vector
    else:
        options = []
        for mapping in mapping_data:
            for key in mapping.keys():
                if key not in options:
                    options.append(key)
        print(f"  Warning: Please provide id_type among available options:")
        for opt in options:
            print(f" - {opt}")
        return []

def get_gene_vector(db_connection, cohort: str,
                     id_type: str = "entrez_id",
                     omic_type: str = "gene_expression") -> List[str]:
    """
    Get ordered vector of gene identifiers from expression index for a cohort.
    
    Args:
        db_connection: ArangoDB database connection
        cohort: TCGA cohort (e.g., "TCGA-BRCA")
        id_type: Type of gene identifier to retrieve. Options include:
                 - "entrez_id"
                 - "hgnc_symbol" 
                 - "gene_id_base" (Ensembl base ID)
                 - "gene_id_ensembl" (full Ensembl ID)
        omic_type: Omic type ("gene_expression", "cnv", "methylation", "mirna", "protein")
    
    Returns:
        List of gene identifiers ordered by position in expression vectors
    """
    if omic_type not in OMIC_CONFIG:
        print(f"  Warning: Unknown omic_type '{omic_type}'. Available: {list(OMIC_CONFIG.keys())}")
        return []
    
    config = OMIC_CONFIG[omic_type]
    index_key = f"{config['index_key_prefix']}_{cohort}"
    mapping_field = config['mapping_field']
    
    aql = f"""
    FOR doc IN {config['collection']}
        FILTER doc.data_type == @indexDataType
        FILTER doc._key == @indexKey
        RETURN doc.{mapping_field}
    """
    
    bind_vars = {
        "indexKey": index_key,
        "indexDataType": config['index_data_type']
    }
    result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    
    if not result or not result[0]:
        print(f"  Warning: No index found for {omic_type} in cohort {cohort}")
        return []
    
    mapping_data = result[0]
    
    # Check if requested id_type exists in the data
    if mapping_data and id_type in mapping_data[0]:
        # Sort by position and extract the requested id type
        sorted_mappings = sorted(mapping_data, key=lambda x: x.get('position', 0))
        genes_vector = [mapping[id_type] for mapping in sorted_mappings]
        return genes_vector
    else:
        # Collect available options
        options = set()
        if mapping_data:
            for mapping in mapping_data:
                options.update(mapping.keys())
        
        print(f"  Warning: id_type '{id_type}' not found in gene mappings.")
        print(f"  Available options: {', '.join(sorted(options))}")
        return []


# =============================================================================
# CNV QUERIES
# =============================================================================

def resolve_cnv_position(db_connection, cohort: str,
                         entrez_id: str = None,
                         ensembl_id: str = None) -> Optional[int]:
    """Resolve gene to position in CNV vector."""
    config = OMIC_CONFIG["cnv"]
    index_key = f"{config['index_key_prefix']}_{cohort}"
    
    aql = f"""
    LET idx = FIRST(
        FOR doc IN {config['collection']}
            FILTER doc.data_type == @indexDataType
            FILTER doc._key == @indexKey
            RETURN doc
    )
    LET pos = FIRST(
        FOR m IN idx.{config['mapping_field']}
            FILTER (@entrezId != null AND m.entrez_id == @entrezId)
                OR (@ensemblId != null AND m.gene_id_base == @ensemblId)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
        "indexDataType": config['index_data_type'],
        "entrezId": entrez_id,
        "ensemblId": ensembl_id
    }
    
    result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    return result[0] if result and result[0] is not None else None


def get_cnv_by_gene(db_connection, cohort: str,
                    entrez_id: str = None,
                    ensembl_id: str = None,
                    sample_ids: List[str] = None) -> List[Dict]:
    """Get CNV values for a specific gene across samples."""
    pos = resolve_cnv_position(db_connection, cohort, entrez_id=entrez_id, ensembl_id=ensembl_id)
    
    if pos is None:
        return []
    
    config = OMIC_CONFIG["cnv"]
    
    aql = f"""
    FOR s IN {config['collection']}
        FILTER s.data_type == @vectorDataType
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s._key IN @sampleIds
        RETURN {{
            sample_id: s._key,
            copy_number: s.values_copy_number[@pos]
        }}
    """
    
    bind_vars = {
        "vectorDataType": config['vector_data_type'],
        "cohort": cohort, 
        "pos": pos, 
        "sampleIds": sample_ids
    }
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


# =============================================================================
# METHYLATION QUERIES
# =============================================================================

def resolve_methylation_positions_by_gene(db_connection, cohort: str,
                                          ensembl_id: str = None,
                                          gene_symbol: str = None) -> List[Dict]:
    """
    Resolve gene to methylation probe positions.
    Returns list of probes with their positions.
    """
    config = OMIC_CONFIG["methylation"]
    index_key = f"{config['index_key_prefix']}_{cohort}"
    
    aql = f"""
    LET idx = FIRST(
        FOR doc IN {config['collection']}
            FILTER doc.data_type == @indexDataType
            FILTER doc._key == @indexKey
            RETURN doc
    )
    FOR m IN idx.{config['mapping_field']}
        FILTER (@ensemblId != null AND @ensemblId IN m.gene_ids)
            OR (@geneSymbol != null AND @geneSymbol IN m.gene_symbols)
        RETURN {{
            position: m.position,
            probe_id: m.probe_id,
            chromosome: m.chromosome,
            gene_symbols: m.gene_symbols
        }}
    """
    
    bind_vars = {
        "indexKey": index_key,
        "indexDataType": config['index_data_type'],
        "ensemblId": ensembl_id,
        "geneSymbol": gene_symbol
    }
    
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


def get_methylation_by_gene(db_connection, cohort: str,
                            ensembl_id: str = None,
                            gene_symbol: str = None,
                            sample_ids: List[str] = None) -> Dict:
    """Get methylation beta values for all probes associated with a gene."""
    probes = resolve_methylation_positions_by_gene(db_connection, cohort, 
                                                    ensembl_id=ensembl_id, 
                                                    gene_symbol=gene_symbol)
    if not probes:
        return {"probes": [], "samples": []}
    
    positions = [p["position"] for p in probes]
    probe_ids = [p["probe_id"] for p in probes]
    
    config = OMIC_CONFIG["methylation"]
    
    aql = f"""
    FOR s IN {config['collection']}
        FILTER s.data_type == @vectorDataType
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s._key IN @sampleIds
        RETURN {{
            sample_id: s._key,
            beta_values: (FOR p IN @positions RETURN s.values_beta[p])
        }}
    """
    
    bind_vars = {
        "vectorDataType": config['vector_data_type'],
        "cohort": cohort, 
        "positions": positions, 
        "sampleIds": sample_ids
    }
    results = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    
    # Add probe info to results
    return {
        "probes": probe_ids,
        "samples": [{"sample_id": r["sample_id"], "values": r["beta_values"]} for r in results]
    }


# =============================================================================
# MIRNA QUERIES
# =============================================================================

def resolve_mirna_position(db_connection, cohort: str,
                           mirna_id: str = None,
                           mirbase_id: str = None) -> Optional[int]:
    """Resolve miRNA to position in expression vector."""
    config = OMIC_CONFIG["mirna"]
    index_key = f"{config['index_key_prefix']}_{cohort}"
    
    aql = f"""
    LET idx = FIRST(
        FOR doc IN {config['collection']}
            FILTER doc.data_type == @indexDataType
            FILTER doc._key == @indexKey
            RETURN doc
    )
    LET pos = FIRST(
        FOR m IN idx.{config['mapping_field']}
            FILTER (@mirnaId != null AND m.mirna_id == @mirnaId)
                OR (@mirbaseId != null AND m.mirbase_id == @mirbaseId)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
        "indexDataType": config['index_data_type'],
        "mirnaId": mirna_id,
        "mirbaseId": mirbase_id
    }
    
    result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    return result[0] if result and result[0] is not None else None


def get_mirna_expression(db_connection, cohort: str,
                         mirna_id: str = None,
                         mirbase_id: str = None,
                         sample_ids: List[str] = None) -> List[Dict]:
    """Get miRNA expression values across samples."""
    pos = resolve_mirna_position(db_connection, cohort, mirna_id=mirna_id, mirbase_id=mirbase_id)
    
    if pos is None:
        return []
    
    config = OMIC_CONFIG["mirna"]
    
    aql = f"""
    FOR s IN {config['collection']}
        FILTER s.data_type == @vectorDataType
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s._key IN @sampleIds
        RETURN {{
            sample_id: s._key,
            expression: s.values_expression[@pos]
        }}
    """
    
    bind_vars = {
        "vectorDataType": config['vector_data_type'],
        "cohort": cohort, 
        "pos": pos, 
        "sampleIds": sample_ids
    }
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


# =============================================================================
# PROTEIN QUERIES
# =============================================================================

def resolve_protein_position(db_connection, cohort: str,
                             entrez_id: str = None,
                             gene_symbol: str = None,
                             peptide_target: str = None) -> Optional[int]:
    """Resolve protein/peptide to position in RPPA vector."""
    config = OMIC_CONFIG["protein"]
    index_key = f"{config['index_key_prefix']}_{cohort}"
    
    aql = f"""
    LET idx = FIRST(
        FOR doc IN {config['collection']}
            FILTER doc.data_type == @indexDataType
            FILTER doc._key == @indexKey
            RETURN doc
    )
    LET pos = FIRST(
        FOR m IN idx.{config['mapping_field']}
            FILTER (@entrezId != null AND m.entrez_id == @entrezId)
                OR (@geneSymbol != null AND m.gene_symbol == @geneSymbol)
                OR (@peptideTarget != null AND m.peptide_target == @peptideTarget)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
        "indexDataType": config['index_data_type'],
        "entrezId": entrez_id,
        "geneSymbol": gene_symbol,
        "peptideTarget": peptide_target
    }
    
    result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    return result[0] if result and result[0] is not None else None


def get_protein_abundance(db_connection, cohort: str,
                          entrez_id: str = None,
                          gene_symbol: str = None,
                          sample_ids: List[str] = None) -> List[Dict]:
    """Get protein abundance values across samples."""
    pos = resolve_protein_position(db_connection, cohort, 
                                   entrez_id=entrez_id, 
                                   gene_symbol=gene_symbol)
    if pos is None:
        return []
    
    config = OMIC_CONFIG["protein"]
    
    aql = f"""
    FOR s IN {config['collection']}
        FILTER s.data_type == @vectorDataType
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s._key IN @sampleIds
        RETURN {{
            sample_id: s._key,
            abundance: s.values_abundance[@pos]
        }}
    """
    
    bind_vars = {
        "vectorDataType": config['vector_data_type'],
        "cohort": cohort, 
        "pos": pos, 
        "sampleIds": sample_ids
    }
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


# =============================================================================
# MULTI-OMICS INTEGRATION
# =============================================================================

def get_multiomics_by_gene(db_connection, cohort: str,
                           entrez_id: str = None,
                           hgnc_symbol: str = None,
                           sample_ids: List[str] = None) -> Dict:
    """
    Get all available omics data for a gene across samples.
    
    Returns dict with expression, cnv, methylation, and protein data.
    """
    results = {
        "gene": entrez_id or hgnc_symbol,
        "cohort": cohort,
        "expression": [],
        "cnv": [],
        "methylation": {},
        "protein": []
    }
    
    # Gene Expression
    results["expression"] = get_gene_expression_by_gene(
        db_connection, cohort, 
        entrez_id=entrez_id, hgnc_symbol=hgnc_symbol,
        sample_ids=sample_ids
    )
    
    # CNV
    results["cnv"] = get_cnv_by_gene(
        db_connection, cohort,
        entrez_id=entrez_id,
        sample_ids=sample_ids
    )
    
    # Methylation
    results["methylation"] = get_methylation_by_gene(
        db_connection, cohort,
        gene_symbol=hgnc_symbol,
        sample_ids=sample_ids
    )
    
    # Protein
    results["protein"] = get_protein_abundance(
        db_connection, cohort,
        entrez_id=entrez_id, gene_symbol=hgnc_symbol,
        sample_ids=sample_ids
    )
    
    return results


def get_sample_full_profile(db_connection, sample_id: str, cohort: str) -> Dict:
    """Get complete omics profile for a single sample."""
    
    expr_config = OMIC_CONFIG["gene_expression"]
    cnv_config = OMIC_CONFIG["cnv"]
    meth_config = OMIC_CONFIG["methylation"]
    mirna_config = OMIC_CONFIG["mirna"]
    prot_config = OMIC_CONFIG["protein"]
    
    aql = f"""
    LET expr = FIRST(FOR s IN {expr_config['collection']} 
                     FILTER s.data_type == @exprVectorType
                     FILTER s._key == @sampleId RETURN s)
    LET cnv = FIRST(FOR s IN {cnv_config['collection']} 
                    FILTER s.data_type == @cnvVectorType
                    FILTER s._key == @sampleId RETURN s)
    LET meth = FIRST(FOR s IN {meth_config['collection']} 
                     FILTER s.data_type == @methVectorType
                     FILTER s._key == @sampleId RETURN s)
    LET mirna = FIRST(FOR s IN {mirna_config['collection']} 
                      FILTER s.data_type == @mirnaVectorType
                      FILTER s._key == @sampleId RETURN s)
    LET prot = FIRST(FOR s IN {prot_config['collection']} 
                     FILTER s.data_type == @protVectorType
                     FILTER s._key == @sampleId RETURN s)
    
    RETURN {{
        sample_id: @sampleId,
        cohort: @cohort,
        expression: expr.values_tpm,
        cnv: cnv.values_copy_number,
        methylation: meth.values_beta,
        mirna: mirna.values_expression,
        protein: prot.values_abundance
    }}
    """
    
    bind_vars = {
        "sampleId": sample_id, 
        "cohort": cohort,
        "exprVectorType": expr_config['vector_data_type'],
        "cnvVectorType": cnv_config['vector_data_type'],
        "methVectorType": meth_config['vector_data_type'],
        "mirnaVectorType": mirna_config['vector_data_type'],
        "protVectorType": prot_config['vector_data_type']
    }
    result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    return result[0] if result else {}


def get_sample_list_intersection(multiomics_data):
    """
    Get the intersection of sample IDs across all omics types.
    Handles different data structures and skips empty/invalid data.
    """
    sets = []
    for omic_type, data in multiomics_data.items():
        # Skip non-omics keys like 'gene', 'cohort', etc.
        if omic_type in ['gene', 'cohort']:
            continue
            
        # Handle different data structures
        if isinstance(data, dict):
            # Methylation case: {'samples': [...], ...}
            if 'samples' in data and isinstance(data['samples'], list) and data['samples']:
                # Extract sample_id from each dict in the samples list
                sample_ids = [entry.get('sample_id') for entry in data['samples'] if isinstance(entry, dict) and 'sample_id' in entry]
                if sample_ids:
                    sets.append(set(sample_ids))
        elif isinstance(data, list):
            # Expression, CNV, Protein case: list of dicts with 'sample_id'
            if data:  # Only process non-empty lists
                sample_ids = [entry.get('sample_id') for entry in data if isinstance(entry, dict) and 'sample_id' in entry]
                if sample_ids:
                    sets.append(set(sample_ids))
    
    # Return intersection only if we have at least one set
    if sets:
        return set.intersection(*sets)
    return set()


# =============================================================================
# TRANSOMIC NETWORKS QUERIES
# la costruizione delle reti trasomiche è lo scopo finale queste funzioni servono a recuperare i dati necessari per costruire le reti trasomiche e fare analisi integrative, e rilvare pattern causalità tra i diversi livelli di dati omici e il fenotipo.
# =============================================================================