"""
Utility functions for querying omics data from ArangoDB.
Supports: Gene Expression, CNV, Methylation, miRNA, Protein data.
"""

from typing import Optional, List, Dict, Any, Union

# GENE DATA QUERIES
def get_gene_info(db_connection, entrez_id: str = None, 
                  hgnc_symbol: str = None,
                    ensembl_id: str = None) -> Optional[Dict]:
    """
    Retrieve gene information by various identifiers.
    
    Args:
        db_connection: ArangoDB database connection
        entrez_id: Entrez gene ID
        hgnc_symbol: HGNC gene symbol
        ensembl_id: Ensembl gene ID base
    
    Returns:
        Gene information dict or None if not found
    """
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
    index_key = f"expr_index_{cohort}"
    
    aql = """
    LET idx = DOCUMENT("GENE_EXPRESSION_INDEX", @indexKey)
    LET pos = FIRST(
        FOR m IN idx.gene_mappings
            FILTER (@entrezId != null AND m.entrez_id == @entrezId)
                OR (@hgncSymbol != null AND m.hgnc_symbol == @hgncSymbol)
                OR (@ensemblId != null AND m.gene_id_base == @ensemblId)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
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
    
    value_field = f"values_{value_type}"
    
    aql = f"""
    FOR s IN GENE_EXPRESSION_SAMPLES
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s.sample_id IN @sampleIds
        FILTER s.{value_field} != null
        RETURN {{
            sample_id: s.sample_id,
            value: s.{value_field}[@pos]
        }}
    """
    
    bind_vars = {
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
    value_field = f"values_{value_type}"
    
    aql = f"""
    FOR s IN GENE_EXPRESSION_SAMPLES
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s.sample_id IN @sampleIds
        FILTER s.{value_field} != null
        RETURN {{
            sample_id: s.sample_id,
            values: (FOR p IN @positions RETURN s.{value_field}[p])
        }}
    """
    
    bind_vars = {
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
                     INDEX = "GENE_EXPRESSION_INDEX") -> List[str]:
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
    
    Returns:
        List of gene identifiers ordered by position in expression vectors
    """
    # Build index key based on the INDEX collection type
    # Map INDEX collection names to their key prefixes and mapping field names
    index_config_map = {
        "GENE_EXPRESSION_INDEX": {"prefix": "expr", "mapping_field": "gene_mappings"},
        "CNV_INDEX": {"prefix": "cnv", "mapping_field": "gene_mappings"},
        "METHYLATION_INDEX": {"prefix": "methylation", "mapping_field": "probe_mappings"},
        "MIRNA_INDEX": {"prefix": "mirna", "mapping_field": "mirna_mappings"},
        "PROTEIN_INDEX": {"prefix": "protein", "mapping_field": "protein_mappings"}
    }
    
    config = index_config_map.get(INDEX, {"prefix": "expr", "mapping_field": "gene_mappings"})
    index_key = f"{config['prefix']}_index_{cohort}"
    mapping_field = config['mapping_field']
    
    aql = f"""
    LET idx = DOCUMENT("{INDEX}", @indexKey)
    RETURN idx.{mapping_field}
    """
    
    bind_vars = {"indexKey": index_key}
    result = list(db_connection.aql.execute(aql, bind_vars=bind_vars))
    
    if not result or not result[0]:
        print(f"  Warning: No index found for {INDEX} in cohort {cohort}")
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
    index_key = f"cnv_index_{cohort}"
    
    aql = """
    LET idx = DOCUMENT("CNV_INDEX", @indexKey)
    LET pos = FIRST(
        FOR m IN idx.gene_mappings
            FILTER (@entrezId != null AND m.entrez_id == @entrezId)
                OR (@ensemblId != null AND m.gene_id_base == @ensemblId)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
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
    
    aql = """
    FOR s IN CNV_SAMPLES
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s.sample_id IN @sampleIds
        RETURN {
            sample_id: s.sample_id,
            copy_number: s.values_copy_number[@pos]
        }
    """
    
    bind_vars = {"cohort": cohort, "pos": pos, "sampleIds": sample_ids}
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
    index_key = f"methylation_index_{cohort}"
    
    aql = """
    LET idx = DOCUMENT("METHYLATION_INDEX", @indexKey)
    FOR m IN idx.probe_mappings
        FILTER (@ensemblId != null AND @ensemblId IN m.gene_ids)
            OR (@geneSymbol != null AND @geneSymbol IN m.gene_symbols)
        RETURN {
            position: m.position,
            probe_id: m.probe_id,
            chromosome: m.chromosome,
            gene_symbols: m.gene_symbols
        }
    """
    
    bind_vars = {
        "indexKey": index_key,
        "ensemblId": ensembl_id,
        "geneSymbol": gene_symbol
    }
    
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


def get_methylation_by_gene(db_connection, cohort: str,
                            ensembl_id: str = None,
                            gene_symbol: str = None,
                            sample_ids: List[str] = None) -> List[Dict]:
    """Get methylation beta values for all probes associated with a gene."""
    probes = resolve_methylation_positions_by_gene(db_connection, cohort, 
                                                    ensembl_id=ensembl_id, 
                                                    gene_symbol=gene_symbol)
    if not probes:
        return []
    
    positions = [p["position"] for p in probes]
    probe_ids = [p["probe_id"] for p in probes]
    
    aql = """
    FOR s IN METHYLATION_SAMPLES
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s.sample_id IN @sampleIds
        RETURN {
            sample_id: s.sample_id,
            beta_values: (FOR p IN @positions RETURN s.values_beta[p])
        }
    """
    
    bind_vars = {"cohort": cohort, "positions": positions, "sampleIds": sample_ids}
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
    index_key = f"mirna_index_{cohort}"
    
    aql = """
    LET idx = DOCUMENT("MIRNA_INDEX", @indexKey)
    LET pos = FIRST(
        FOR m IN idx.mirna_mappings
            FILTER (@mirnaId != null AND m.mirna_id == @mirnaId)
                OR (@mirbaseId != null AND m.mirbase_id == @mirbaseId)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
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
    
    aql = """
    FOR s IN MIRNA_SAMPLES
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s.sample_id IN @sampleIds
        RETURN {
            sample_id: s.sample_id,
            expression: s.values_expression[@pos]
        }
    """
    
    bind_vars = {"cohort": cohort, "pos": pos, "sampleIds": sample_ids}
    return list(db_connection.aql.execute(aql, bind_vars=bind_vars))


# =============================================================================
# PROTEIN QUERIES
# =============================================================================

def resolve_protein_position(db_connection, cohort: str,
                             entrez_id: str = None,
                             gene_symbol: str = None,
                             peptide_target: str = None) -> Optional[int]:
    """Resolve protein/peptide to position in RPPA vector."""
    index_key = f"protein_index_{cohort}"
    
    aql = """
    LET idx = DOCUMENT("PROTEIN_INDEX", @indexKey)
    LET pos = FIRST(
        FOR m IN idx.protein_mappings
            FILTER (@entrezId != null AND m.entrez_id == @entrezId)
                OR (@geneSymbol != null AND m.gene_symbol == @geneSymbol)
                OR (@peptideTarget != null AND m.peptide_target == @peptideTarget)
            RETURN m.position
    )
    RETURN pos
    """
    
    bind_vars = {
        "indexKey": index_key,
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
    
    aql = """
    FOR s IN PROTEIN_SAMPLES
        FILTER s.cohort == @cohort
        FILTER @sampleIds == null OR s.sample_id IN @sampleIds
        RETURN {
            sample_id: s.sample_id,
            abundance: s.values_abundance[@pos]
        }
    """
    
    bind_vars = {"cohort": cohort, "pos": pos, "sampleIds": sample_ids}
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
    
    aql = """
    LET expr = FIRST(FOR s IN GENE_EXPRESSION_SAMPLES 
                     FILTER s.sample_id == @sampleId RETURN s)
    LET cnv = FIRST(FOR s IN CNV_SAMPLES 
                    FILTER s.sample_id == @sampleId RETURN s)
    LET meth = FIRST(FOR s IN METHYLATION_SAMPLES 
                     FILTER s.sample_id == @sampleId RETURN s)
    LET mirna = FIRST(FOR s IN MIRNA_SAMPLES 
                      FILTER s.sample_id == @sampleId RETURN s)
    LET prot = FIRST(FOR s IN PROTEIN_SAMPLES 
                     FILTER s.sample_id == @sampleId RETURN s)
    
    RETURN {
        sample_id: @sampleId,
        cohort: @cohort,
        expression: expr.values_tpm,
        cnv: cnv.values_copy_number,
        methylation: meth.values_beta,
        mirna: mirna.values_expression,
        protein: prot.values_abundance
    }
    """
    
    bind_vars = {"sampleId": sample_id, "cohort": cohort}
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
