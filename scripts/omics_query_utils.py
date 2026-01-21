"""
Utility functions for querying omics data from ArangoDB.
Supports: Gene Expression, CNV, Methylation, miRNA, Protein data.
"""

from typing import Optional, List, Dict, Any, Union

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