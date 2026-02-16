"""
Multi-strategy KG mapping per omic type:

Omic	Strategy 1	Fallback
gene_expression	entrez_id → KG(EntrezID)	ENSG → GENES → entrez_id → KG
CNV	entrez_id → KG(EntrezID)	ENSG → GENES → entrez_id → KG
methylation	gene_ids(ENSG) → GENES → entrez_id	gene_symbols → GENES → entrez_id
miRNA	hgnc_symbol → GENES (+ normalized no dash) → entrez_id	-
protein	entrez_id → KG(EntrezID)	gene_symbol → GENES → entrez_id
"""


#%%
"""
Build transomic omic-node layer for a sample with comprehensive KG mapping.

This version focuses ONLY on:
1. Retrieving omic data (gene_expression, CNV, methylation, miRNA, protein)
2. Mapping omic features to KG nodes with multiple fallback strategies
3. Tracking mapping statistics (mapped/missing/alternative)

Edges are NOT included - they will be added in a subsequent step after 
all omic nodes are properly mapped.

Mapping Strategies by Omic Type:
- gene_expression: entrez_id → KG(EntrezID) | ENSG → GENES → entrez_id → KG
- CNV: entrez_id → KG(EntrezID) | ENSG → GENES → entrez_id → KG
- methylation: gene_ids(ENSG) → GENES → entrez_id → KG | gene_symbols → GENES → entrez_id → KG
- miRNA: hgnc_symbol → GENES → entrez_id → KG (with dash normalization)
- protein: entrez_id → KG(EntrezID) | gene_symbol → GENES → entrez_id → KG
"""

#%% Imports
import json
import os
from typing import Dict, List, Optional, Any, Tuple
from collections import defaultdict

from arangodb_utils import setup_arangodb_connection
from query_utils import OMIC_CONFIG, _fetch_index_and_vector

# =============================================================================
# MAPPING CONFIGURATION
# =============================================================================

ID_PRIORITY = {
    "gene_expression": ["entrez_id", "hgnc_symbol", "gene_id_base", "gene_id_ensembl"],
    "cnv": ["entrez_id", "gene_id_base", "gene_id_ensembl"],
    "methylation": ["probe_id"],
    "mirna": ["mirna_id", "hgnc_symbol", "mirbase_id"],
    "protein": ["gene_symbol", "entrez_id", "peptide_target"]
}

# =============================================================================
# EDGE LAYER CONFIGURATION
# Knowledge Graph edges are organized into 3 semantic layers:
# =============================================================================

EDGE_LAYERS = {
    # Layer 1: Ontological/Semantic Layer
    # Basic RDF/OWL type relationships and class hierarchy
    "ontological": {
        "description": "Ontological type relationships and class hierarchy",
        "predicate_labels": [
            "type",
            "subclass of",
            "instance of",
            "equivalent class"
        ],
        "predicate_class_codes": [
            "RDF",      # Resource Description Framework
            "OWL",      # Web Ontology Language
            "RDFS"      # RDF Schema
        ]
    },
    
    # Layer 2: Biological Layer
    # Molecular interactions, biological processes, regulation, localization
    "biological": {
        "description": "Molecular interactions, biological processes, and regulation",
        "predicate_labels": [
            "molecularly interacts with",
            "interacts with",
            "located_in",
            "location_of",
            "participates in",
            "has participant",
            "transcribed to",
            "transcribed from",
            "causally influences",
            "causally influenced by",
            "has function",
            "function of",
            "part of",
            "has part",
            "regulates",
            "regulated by",
            "positively regulates",
            "negatively regulates",
            "enables",
            "has gene product",
            "gene product of",
            "only_in_taxon"
        ],
        "predicate_class_codes": [
            "RO",       # Relation Ontology - biological relations
            "BFO",      # Basic Formal Ontology - foundational relations
            "chebi",    # Chemical Entities of Biological Interest
            "pr",       # Protein Ontology
            "CLO",      # Cell Line Ontology
            "VO",       # Vaccine Ontology
            "so",       # Sequence Ontology
            "BSPO",     # Biological Spatial Ontology
            "core",     # Core Ontology
            "pw",       # Pathway Ontology
            "cl",       # Cell Ontology
            "GO"        # Gene Ontology (if present)
        ]
    },
    
    # Layer 3: Medical/Clinical Layer
    # Disease associations, phenotypes, clinical manifestations
    "medical": {
        "description": "Disease associations, phenotypes, and clinical relationships",
        "predicate_labels": [
            "causes or contributes to condition",
            "has phenotype",
            "phenotype of",
            "is substance that treats",
            "is treated by substance",
            "risk factor for",
            "associated with",
            "manifestation of",
            "has manifestation"
        ],
        "predicate_class_codes": [
            "MONDO",    # Mondo Disease Ontology
            "mondo",    # Mondo (lowercase variant)
            "pato",     # Phenotype And Trait Ontology
            "PATO",     # PATO (uppercase variant)
            "HP",       # Human Phenotype Ontology
            "DOID",     # Disease Ontology
            "EFO",      # Experimental Factor Ontology
            "SIO"       # Semanticscience Integrated Ontology
        ]
    }
}


def classify_edge_layer(edge: Dict) -> str:
    """
    Classify an edge into one of the 3 layers based on predicate_label or predicate_class_code.
    Returns: 'ontological', 'biological', 'medical', or 'unknown'
    """
    pred_label = edge.get("predicate_label", "").lower()
    pred_class = edge.get("predicate_class_code", "")
    
    # Check each layer
    for layer_name, layer_config in EDGE_LAYERS.items():
        # Check predicate_label (case-insensitive)
        for label in layer_config["predicate_labels"]:
            if label.lower() == pred_label:
                return layer_name
        # Check predicate_class_code
        if pred_class in layer_config["predicate_class_codes"]:
            return layer_name
    
    return "unknown"


def get_layer_statistics(edges: List[Dict]) -> Dict[str, Any]:
    """
    Compute edge distribution across the 3 semantic layers.
    """
    layer_counts = defaultdict(int)
    predicate_by_layer = defaultdict(lambda: defaultdict(int))
    
    for edge in edges:
        edge_data = edge.get("edge", edge)  # Handle both formats
        layer = classify_edge_layer(edge_data)
        layer_counts[layer] += 1
        pred_label = edge_data.get("predicate_label", "unknown")
        predicate_by_layer[layer][pred_label] += 1
    
    return {
        "layer_counts": dict(layer_counts),
        "predicate_by_layer": {k: dict(v) for k, v in predicate_by_layer.items()}
    }


# =============================================================================
# KG LOOKUP FUNCTIONS
# =============================================================================

def resolve_kg_nodes_by_entrez(db_connection, entrez_ids: List[str]) -> Dict[str, Dict]:
    """Resolve entrez IDs to KG nodes (class_code=EntrezID)."""
    if not entrez_ids:
        return {}
    aql = """
    FOR n IN nodes
        FILTER n.class_code == "EntrezID"
        FILTER n.entity_id IN @ids
        RETURN { id: n.entity_id, _key: n._key, uri: n.uri, label: n.label, class_code: n.class_code }
    """
    results = list(db_connection.aql.execute(aql, bind_vars={"ids": entrez_ids}))
    return {r["id"]: r for r in results}


def resolve_kg_nodes_by_protein(db_connection, uniprot_ids: List[str]) -> Dict[str, Dict]:
    """Resolve UniProt IDs to KG protein nodes (class_code=PR)."""
    if not uniprot_ids:
        return {}
    pr_ids = [f"PR_{uid}" if not uid.startswith("PR_") else uid for uid in uniprot_ids]
    aql = """
    FOR n IN nodes
        FILTER n.class_code == "PR"
        FILTER n.entity_id IN @ids
        RETURN { id: n.entity_id, _key: n._key, uri: n.uri, label: n.label, class_code: n.class_code }
    """
    results = list(db_connection.aql.execute(aql, bind_vars={"ids": pr_ids}))
    out = {}
    for r in results:
        uid = r["id"].replace("PR_", "") if r["id"].startswith("PR_") else r["id"]
        out[uid] = r
    return out


def lookup_genes_by_ensg(db_connection, ensg_ids: List[str]) -> Dict[str, Dict]:
    """Lookup GENES collection by ENSG IDs to get entrez_id."""
    if not ensg_ids:
        return {}
    base_ids = [eid.split('.')[0] for eid in ensg_ids]
    aql = """
    FOR g IN GENES
        FILTER g._key IN @ids
        RETURN { ensg: g._key, entrez_id: g.entrez_id, hgnc_symbol: g.hgnc_symbol, uniprot_id: g.uniprot_id }
    """
    results = list(db_connection.aql.execute(aql, bind_vars={"ids": base_ids}))
    return {r["ensg"]: r for r in results}


def lookup_genes_by_symbol(db_connection, symbols: List[str]) -> Dict[str, Dict]:
    """Lookup GENES collection by HGNC symbols to get entrez_id."""
    if not symbols:
        return {}
    aql = """
    FOR g IN GENES
        FILTER g.hgnc_symbol IN @symbols
        RETURN { symbol: g.hgnc_symbol, entrez_id: g.entrez_id, ensg: g._key, uniprot_id: g.uniprot_id, gene_type: g.gene_type }
    """
    results = list(db_connection.aql.execute(aql, bind_vars={"symbols": symbols}))
    return {r["symbol"]: r for r in results}


def lookup_genes_with_biotype(db_connection, ensg_ids: List[str]) -> Dict[str, Dict]:
    """Lookup GENES collection by ENSG IDs including gene_type (biotype)."""
    if not ensg_ids:
        return {}
    base_ids = [eid.split('.')[0] for eid in ensg_ids]
    aql = """
    FOR g IN GENES
        FILTER g._key IN @ids
        RETURN { ensg: g._key, entrez_id: g.entrez_id, hgnc_symbol: g.hgnc_symbol, gene_type: g.gene_type }
    """
    results = list(db_connection.aql.execute(aql, bind_vars={"ids": base_ids}))
    return {r["ensg"]: r for r in results}


# =============================================================================
# ANALYTICS FUNCTIONS
# =============================================================================

def analyze_gene_expression_mapping(db_connection, features: List[Dict]) -> Dict[str, Any]:
    """
    Analyze gene expression mapping failures and bioentity types.
    Returns detailed breakdown of why genes did/didn't map.
    """
    # Collect all ENSG IDs
    ensg_ids = set()
    for f in features:
        ensg = f.get("gene_id_base") or (f.get("gene_id_ensembl", "").split('.')[0] if f.get("gene_id_ensembl") else None)
        if ensg:
            ensg_ids.add(ensg)
    
    # Get gene info including biotype
    genes_info = lookup_genes_with_biotype(db_connection, list(ensg_ids))
    
    # Analyze mapped vs unmapped
    mapped_biotypes = defaultdict(int)
    unmapped_biotypes = defaultdict(int)
    unmapped_reasons = defaultdict(int)
    
    for f in features:
        ensg = f.get("gene_id_base") or (f.get("gene_id_ensembl", "").split('.')[0] if f.get("gene_id_ensembl") else None)
        gene_info = genes_info.get(ensg, {}) if ensg else {}
        biotype = gene_info.get("gene_type", "unknown")
        
        if f.get("kg_node_key"):
            mapped_biotypes[biotype] += 1
        else:
            unmapped_biotypes[biotype] += 1
            # Analyze why it didn't map
            if not ensg:
                unmapped_reasons["no_ensg_id"] += 1
            elif ensg not in genes_info:
                unmapped_reasons["ensg_not_in_GENES"] += 1
            elif not gene_info.get("entrez_id"):
                unmapped_reasons["no_entrez_in_GENES"] += 1
            else:
                unmapped_reasons["entrez_not_in_KG"] += 1
    
    return {
        "mapped_by_biotype": dict(sorted(mapped_biotypes.items(), key=lambda x: -x[1])),
        "unmapped_by_biotype": dict(sorted(unmapped_biotypes.items(), key=lambda x: -x[1])),
        "unmapped_reasons": dict(sorted(unmapped_reasons.items(), key=lambda x: -x[1]))
    }


# =============================================================================
# FEATURE EXTRACTION
# =============================================================================

def extract_features_from_omic(db_connection, cohort: str, sample_id: str, omic_type: str,
                               value_field: Optional[str] = None,
                               value_abs_threshold: Optional[float] = None,
                               top_n: Optional[int] = None) -> List[Dict]:
    """Extract features from an omic vector with raw mapping info."""
    data = _fetch_index_and_vector(db_connection, cohort, sample_id, omic_type)
    idx_doc, vec_doc = data.get("index"), data.get("vector")
    
    if not idx_doc or not vec_doc:
        return []
    
    config = OMIC_CONFIG[omic_type]
    mapping_field = config["mapping_field"]
    mappings = idx_doc.get(mapping_field, []) or []
    
    chosen_value_field = value_field or config["value_fields"][0]
    values = vec_doc.get(chosen_value_field, []) or []
    
    features = []
    for mapping in mappings:
        pos = mapping.get("position")
        if pos is None or pos >= len(values):
            continue
        val = values[pos]
        
        if value_abs_threshold is not None and val is not None:
            if abs(val) < value_abs_threshold:
                continue
        
        feat = {
            "position": pos,
            "value": val,
            "value_field": chosen_value_field,
            "omic_type": omic_type,
            **mapping
        }
        features.append(feat)
    
    if top_n is not None:
        features = features[:top_n]
    
    return features


# =============================================================================
# MULTI-STRATEGY KG MAPPING
# =============================================================================

def map_gene_expression_features(db_connection, features: List[Dict]) -> Tuple[List[Dict], Dict]:
    """Map gene expression features: entrez_id → KG | ENSG → GENES → entrez_id → KG"""
    stats = {"total": len(features), "mapped_direct": 0, "mapped_via_genes": 0, "missing": 0}
    
    entrez_ids = set()
    ensg_ids = set()
    for f in features:
        if f.get("entrez_id"):
            entrez_ids.add(str(f["entrez_id"]))
        if f.get("gene_id_base"):
            ensg_ids.add(f["gene_id_base"])
        elif f.get("gene_id_ensembl"):
            ensg_ids.add(f["gene_id_ensembl"].split('.')[0])
    
    kg_by_entrez = resolve_kg_nodes_by_entrez(db_connection, list(entrez_ids))
    genes_by_ensg = lookup_genes_by_ensg(db_connection, list(ensg_ids))
    
    fallback_entrez = set()
    for ensg, gene_info in genes_by_ensg.items():
        if gene_info.get("entrez_id") and gene_info["entrez_id"] not in kg_by_entrez:
            fallback_entrez.add(str(gene_info["entrez_id"]))
    
    kg_by_entrez_fallback = resolve_kg_nodes_by_entrez(db_connection, list(fallback_entrez))
    kg_by_entrez.update(kg_by_entrez_fallback)
    
    for f in features:
        kg_hit = None
        mapping_method = None
        
        eid = str(f.get("entrez_id", "")) if f.get("entrez_id") else None
        if eid and eid in kg_by_entrez:
            kg_hit = kg_by_entrez[eid]
            mapping_method = "direct_entrez"
            stats["mapped_direct"] += 1
        else:
            ensg = f.get("gene_id_base") or (f.get("gene_id_ensembl", "").split('.')[0] if f.get("gene_id_ensembl") else None)
            if ensg and ensg in genes_by_ensg:
                gene_info = genes_by_ensg[ensg]
                fallback_eid = str(gene_info.get("entrez_id", "")) if gene_info.get("entrez_id") else None
                if fallback_eid and fallback_eid in kg_by_entrez:
                    kg_hit = kg_by_entrez[fallback_eid]
                    mapping_method = "via_genes_ensg"
                    stats["mapped_via_genes"] += 1
        
        if kg_hit:
            f["kg_node_key"] = kg_hit["_key"]
            f["kg_uri"] = kg_hit.get("uri")
            f["kg_label"] = kg_hit.get("label")
            f["kg_class_code"] = kg_hit.get("class_code")
            f["mapping_method"] = mapping_method
        else:
            f["kg_node_key"] = None
            f["mapping_method"] = "missing"
            stats["missing"] += 1
    
    return features, stats


def map_cnv_features(db_connection, features: List[Dict]) -> Tuple[List[Dict], Dict]:
    """Map CNV features - same strategy as gene_expression."""
    return map_gene_expression_features(db_connection, features)


def map_methylation_features(db_connection, features: List[Dict]) -> Tuple[List[Dict], Dict]:
    """Map methylation: gene_ids(ENSG) → GENES → entrez_id | gene_symbols → GENES → entrez_id"""
    stats = {"total": len(features), "mapped_via_ensg": 0, "mapped_via_symbol": 0, "missing": 0}
    
    ensg_ids = set()
    symbols = set()
    for f in features:
        for gid in (f.get("gene_ids") or []):
            ensg_ids.add(gid.split('.')[0])
        for sym in (f.get("gene_symbols") or []):
            symbols.add(sym)
    
    genes_by_ensg = lookup_genes_by_ensg(db_connection, list(ensg_ids))
    genes_by_symbol = lookup_genes_by_symbol(db_connection, list(symbols))
    
    all_entrez = set()
    for g in genes_by_ensg.values():
        if g.get("entrez_id"):
            all_entrez.add(str(g["entrez_id"]))
    for g in genes_by_symbol.values():
        if g.get("entrez_id"):
            all_entrez.add(str(g["entrez_id"]))
    
    kg_by_entrez = resolve_kg_nodes_by_entrez(db_connection, list(all_entrez))
    
    for f in features:
        kg_hit = None
        mapping_method = None
        mapped_gene = None
        
        for gid in (f.get("gene_ids") or []):
            ensg = gid.split('.')[0]
            if ensg in genes_by_ensg:
                gene_info = genes_by_ensg[ensg]
                eid = str(gene_info.get("entrez_id", "")) if gene_info.get("entrez_id") else None
                if eid and eid in kg_by_entrez:
                    kg_hit = kg_by_entrez[eid]
                    mapping_method = "via_ensg"
                    mapped_gene = gene_info.get("hgnc_symbol") or ensg
                    stats["mapped_via_ensg"] += 1
                    break
        
        if not kg_hit:
            for sym in (f.get("gene_symbols") or []):
                if sym in genes_by_symbol:
                    gene_info = genes_by_symbol[sym]
                    eid = str(gene_info.get("entrez_id", "")) if gene_info.get("entrez_id") else None
                    if eid and eid in kg_by_entrez:
                        kg_hit = kg_by_entrez[eid]
                        mapping_method = "via_symbol"
                        mapped_gene = sym
                        stats["mapped_via_symbol"] += 1
                        break
        
        if kg_hit:
            f["kg_node_key"] = kg_hit["_key"]
            f["kg_uri"] = kg_hit.get("uri")
            f["kg_label"] = kg_hit.get("label")
            f["kg_class_code"] = kg_hit.get("class_code")
            f["mapping_method"] = mapping_method
            f["mapped_gene"] = mapped_gene
        else:
            f["kg_node_key"] = None
            f["mapping_method"] = "missing"
            stats["missing"] += 1
    
    return features, stats


def map_mirna_features(db_connection, features: List[Dict]) -> Tuple[List[Dict], Dict]:
    """Map miRNA: hgnc_symbol → GENES → entrez_id (with dash normalization)"""
    stats = {"total": len(features), "mapped_via_symbol": 0, "missing": 0}
    
    symbols_raw = set()
    symbols_normalized = {}
    for f in features:
        sym = f.get("hgnc_symbol")
        if sym:
            symbols_raw.add(sym)
            norm = sym.replace("-", "")
            symbols_normalized[sym] = norm
    
    all_symbols = list(symbols_raw) + list(set(symbols_normalized.values()))
    genes_by_symbol = lookup_genes_by_symbol(db_connection, all_symbols)
    
    all_entrez = set()
    for g in genes_by_symbol.values():
        if g.get("entrez_id"):
            all_entrez.add(str(g["entrez_id"]))
    
    kg_by_entrez = resolve_kg_nodes_by_entrez(db_connection, list(all_entrez))
    
    for f in features:
        kg_hit = None
        mapping_method = None
        
        sym = f.get("hgnc_symbol")
        if sym:
            if sym in genes_by_symbol:
                gene_info = genes_by_symbol[sym]
                eid = str(gene_info.get("entrez_id", "")) if gene_info.get("entrez_id") else None
                if eid and eid in kg_by_entrez:
                    kg_hit = kg_by_entrez[eid]
                    mapping_method = "via_symbol"
            
            if not kg_hit:
                norm = symbols_normalized.get(sym)
                if norm and norm in genes_by_symbol:
                    gene_info = genes_by_symbol[norm]
                    eid = str(gene_info.get("entrez_id", "")) if gene_info.get("entrez_id") else None
                    if eid and eid in kg_by_entrez:
                        kg_hit = kg_by_entrez[eid]
                        mapping_method = "via_symbol_normalized"
        
        if kg_hit:
            f["kg_node_key"] = kg_hit["_key"]
            f["kg_uri"] = kg_hit.get("uri")
            f["kg_label"] = kg_hit.get("label")
            f["kg_class_code"] = kg_hit.get("class_code")
            f["mapping_method"] = mapping_method
            stats["mapped_via_symbol"] += 1
        else:
            f["kg_node_key"] = None
            f["mapping_method"] = "missing"
            stats["missing"] += 1
    
    return features, stats


def map_protein_features(db_connection, features: List[Dict]) -> Tuple[List[Dict], Dict]:
    """Map protein: entrez_id → KG | gene_symbol → GENES → entrez_id"""
    stats = {"total": len(features), "mapped_entrez": 0, "mapped_via_symbol": 0, "missing": 0}
    
    entrez_ids = set()
    symbols = set()
    for f in features:
        if f.get("entrez_id"):
            entrez_ids.add(str(f["entrez_id"]))
        if f.get("gene_symbol"):
            symbols.add(f["gene_symbol"])
    
    kg_by_entrez = resolve_kg_nodes_by_entrez(db_connection, list(entrez_ids))
    genes_by_symbol = lookup_genes_by_symbol(db_connection, list(symbols))
    
    for g in genes_by_symbol.values():
        if g.get("entrez_id"):
            eid = str(g["entrez_id"])
            if eid not in kg_by_entrez:
                entrez_ids.add(eid)
    
    kg_by_entrez = resolve_kg_nodes_by_entrez(db_connection, list(entrez_ids))
    
    for f in features:
        kg_hit = None
        mapping_method = None
        
        eid = str(f.get("entrez_id", "")) if f.get("entrez_id") else None
        if eid and eid in kg_by_entrez:
            kg_hit = kg_by_entrez[eid]
            mapping_method = "direct_entrez"
            stats["mapped_entrez"] += 1
        else:
            sym = f.get("gene_symbol")
            if sym and sym in genes_by_symbol:
                gene_info = genes_by_symbol[sym]
                fallback_eid = str(gene_info.get("entrez_id", "")) if gene_info.get("entrez_id") else None
                if fallback_eid and fallback_eid in kg_by_entrez:
                    kg_hit = kg_by_entrez[fallback_eid]
                    mapping_method = "via_symbol"
                    stats["mapped_via_symbol"] += 1
        
        if kg_hit:
            f["kg_node_key"] = kg_hit["_key"]
            f["kg_uri"] = kg_hit.get("uri")
            f["kg_label"] = kg_hit.get("label")
            f["kg_class_code"] = kg_hit.get("class_code")
            f["mapping_method"] = mapping_method
        else:
            f["kg_node_key"] = None
            f["mapping_method"] = "missing"
            stats["missing"] += 1
    
    return features, stats


MAPPING_FUNCTIONS = {
    "gene_expression": map_gene_expression_features,
    "cnv": map_cnv_features,
    "methylation": map_methylation_features,
    "mirna": map_mirna_features,
    "protein": map_protein_features
}


# =============================================================================
# MAIN BUILD FUNCTION
# =============================================================================

def build_transomic_omic_layer(
    db_connection,
    cohort: str,
    sample_id: str,
    omic_types: Optional[List[str]] = None,
    value_field_override: Optional[Dict[str, str]] = None,
    value_abs_threshold: Optional[float] = None,
    top_n: Optional[int] = None,
    include_analytics: bool = True
) -> Dict[str, Any]:
    """Build the omic layer of a transomic network for a sample."""
    omic_types = omic_types or list(OMIC_CONFIG.keys())
    value_field_override = value_field_override or {}
    
    all_nodes = []
    all_stats = {}
    analytics = {}
    
    for omic in omic_types:
        print(f"Processing {omic}...")
        
        features = extract_features_from_omic(
            db_connection,
            cohort=cohort,
            sample_id=sample_id,
            omic_type=omic,
            value_field=value_field_override.get(omic),
            value_abs_threshold=value_abs_threshold,
            top_n=top_n
        )
        
        if not features:
            print(f"  No features found for {omic}")
            all_stats[omic] = {"total": 0, "missing": 0}
            continue
        
        map_fn = MAPPING_FUNCTIONS.get(omic)
        if map_fn:
            features, stats = map_fn(db_connection, features)
            all_stats[omic] = stats
            print(f"  {omic}: {stats}")
            
            # Run analytics for gene_expression
            if include_analytics and omic == "gene_expression":
                analytics["gene_expression"] = analyze_gene_expression_mapping(db_connection, features)
        else:
            all_stats[omic] = {"total": len(features), "missing": len(features)}
        
        all_nodes.extend(features)
    
    total_stats = {
        "total_features": sum(s.get("total", 0) for s in all_stats.values()),
        "total_mapped": sum(s.get("total", 0) - s.get("missing", 0) for s in all_stats.values()),
        "total_missing": sum(s.get("missing", 0) for s in all_stats.values())
    }
    total_stats["mapping_rate"] = (
        round(100 * total_stats["total_mapped"] / total_stats["total_features"], 2)
        if total_stats["total_features"] > 0 else 0
    )
    
    return {
        "metadata": {
            "cohort": cohort,
            "sample_id": sample_id,
            "omic_types": omic_types,
            "value_abs_threshold": value_abs_threshold,
            "top_n": top_n
        },
        "nodes": all_nodes,
        "edges": [],
        "mapping_stats": {
            "per_omic": all_stats,
            "total": total_stats
        },
        "analytics": analytics
    }


# =============================================================================
# CLI / NOTEBOOK INTERFACE
# =============================================================================

def build_and_save(
    db_name: str,
    cohort: str,
    sample_id: str,
    omic_types: Optional[List[str]] = None,
    value_abs_threshold: Optional[float] = None,
    top_n: Optional[int] = None,
    output_path: Optional[str] = None
) -> Tuple[Dict, str]:
    """Build transomic omic layer and save to JSON."""
    db = setup_arangodb_connection(db_name)
    if db is None:
        raise RuntimeError("Failed to connect to ArangoDB")
    
    graph = build_transomic_omic_layer(
        db_connection=db,
        cohort=cohort,
        sample_id=sample_id,
        omic_types=omic_types,
        value_abs_threshold=value_abs_threshold,
        top_n=top_n
    )
    
    if output_path is None:
        safe_sample = sample_id.replace("/", "-")
        output_path = os.path.join("..", "transomic-networks", f"transomic_omics_{cohort}_{safe_sample}.json")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2)
    
    return graph, output_path


#%% Quick-start cell
if __name__ == "__main__":
    DB_NAME = "PKT_test10000"
    COHORT = "TCGA-BRCA"
    SAMPLE_ID = 'TCGA-A7-A0CD-01A' # with comprehensive data and good mapping rates
    # SAMPLE_ID = "TCGA-BH-A1F2-01A" # no methilation data
    
    graph, path = build_and_save(
        db_name=DB_NAME,
        cohort=COHORT,
        sample_id=SAMPLE_ID,
        omic_types=None,
        value_abs_threshold=None,
        top_n=None
    )
    
    print(f"\n{'='*60}")
    print("MAPPING SUMMARY")
    print('='*60)
    stats = graph["mapping_stats"]
    for omic, s in stats["per_omic"].items():
        mapped = s.get("total", 0) - s.get("missing", 0)
        total = s.get("total", 0)
        pct = round(100 * mapped / total, 1) if total > 0 else 0
        print(f"{omic:20s}: {mapped:6d}/{total:6d} mapped ({pct}%)")
    
    print("-"*60)
    t = stats["total"]
    print(f"{'TOTAL':20s}: {t['total_mapped']:6d}/{t['total_features']:6d} mapped ({t['mapping_rate']}%)")
    
    # Print biotype analytics for gene_expression
    if "analytics" in graph and "gene_expression" in graph["analytics"]:
        analytics = graph["analytics"]["gene_expression"]
        
        print(f"\n{'='*60}")
        print("GENE EXPRESSION BIOTYPE ANALYSIS")
        print('='*60)
        
        print("\n--- MAPPED genes by biotype ---")
        for biotype, count in list(analytics["mapped_by_biotype"].items())[:15]:
            print(f"  {biotype:30s}: {count:6d}")
        
        print("\n--- UNMAPPED genes by biotype ---")
        for biotype, count in list(analytics["unmapped_by_biotype"].items())[:15]:
            print(f"  {biotype:30s}: {count:6d}")
        
        print("\n--- UNMAPPED reasons ---")
        for reason, count in analytics["unmapped_reasons"].items():
            print(f"  {reason:30s}: {count:6d}")
    
    print(f"\nSaved to: {path}")

#%%

"""
============================================================
MAPPING SUMMARY
============================================================
gene_expression     :  25695/ 60660 mapped (42.4%)
cnv                 :  25659/ 60623 mapped (42.3%)
methylation         :  25725/ 27578 mapped (93.3%)
mirna               :   1831/  1881 mapped (97.3%)
protein             :    478/   487 mapped (98.2%)
------------------------------------------------------------
TOTAL               :  79388/151229 mapped (52.5%)

============================================================
GENE EXPRESSION BIOTYPE ANALYSIS
============================================================

--- MAPPED genes by biotype ---
  protein_coding                :  19290
  lncRNA                        :   3027
  miRNA                         :   1836
  transcribed_unprocessed_pseudogene:    405
  snoRNA                        :    377
  IG_V_gene                     :    128
  transcribed_processed_pseudogene:    107
  TR_V_gene                     :    103
  transcribed_unitary_pseudogene:     99
  processed_pseudogene          :     80
  snRNA                         :     44
  unprocessed_pseudogene        :     41
  Mt_tRNA                       :     21
  rRNA                          :     20
  scaRNA                        :     20

--- UNMAPPED genes by biotype ---
  unknown                       :  18099
  processed_pseudogene          :   6345
  lncRNA                        :   2739
  misc_RNA                      :   1952
  snRNA                         :   1854
  unprocessed_pseudogene        :   1153
  transcribed_unprocessed_pseudogene:    747
  transcribed_processed_pseudogene:    652
  rRNA_pseudogene               :    497
  snoRNA                        :    229
  IG_V_pseudogene               :    182
  protein_coding                :    141
  transcribed_unitary_pseudogene:     67
  unitary_pseudogene            :     67
  TR_J_gene                     :     61

--- UNMAPPED reasons ---
  ensg_not_in_GENES             :  18099
  entrez_not_in_KG              :  15909
  no_entrez_in_GENES            :    957


Saved to: ..\transomic-networks\transomic_omics_TCGA-BRCA_TCGA-A7-A0CD-01A.json
"""


#%%
#=============================================================================
# SECONBD LEVEL - ADD EDGES BASED ON KG RELATIONSHIPS
#=============================================================================

"""
Total edges: 11082103

Property                  Non-Null     Null         Coverage %
------------------------------------------------------------
target_uri                11082103     0            100.0     
source_uri                11082103     0            100.0     
predicate_uri             11082103     0            100.0     
predicate_class_code      11082103     0            100.0     
predicate_bioentity_type  11082103     0            100.0     
predicate_source          11082103     0            100.0     
predicate_label           11082094     9            100.0  

--- predicate_label value distribution (top 20) ---
  type: 4521542
  molecularly interacts with: 1384824
  located_in: 689135
  location_of: 688940
  phenotype of: 428374
  has phenotype: 428374
  has participant: 381760
  participates in: 381721
  is substance that treats: 279052
  is treated by substance: 279052
  transcribed to: 182692
  transcribed from: 182692
  interacts with: 176638
  causally influenced by: 145129
  causally influences: 145129
  only_in_taxon: 86315
  causes or contributes to condition: 83530
  has function: 72230
  function of: 72227
  part of: 67365

--- predicate_class_code value distribution (top 20) ---
  RO: 6362629
  RDF: 4521542
  BFO: 119618
  chebi: 41661
  pr: 19845
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
  envo: 8
  exo.obo: 5
  OGG: 4

"""

"""
EDGE LAYER CLASSIFICATION (see EDGE_LAYERS config above):

1. ONTOLOGICAL EDGES - Type relationships and class hierarchy
   predicate_label: type, subclass of, instance of, equivalent class
   predicate_class_code: RDF, OWL, RDFS

2. BIOLOGICAL EDGES - Molecular interactions, processes, regulation
   predicate_label: molecularly interacts with, interacts with, located_in, location_of, 
                    participates in, has participant, transcribed to/from, 
                    causally influences/influenced by, has function/function of, 
                    part of/has part, regulates, only_in_taxon
   predicate_class_code: RO, BFO, chebi, pr, CLO, VO, so, BSPO, core, pw, cl, GO

3. MEDICAL EDGES - Disease associations and phenotypes
   predicate_label: causes or contributes to condition, has phenotype, phenotype of,
                    is substance that treats, is treated by substance
   predicate_class_code: MONDO, mondo, pato, PATO, HP, DOID, EFO, SIO
"""

"""
EDGES RELEVANT FOR TRANSOMIC NETWORKS:

In una rete transomica ci interessano relazioni che collegano entità omiche 
(geni, proteine, metaboliti) con significato biologico/funzionale.

✓ INCLUDE (high priority):
  - molecularly interacts with  → interazioni proteina-proteina, gene-proteina
  - interacts with              → interazioni generiche
  - transcribed to/from         → relazioni gene → trascritto
  - regulates / regulated by    → regolazione genica
  - positively/negatively regulates → direzione della regolazione
  - causally influences/influenced by → relazioni causali
  - has gene product / gene product of → gene → proteina
  - participates in / has participant → coinvolgimento in pathway

✓ INCLUDE (medium priority - context):
  - has function / function of  → annotazioni funzionali (GO terms)
  - part of / has part          → relazioni strutturali
  - causes or contributes to condition → link gene-malattia
  - has phenotype / phenotype of → link gene-fenotipo

✗ EXCLUDE (not useful for transomic):
  - type / subclass of          → solo classificazione ontologica
  - instance of / equivalent class → metadati RDF
  - located_in / location_of    → localizzazione subcellulare (meno rilevante)
  - only_in_taxon               → vincolo tassonomico
  - is substance that treats    → farmacologia (opzionale)

"""

PREDICATE_LABELS_FOR_TRANSOMIC_NETWORK={
  "CORE": [
      "molecularly interacts with", # CHEBI:PR:PR edges
      "interacts with",             # CHEBI:PR:PR edges
      "transcribed to",     # Gene:ENST edges 
      "transcribed from",   # Gene:ENST edges   
      "regulates", 
      "regulated by", 
      "positively regulates",  # GO-BP related edges
      "negatively regulates",  # GO-BP related edges
      "causally influences",    # SNP related edges
      "causally influenced by", # SNP related edges
      "has gene product", 
      "gene product of", 
      "participates in",   # Reactome-GO-BP:PR edges 
      "has participant"    # Reactome-GO-BP:PR edges
      ],
  
  "EXTENDED": [
      "has function", 
      "function of", 
      "part of", 
      "has part", 
      "causes or contributes to condition", 
      "has phenotype", 
      "phenotype of"
      ]
             }



#%%
# TEST:show example edges  for a "class_code" of RDF
from query_utils import search_edges
db = setup_arangodb_connection(DB_NAME)
edges = search_edges(db, predicate_class_code = None, 
                    predicate_label_contains = "interacts with",
                     limit=50)
for e in edges:
    print(f"{e['source_uri']} --[{e['predicate_label']} ({e['predicate_class_code']})]--> {e['target_uri']}")
    # print(e['predicate_label'])


#%%







#%%



# %%
