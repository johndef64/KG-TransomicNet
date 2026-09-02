"""
load_graph_to_arangodb.py
=========================
Load the property-graph JSON files produced by build_property_graph.py
(nodes.json + edges.json) into an ArangoDB database.

The script always loads the full graph (both nodes and edges) -- there is no
selective loading because the semantic backbone is a single coherent unit.

Usage
-----
    # Default: read from data/pkt/builds/v3.0.2/property_graph/, db = PKT_main
    python scripts/load_graph_to_arangodb.py

    # Custom database name and JSON directory
    python scripts/load_graph_to_arangodb.py --db PKT_main \\
        --input-dir data/pkt/builds/v3.0.2/property_graph/sample_20260514_120000

    # Load from a graph-format JSON file (single file with 'nodes' and 'edges' keys)
    python scripts/load_graph_to_arangodb.py --graph-file path/to/graph.json

    # Skip the automatic fix of edges referencing the nodes/Summary placeholder
    python scripts/load_graph_to_arangodb.py --no-fix-summary-edges

    # Custom ArangoDB endpoint / credentials
    python scripts/load_graph_to_arangodb.py \\
        --host http://my-server:8529 --user myuser --password mypass --db MyDB

Default connection settings come from scripts/arangodb_utils.py
(arangodb_hosts / arangodb_user / arangodb_password). The CLI flags override
those module-level defaults at runtime.
"""

import argparse
import json
import os
import re
import sys
import traceback
from collections import defaultdict
from pathlib import Path

from arango import ArangoClient

import plotly  # noqa: F401  (imported for downstream plotting utilities)

# --- Default configuration (override via CLI) -----------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_INPUT_DIR = PROJECT_ROOT / "data" / "pkt" / "builds" / "v3.0.2" / "property_graph"
DEFAULT_DB_NAME = "PKT_main"

# Connection defaults (kept for backwards compatibility with importing scripts)
db_name = DEFAULT_DB_NAME
arangodb_hosts = "http://localhost:8529"
arangodb_user = "root"
arangodb_password = "avocadodb"
data_dir = str(DEFAULT_INPUT_DIR) + "/"
NODES_FILE = str(DEFAULT_INPUT_DIR / "nodes.json")
EDGES_FILE = str(DEFAULT_INPUT_DIR / "edges.json")
BATCH_SIZE = 1000  # Batch size for ArangoDB insert_many operations

""" The structure of the JSON data to be loaded are this

edge.json
[
  {
    "edge_id": "edge_1518276",
    "source_uri": "http://purl.obolibrary.org/obo/PR_Q9H609",
    "target_uri": "http://www.ncbi.nlm.nih.gov/gene/79177",
    "predicate_uri": "http://purl.obolibrary.org/obo/RO_0002204",
    "predicate_label": "gene product of",
    "predicate_class_code": "RO",
    "predicate_bioentity_type": "phenotype",
    "predicate_source": "Relation Ontology"
  },
  {
    "edge_id": "edge_7453030",
    "source_uri": "http://purl.obolibrary.org/obo/CLO_0023083",
    "target_uri": "http://purl.obolibrary.org/obo/CLO_0000652",
    "predicate_uri": "http://www.w3.org/1999/02/22-rdf-syntax-ns#type",
    "predicate_label": "type",
    "predicate_class_code": "RDF",
    "predicate_bioentity_type": "go",
    "predicate_source": "Resource Description Framework"
  },
  {

node.json
[
  {
    "_key": "PR_Q9H609",
    "uri": "http://purl.obolibrary.org/obo/PR_Q9H609",
    "namespace": "purl.obolibrary.org",
    "entity_id": "PR_Q9H609",
    "class_code": "PR",
    "label": "zinc finger protein 576 (human)",
    "bioentity_type": "protein",
    "description": "A zinc finger protein 576 that is encoded in the genome of human.",
    "synonym": "hZNF576|ZNF576",
    "source": "Protein Ontology",
    "source_type": "Ontology",
    "integer_id": 39530
  },
  {
    "_key": "79177",
    "uri": "http://www.ncbi.nlm.nih.gov/gene/79177",
    "namespace": "www.ncbi.nlm.nih.gov",
    "entity_id": "79177",
    "class_code": "EntrezID",
    "label": "ZNF576 (human)",
    "bioentity_type": "gene",
    "description": "A protein coding gene ZNF576 in human.",
    "synonym": "",
    "source": "NCBI Entrez Gene",
    "source_type": "Database",
    "integer_id": 216413
  },
"""


# --- ArangoDB core functions ---

# from arangodb_utils import *
from arangodb_utils import setup_arangodb_connection, get_collections_data, get_node_centric_graph, plot_subgraph

# Create collections with appropriate types and indexes

def create_arangodb_collections(db_connection,
                                nodes_collection='nodes',
                                edges_collection='edges'):
    """Create nodes (vertex) and edges (edge) collections"""
    if db_connection is None: return

    try:
        # Drop existing collections to allow a clean reload
        for collection_name in [nodes_collection, edges_collection]:
            if db_connection.has_collection(collection_name):
                db_connection.delete_collection(collection_name)
                print(f"Deleted existing collection: {collection_name}")

        # Vertex collection for nodes
        db_connection.create_collection(nodes_collection)
        print("[OK] Created nodes (vertex) collection")

        # Edge collection for relations
        db_connection.create_collection(edges_collection, edge=True)
        print("[OK] Created edges (edge) collection")

        # Hash indexes for common lookups
        try:
            db_connection.collection(nodes_collection).add_index({'fields': ['label'], 'type': 'hash'})
            db_connection.collection(nodes_collection).add_index({'fields': ['entity_id'], 'type': 'hash'})
            db_connection.collection(edges_collection).add_index({'fields': ['relationship'], 'type': 'hash'})
            print("✓ Created indexes for better query performance")
        except Exception as e:
            print(f"Warning: Could not create indexes: {e}")
            
    except Exception as e:
        print(f"✗ Error creating ArangoDB collections: {e}")
        sys.exit(1)


# --- Data import function ---

def import_data_to_arangodb(db_connection,
                            nodes_collection='nodes',
                            edges_collection='edges',
                            nodes_file=NODES_FILE,
                            edges_file=EDGES_FILE,
                            graph_file=None):
    """Load JSON files and import data into ArangoDB"""
    if db_connection is None:
        print("✗ No database connection available")
        return
    
    try:
        if not graph_file:
            # Load nodes
            print(f"\n📂 Loading nodes from {nodes_file}...")
            with open(nodes_file, 'r', encoding='utf-8') as f:
                nodes_data = json.load(f)

            # Load edges
            print(f"📂 Loading edges from {edges_file}...")
            with open(edges_file, 'r', encoding='utf-8') as f:
                edges_data = json.load(f)
            skipt_transform = False
        else:
            # Load from a graph-format JSON file: nodes/edges are top-level keys
            print(f"\n[load] Loading graph from {graph_file}...")
            with open(graph_file, 'r', encoding='utf-8') as f:
                graph_data = json.load(f)
                nodes_data = graph_data.get('nodes', [])
                edges_data = graph_data.get('edges', [])
                # Some graph dumps wrap each edge as {edge: {...}, neighbor: {...}} -- unwrap it
                if edges_data and isinstance(edges_data, list) and isinstance(edges_data[0], dict) and 'edge' in edges_data[0]:
                    print("[warn] Detected nested {edge, neighbor} format; unwrapping edge properties.")
                    edges_data = [e.get('edge', {}) for e in edges_data]
                skipt_transform = True
        
        print(f"✓ Loaded {len(nodes_data)} nodes and {len(edges_data)} edges from JSON files")
        
        # Get collections
        nodes_collection = db_connection.collection(nodes_collection)
        edges_collection = db_connection.collection(edges_collection)

          
        # Import nodes in batches
        print(f"\n🔄 Importing nodes in batches of {BATCH_SIZE}...")
        nodes_imported = 0
        nodes_failed = 0
        
        for i in range(0, len(nodes_data), BATCH_SIZE):
            batch = nodes_data[i:i + BATCH_SIZE]
            try:
                result = nodes_collection.insert_many(batch, overwrite=False, silent=True)
                nodes_imported += len(batch)
                if (i // BATCH_SIZE + 1) % 10 == 0 or i + BATCH_SIZE >= len(nodes_data):
                    print(f"  Processed {min(i + BATCH_SIZE, len(nodes_data))}/{len(nodes_data)} nodes...")
            except Exception as e:
                nodes_failed += len(batch)
                print(f"  ⚠ Warning: Batch {i//BATCH_SIZE + 1} failed: {e}")
        
        print(f"✓ Imported {nodes_imported} nodes (failed: {nodes_failed})")
            
        # Import edges in batches
        print(f"\n🔄 Importing edges in batches of {BATCH_SIZE}...")
        edges_imported = 0
        edges_failed = 0
        
        for i in range(0, len(edges_data), BATCH_SIZE):
            batch = edges_data[i:i + BATCH_SIZE]
            
            # Transform edges to ArangoDB format while preserving ALL properties
            if not skipt_transform:
                transformed_batch = []
                for edge in batch: 
                    # Extract entity_id from URIs (last component after /)
                    source_key = edge['source_uri'].split('/')[-1].split('?t=')[-1]  # Handle query params
                    target_key = edge['target_uri'].split('/')[-1].split('?t=')[-1]
                    
                    # Create edge document with _from/_to and ALL original properties
                    edge_doc = {
                        '_key': edge['edge_id'],
                        '_from': f"nodes/{source_key}",
                        '_to': f"nodes/{target_key}",
                        # Preserve ALL edge properties from the original JSON
                        **{k: v for k, v in edge.items() if k != 'edge_id'}
                    }
                    transformed_batch.append(edge_doc)
            else:
                print("⚠ Skipping edge transformation as data appears to already be in ArangoDB format")
                transformed_batch = batch
                
            try:
                result = edges_collection.insert_many(transformed_batch, overwrite=False, silent=True)
                edges_imported += len(transformed_batch)
                if (i // BATCH_SIZE + 1) % 10 == 0 or i + BATCH_SIZE >= len(edges_data):
                    print(f"  Processed {min(i + BATCH_SIZE, len(edges_data))}/{len(edges_data)} edges...")
            except Exception as e:
                edges_failed += len(transformed_batch)
                print(f"  ⚠ Warning: Batch {i//BATCH_SIZE + 1} failed: {e}")
        
        print(f"✓ Imported {edges_imported} edges (failed: {edges_failed})")
        print(f"\n✅ Import completed successfully!")
        
    except FileNotFoundError as e:
        print(f"✗ Error: Could not find data files: {e}")
        print(f"  Please ensure {NODES_FILE} and {EDGES_FILE} exist")
    except json.JSONDecodeError as e:
        print(f"✗ Error: Invalid JSON format: {e}")
    except Exception as e:
        print(f"✗ Unexpected error during import: {e}")
        traceback.print_exc()


def _extract_enst_id(uri):
    """Return ENST identifier found after the Summary query parameter in a URI."""
    if not uri:
        return None

    marker = "Summary?t="
    candidate = uri.split(marker, 1)[-1] if marker in uri else uri
    candidate = candidate.split('&', 1)[0]
    candidate = candidate.split('#', 1)[0]
    candidate = candidate.split('/', 1)[0]

    match = re.search(r"ENST\d+", candidate)
    return match.group(0) if match else None

def fix_summary_edges(db_connection):
    """Replace placeholder nodes/Summary references in edges with the actual ENST ids."""
    if db_connection is None:
        print("✗ No database connection available for ENST edge fix")
        return

    print("\n🔧 Fixing edges referencing nodes/Summary ...")
    query = """
    FOR edge IN edges
        FILTER edge._from == "nodes/Summary" OR edge._to == "nodes/Summary"
        RETURN {
            _key: edge._key,
            _from: edge._from,
            _to: edge._to,
            source_uri: edge.source_uri,
            target_uri: edge.target_uri
        }
    """

    edges_to_fix = list(db_connection.aql.execute(query))
    print(f"Found {len(edges_to_fix)} edges with placeholder references to nodes/Summary")
    if not edges_to_fix:
        print("✓ No edges with placeholder references detected")
        return

    edges_collection = db_connection.collection('edges')
    updated = 0
    skipped = 0
    from tqdm import tqdm
    for edge in tqdm(edges_to_fix, desc="Fixing edges"):
        update_doc = {'_key': edge['_key']}
        missing_identifier = False

        if edge['_from'] == 'nodes/Summary':
            enst_id = _extract_enst_id(edge.get('source_uri')) or _extract_enst_id(edge.get('target_uri'))
            if enst_id:
                update_doc['_from'] = f"nodes/{enst_id}"
            else:
                missing_identifier = True

        if edge['_to'] == 'nodes/Summary':
            enst_id = _extract_enst_id(edge.get('target_uri')) or _extract_enst_id(edge.get('source_uri'))
            if enst_id:
                update_doc['_to'] = f"nodes/{enst_id}"
            else:
                missing_identifier = True

        if len(update_doc) == 1 or missing_identifier:
            skipped += 1
            if missing_identifier:
                print(f"  ⚠ Skipped edge {edge['_key']} - unable to determine ENST identifier")
            continue

        try:
            edges_collection.update(update_doc)
            updated += 1
        except Exception as e:
            skipped += 1
            print(f"  ⚠ Failed to update edge {edge['_key']}: {e}")

    print(f"✓ Fixed {updated}/{len(edges_to_fix)} edges (skipped: {skipped})")



# --- Main execution / CLI -------------------------------------------------

def _parse_args():
    p = argparse.ArgumentParser(
        description="Load the property-graph JSON files into ArangoDB.",
    )
    p.add_argument("--db", default=DEFAULT_DB_NAME,
                   help=f"ArangoDB database name (default: {DEFAULT_DB_NAME}).")
    p.add_argument("--host", default=arangodb_hosts,
                   help=f"ArangoDB host URL (default: {arangodb_hosts}).")
    p.add_argument("--user", default=arangodb_user,
                   help=f"ArangoDB user (default: {arangodb_user}).")
    p.add_argument("--password", default=arangodb_password,
                   help="ArangoDB password (default: same as the project default).")
    p.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR,
                   help="Directory containing nodes.json + edges.json.")
    p.add_argument("--graph-file", type=Path, default=None,
                   help="Single-file graph dump (with 'nodes' and 'edges' keys). "
                        "Overrides --input-dir.")
    p.add_argument("--nodes-collection", default="nodes",
                   help="Vertex collection name (default: nodes).")
    p.add_argument("--edges-collection", default="edges",
                   help="Edge collection name (default: edges).")
    p.add_argument("--no-fix-summary-edges", action="store_true",
                   help="Skip the post-import fix of edges referencing the nodes/Summary placeholder (default: fix is applied).")
    p.add_argument("--no-schema-report", action="store_true",
                   help="Skip the final KG schema report.")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    print("=" * 60)
    print("load_graph_to_arangodb")
    print(f"  host      : {args.host}")
    print(f"  user      : {args.user}")
    print(f"  database  : {args.db}")
    if args.graph_file:
        print(f"  graph file: {args.graph_file}")
    else:
        print(f"  input dir : {args.input_dir}")
    print("=" * 60)

    # Override arangodb_utils module globals with the CLI values so that
    # setup_arangodb_connection() (which reads them) honours --host/--user/--password.
    import arangodb_utils
    arangodb_utils.arangodb_hosts = args.host
    arangodb_utils.arangodb_user = args.user
    arangodb_utils.arangodb_password = args.password
    from arangodb_utils import setup_arangodb_connection
    # This is a loader: it legitimately populates a fresh instance, so it is the
    # one place allowed to create the database if it is not there yet.
    db = setup_arangodb_connection(args.db, create=True)
    if db is None:
        sys.exit("[ERROR] Failed to connect to ArangoDB. Is the server running?")

    create_arangodb_collections(
        db,
        nodes_collection=args.nodes_collection,
        edges_collection=args.edges_collection,
    )

    if args.graph_file:
        import_data_to_arangodb(
            db,
            nodes_collection=args.nodes_collection,
            edges_collection=args.edges_collection,
            graph_file=str(args.graph_file),
        )
    else:
        import_data_to_arangodb(
            db,
            nodes_collection=args.nodes_collection,
            edges_collection=args.edges_collection,
            nodes_file=str(args.input_dir / "nodes.json"),
            edges_file=str(args.input_dir / "edges.json"),
        )

    if not args.no_fix_summary_edges:
        fix_summary_edges(db)

    print("\nFinal verification:")
    get_collections_data(db)

    print("\n" + "=" * 60)
    print("[OK] load_graph_to_arangodb finished.")
    print("=" * 60)



"""

Knowledge Graph in PKT_main:
- collection: nodes
  properties: _key, uri, namespace, entity_id, class_code, label, bioentity_type, description, synonym, source, source_type, integer_id, entrez_id, uniprot_id, rsid, mondo_id, hp_id

- collection: edges
  properties: _from, _to, edge_id, source_uri, target_uri, predicate_uri, predicate_label, predicate_class_code, predicate_bioentity_type, predicate_source

"""

# =============================================================================
# KG SCHEMA ANALYSIS - Property Statistics
# =============================================================================

def get_collection_schema_stats(db_connection, collection_name: str, sample_size: int = 1000) -> dict:
    """
    Get comprehensive schema statistics for a collection.
    Uses sampling for property discovery and individual counts for efficiency.
    """
    # Get total document count
    total_query = f"RETURN LENGTH({collection_name})"
    total_count = list(db_connection.aql.execute(total_query))[0]
    
    # Get property keys from a sample (much faster than scanning all docs)
    sample_aql = f"""
    LET sample = (
        FOR doc IN {collection_name}
            LIMIT {sample_size}
            RETURN ATTRIBUTES(doc, true)
    )
    RETURN UNIQUE(FLATTEN(sample))
    """
    property_keys = list(db_connection.aql.execute(sample_aql))[0]
    
    # Count non-null values for each property individually
    results = []
    for prop in property_keys:
        count_aql = f"""
        RETURN LENGTH(
            FOR doc IN {collection_name}
                FILTER doc.@prop != null AND doc.@prop != ""
                RETURN 1
        )
        """
        non_null_count = list(db_connection.aql.execute(count_aql, bind_vars={"prop": prop}))[0]
        coverage_pct = round(100 * non_null_count / total_count, 2) if total_count > 0 else 0
        results.append({
            "property": prop,
            "non_null_count": non_null_count,
            "null_count": total_count - non_null_count,
            "coverage_pct": coverage_pct
        })
    
    # Sort by non_null_count descending
    results.sort(key=lambda x: x["non_null_count"], reverse=True)
    
    return {
        "collection": collection_name,
        "total_documents": total_count,
        "properties": results
    }


def get_property_value_distribution(db_connection, 
                                    collection_name: str,
                                    property_name: str, 
                                    limit: int = 20) -> dict:
    """
    Get value distribution for a specific property (useful for categorical fields).
    """
    aql = f"""
    FOR doc IN {collection_name}
        FILTER doc.@prop != null AND doc.@prop != ""
        COLLECT value = doc.@prop WITH COUNT INTO count
        SORT count DESC
        LIMIT @limit
        RETURN {{value: value, count: count}}
    """
    results = list(db_connection.aql.execute(aql, bind_vars={"prop": property_name, "limit": limit}))
    return {
        "property": property_name,
        "top_values": results
    }


def analyze_kg_schema(db_connection, verbose: bool = True, skip_nodes: bool = False) -> dict:
    """
    Complete KG schema analysis for nodes and edges collections.
    Returns schema stats and value distributions for key categorical properties.
    """
    report = {}
    
    # Analyze nodes collection
    if not skip_nodes:
        print("=" * 60)
        print("ANALYZING NODES COLLECTION")
        print("=" * 60)
        nodes_stats = get_collection_schema_stats(db_connection, "nodes")
        report["nodes"] = nodes_stats
        
        if verbose:
            print(f"\nTotal nodes: {nodes_stats['total_documents']}")
            print(f"\n{'Property':<25} {'Non-Null':<12} {'Null':<12} {'Coverage %':<10}")
            print("-" * 60)
            for p in nodes_stats["properties"]:
                print(f"{p['property']:<25} {p['non_null_count']:<12} {p['null_count']:<12} {p['coverage_pct']:<10}")
        
        # Key categorical properties for nodes
        node_categorical = ["class_code", "bioentity_type", "source", "source_type", "namespace"]
        report["nodes_distributions"] = {}
        
        for prop in node_categorical:
            dist = get_property_value_distribution(db_connection, "nodes", prop)
            report["nodes_distributions"][prop] = dist
            if verbose:
                print(f"\n--- {prop} value distribution (top 20) ---")
                for v in dist["top_values"]:
                    print(f"  {v['value']}: {v['count']}")
        
    # Analyze edges collection
    print("\n" + "=" * 60)
    print("ANALYZING EDGES COLLECTION")
    print("=" * 60)
    edges_stats = get_collection_schema_stats(db_connection, "edges")
    report["edges"] = edges_stats
    
    if verbose:
        print(f"\nTotal edges: {edges_stats['total_documents']}")
        print(f"\n{'Property':<25} {'Non-Null':<12} {'Null':<12} {'Coverage %':<10}")
        print("-" * 60)
        for p in edges_stats["properties"]:
            print(f"{p['property']:<25} {p['non_null_count']:<12} {p['null_count']:<12} {p['coverage_pct']:<10}")
    
    # Key categorical properties for edges
    edge_categorical = ["predicate_label", "predicate_class_code", "predicate_bioentity_type", "predicate_source"]
    report["edges_distributions"] = {}
    
    for prop in edge_categorical:
        dist = get_property_value_distribution(db_connection, "edges", prop, limit = None)  # Get full distribution for edge properties
        report["edges_distributions"][prop] = dist
        if verbose:
            print(f"\n--- {prop} value distribution ---")
            for v in dist["top_values"]:
                print(f"  {v['value']}: {v['count']}")
    
    print("\n" + "=" * 60)
    print("SCHEMA ANALYSIS COMPLETE")
    print("=" * 60)
    
    return report


# --- Generic class-code schema extraction ---------------------------------

def extract_generic_schema(db_connection,
                           nodes_collection='nodes',
                           edges_collection='edges',
                           exclude_predicate_class_codes=None,
                           as_dataframe=False,
                           verbose=True):
    """
    Extract the generic schema of the KG by aggregating triples on the
    class_code of the source/target nodes and the predicate_class_code
    of the edges.

    Implementation note: instead of calling DOCUMENT(edge._from) for every
    edge (slow), this query uses an explicit FOR...FILTER on _id, which lets
    ArangoDB use its primary index for the join.

    Returns a list of dicts (or a DataFrame if as_dataframe=True) with the
    columns: source_class_code, predicate_class_code, target_class_code, count.
    """
    if exclude_predicate_class_codes is None:
        exclude_predicate_class_codes = []

    aql = f"""
    FOR edge IN {edges_collection}
        FILTER LENGTH(@exclude_codes) == 0
            OR edge.predicate_class_code NOT IN @exclude_codes
        FOR source IN {nodes_collection}
            FILTER source._id == edge._from
            FOR target IN {nodes_collection}
                FILTER target._id == edge._to
                COLLECT
                    src_cc  = source.class_code,
                    pred_cc = edge.predicate_class_code,
                    tgt_cc  = target.class_code
                WITH COUNT INTO count
                SORT count DESC
                RETURN {{
                    "source_class_code":    src_cc,
                    "predicate_class_code": pred_cc,
                    "target_class_code":    tgt_cc,
                    "count":                count
                }}
    """

    bind_vars = {"exclude_codes": exclude_predicate_class_codes}

    if verbose:
        print("[run] Executing generic-schema query (explicit JOIN) ...")

    results = list(db_connection.aql.execute(
        aql,
        bind_vars=bind_vars,
        max_runtime=540  # server-side cap (under the 600s client timeout)
    ))

    if verbose:
        print(f"\n[OK] Extracted generic schema (class_code): {len(results)} distinct triples")
        print(f"\n{'Source':>20}  {'Predicate':>15}  {'Target':>20}  {'Count':>8}")
        print("-" * 70)
        for r in results[:20]:
            print(
                f"  {r['source_class_code']:>20}"
                f"  [{r['predicate_class_code']:>13}]"
                f"  {r['target_class_code']:>20}"
                f"  (n={r['count']})"
            )
        if len(results) > 20:
            print(f"  ... ({len(results) - 20} more triples)")

    if as_dataframe:
        try:
            import pandas as pd
            return pd.DataFrame(results)
        except ImportError:
            print("[warn] pandas not installed; returning list of dicts.")
            return results

    return results




# --- Direct-from-JSON schema extraction (no ArangoDB required) -----------

def extract_generic_schema_from_files(
        nodes_file,
        edges_file,
        exclude_predicate_class_codes=None,
        include_predicate_class_codes=None,
        as_dataframe=False,
        verbose=True):
    """
    Extract the generic class-code schema directly from the local nodes.json
    and edges.json files (faster than going through ArangoDB).

    Args:
        nodes_file: path to nodes.json
        edges_file: path to edges.json
        exclude_predicate_class_codes: predicate_class_codes to skip (e.g. ["RDF"])
        include_predicate_class_codes: if non-empty, only keep these predicate_class_codes
        as_dataframe: return a pandas DataFrame instead of a list of dicts
        verbose: print progress and a summary

    Returns the same data shape as extract_generic_schema().
    """
    if exclude_predicate_class_codes is None:
        exclude_predicate_class_codes = []
    if include_predicate_class_codes is None:
        include_predicate_class_codes = []
    exclude_set = set(exclude_predicate_class_codes)
    include_set = set(include_predicate_class_codes)

    # Step 1: load nodes and build a URI -> class_code lookup
    if verbose:
        print(f"[load] reading nodes from {nodes_file}...")

    with open(nodes_file, 'r', encoding='utf-8') as f:
        nodes_data = json.load(f)

    uri_to_class_code = {}
    for node in nodes_data:
        uri = node.get('uri')
        class_code = node.get('class_code')
        if uri and class_code:
            uri_to_class_code[uri] = class_code

    if verbose:
        print(f"[OK] {len(nodes_data)} nodes loaded -> {len(uri_to_class_code)} URIs indexed")

    # Step 2: load edges and aggregate (src_class, pred_class, tgt_class) triples
    if verbose:
        print(f"[load] reading edges from {edges_file}...")

    with open(edges_file, 'r', encoding='utf-8') as f:
        edges_data = json.load(f)

    if verbose:
        print(f"[OK] {len(edges_data)} edges loaded")
        print("[run] aggregating schema triples ...")

    triple_counts = defaultdict(int)
    skipped = 0

    for edge in edges_data:
        pred_cc = edge.get('predicate_class_code')

        if pred_cc in exclude_set:
            continue
        if include_set and pred_cc not in include_set:
            continue

        src_uri = edge.get('source_uri')
        tgt_uri = edge.get('target_uri')
        src_cc = uri_to_class_code.get(src_uri)
        tgt_cc = uri_to_class_code.get(tgt_uri)

        if not src_cc or not tgt_cc:
            skipped += 1
            continue

        triple_counts[(src_cc, pred_cc, tgt_cc)] += 1

    if verbose and skipped > 0:
        print(f"[warn] {skipped} edges skipped (src or tgt node missing from nodes.json)")

    # Step 3: build the result list, sorted by descending count
    results = [
        {
            "source_class_code":    src_cc,
            "predicate_class_code": pred_cc,
            "target_class_code":    tgt_cc,
            "count":                count
        }
        for (src_cc, pred_cc, tgt_cc), count
        in sorted(triple_counts.items(), key=lambda x: -x[1])
    ]

    if verbose:
        print(f"\n[OK] generic schema extracted: {len(results)} distinct triples\n")
        print(f"  {'Source':<20}  {'Predicate':<15}  {'Target':<20}  {'Count':>8}")
        print("  " + "-" * 68)
        for r in results[:30]:
            print(
                f"  {r['source_class_code']:<20}"
                f"  {r['predicate_class_code']:<15}"
                f"  {r['target_class_code']:<20}"
                f"  {r['count']:>8}"
            )
        if len(results) > 30:
            print(f"  ... ({len(results) - 30} more triples)")

    if as_dataframe:
        try:
            import pandas as pd
            return pd.DataFrame(results)
        except ImportError:
            print("[warn] pandas not installed; returning list of dicts.")
            return results

    return results
