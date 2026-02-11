#%%
"""
Analyzer for transomic property-graph JSON files produced by build_transomic_network.py.

Functions:
- load_graph(path): read a single graph JSON.
- analyze_transomic_graph(graph): compute stats (counts, predicates, omic value summaries).
- top_features_by_abs(nodes, omic_type, n): extract top-N |value| features for an omic.

CLI:
	python step5-analyzie-kg-transomics_v0.py --folder ../temp_dir --sample TCGA-XX-YYYY

Notebook:
	from step5-analyzie-kg-transomics_v0 import load_graph, analyze_transomic_graph
	g = load_graph("../temp_dir/transomic_graph_TCGA-BRCA_TCGA-XX-YYYY.json")
	summary = analyze_transomic_graph(g)
"""

#%% Imports
import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

from matplotlib.pylab import rint


#%% Core helpers
def load_graph(path: str) -> Dict:
	with open(path, "r", encoding="utf-8") as f:
		return json.load(f)


def _numeric_stats(values: List[float]) -> Dict[str, Optional[float]]:
	vals = [v for v in values if isinstance(v, (int, float))]
	if not vals:
		return {"count": 0, "min": None, "max": None, "mean": None}
	return {
		"count": len(vals),
		"min": min(vals),
		"max": max(vals),
		"mean": sum(vals) / len(vals),
	}


def analyze_transomic_graph(graph: Dict) -> Dict:
	nodes = graph.get("nodes", []) or []
	edges = graph.get("edges", []) or []

	# Node-level summaries
	omic_counts = Counter(n.get("omic_type") for n in nodes)
	kg_linked = sum(1 for n in nodes if n.get("kg_node_key"))

	value_stats = {}
	values_by_omic = defaultdict(list)
	for n in nodes:
		if isinstance(n.get("value"), (int, float)):
			values_by_omic[n.get("omic_type")].append(n["value"])
	for omic, vals in values_by_omic.items():
		value_stats[omic] = _numeric_stats(vals)

	# Edge-level summaries
	predicate_counts = Counter()
	neighbor_bio_types = Counter()
	for e in edges:
		edge_doc = e.get("edge", {}) or {}
		predicate_counts[edge_doc.get("predicate_label")] += 1
		neighbor = e.get("neighbor", {}) or {}
		bio = neighbor.get("bioentity_type")
		if bio:
			neighbor_bio_types[bio] += 1

	# Metadata
	meta = graph.get("metadata", {})

	return {
		"metadata": meta,
		"node_count": len(nodes),
		"edge_count": len(edges),
		"omic_counts": dict(omic_counts),
		"kg_linked_nodes": kg_linked,
		"value_stats": value_stats,
		"predicate_counts": dict(predicate_counts),
		"neighbor_bioentity_types": dict(neighbor_bio_types),
	}


def top_features_by_abs(nodes: List[Dict], omic_type: str, n: int = 10) -> List[Tuple[str, float]]:
	feats = []
	for node in nodes:
		if node.get("omic_type") != omic_type:
			continue
		val = node.get("value")
		if isinstance(val, (int, float)):
			feats.append((node.get("id"), val))
	feats.sort(key=lambda x: abs(x[1]), reverse=True)
	return feats[:n]

#%% Notebook version of main analysis 

def analyze_folder(
	folder: str = "..\\transomic-networks",
	sample_substr: Optional[str] = None,
	print_top: int = 0,
) -> List[Dict]:
	"""
	Notebook-friendly wrapper: analyzes all graph JSON in folder (optionally filtered).

	Returns a list of summary dicts (one per file). Optionally prints top-N features per omic.
	"""
	files = [f for f in os.listdir(folder) if f.endswith(".json")]
	if sample_substr:
		files = [f for f in files if sample_substr in f]

	summaries = []
	for fname in sorted(files):
		path = os.path.join(folder, fname)
		graph = load_graph(path)
		summary = analyze_transomic_graph(graph)
		summary["file"] = fname
		summaries.append(summary)

		if print_top > 0:
			nodes = graph.get("nodes", []) or []
			omics = sorted({n.get("omic_type") for n in nodes})
			print(f"\n==== {fname} ====")
			print("Sample:", summary["metadata"].get("sample_id"))
			for omic in omics:
				tops = top_features_by_abs(nodes, omic, n=print_top)
				print(f"Top {print_top} {omic} by |value|:", tops)

	return summaries


def analyze_graph_file(path: str, print_top: int = 0) -> Dict:
	"""Analyze a single graph JSON file (no folder sweep)."""
	graph = load_graph(path)
	summary = analyze_transomic_graph(graph)

	# Print compact summary
	print(f"\n==== {os.path.basename(path)} ====")
	print("Sample:", summary["metadata"].get("sample_id"))
	print("Cohort:", summary["metadata"].get("cohort"))
	print("Nodes:", summary["node_count"], "(KG-linked:", summary["kg_linked_nodes"], ")")
	print("Edges:", summary["edge_count"])
	print("Omic counts:")
	for omic, count in summary["omic_counts"].items():
		print(f"  - {omic}: {count}")
	top = 200
	print(f"Predicates (top {top}):")
	for pred, count in Counter(summary["predicate_counts"]).most_common(top):
		print(f"  - {pred}: {count}")

	# if print_top > 0:
	# 	nodes = graph.get("nodes", []) or []
	# 	omics = sorted({n.get("omic_type") for n in nodes})
	# 	for omic in omics:
	# 		tops = top_features_by_abs(nodes, omic, n=print_top)
	# 		print(f"Top {print_top} {omic} by |value|:", tops)


	return summary

# TCGA-GM-A2DD
Case = "TCGA-BH-A18U"
filename =f"transomic_graph_TCGA-BRCA_{Case}-01A.json"
# filename = "transomic_graph_TCGA-BH-A18U-01A_edgelimit-none.json"
summary = analyze_graph_file(f"../transomic-networks\\{filename}", print_top=5)
#%%
summary["predicate_counts"]
"""
Predicates (no edge limit):
  - transcribed from: 176705
  - transcribed to: 176705
  - causally influences: 143789
  - causally influenced by: 143789
  - has participant: 104055
  - participates in: 104055
  - type: 69062
  - causes or contributes to condition: 37152
  - interacts with: 32662
  - only_in_taxon: 19540
  - has_gene_template: 19490
  - gene product of: 19172
  - has gene product: 19172
  - genetically interacts with: 6034
  - disease has basis in dysfunction of: 4471
  - has material basis in gain of function germline mutation in: 4
  - disease has basis in disruption of: 2
  - disease causes dysfunction of: 1

Predicates (edge limit 5000):
  - transcribed from: 1141
  - transcribed to: 1141
  - type: 1069
  - has participant: 312
  - participates in: 312
  - interacts with: 264
  - causally influences: 135
  - causally influenced by: 135
  - causes or contributes to condition: 123
  - has_gene_template: 87
  - gene product of: 86
  - has gene product: 85
  - only_in_taxon: 85
  - disease has basis in dysfunction of: 13
  - genetically interacts with: 12

"""
#%% CLI
def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Analyze transomic graph JSON files")
	parser.add_argument("--file", help="Path to a single graph JSON to analyze")
	parser.add_argument("--folder", default="..\\transomic-networks", help="Folder containing graph JSON files (used if --file is not provided)")
	parser.add_argument("--sample", help="Substring to filter filenames inside folder; the first match is used")
	parser.add_argument("--print-top", type=int, default=0, help="If >0, print top-N features per omic")
	return parser.parse_args()


def main():
	args = parse_args()

	# Prefer a single explicit file if provided
	if args.file:
		if not os.path.isfile(args.file):
			print("File not found:", args.file)
			return
		summary = analyze_graph_file(args.file, print_top=args.print_top)
		print("\n====", os.path.basename(args.file), "====")
		print("Sample:", summary["metadata"].get("sample_id"))
		print("Cohort:", summary["metadata"].get("cohort"))
		print("Nodes:", summary["node_count"], "(KG-linked:", summary["kg_linked_nodes"], ")")
		print("Edges:", summary["edge_count"])
		print("Omic counts:", summary["omic_counts"])
		print("Predicates (top 5):", Counter(summary["predicate_counts"]).most_common(5))
		return

	# Otherwise pick a single file from folder (optionally filtered)
	folder = args.folder
	files = [f for f in os.listdir(folder) if f.endswith(".json")]
	if args.sample:
		files = [f for f in files if args.sample in f]

	if not files:
		print("No graph files found.")
		return

	fname = sorted(files)[0]
	path = os.path.join(folder, fname)
	summary = analyze_graph_file(path, print_top=args.print_top)
	print("\n====", fname, "====")
	print("Sample:", summary["metadata"].get("sample_id"))
	print("Cohort:", summary["metadata"].get("cohort"))
	print("Nodes:", summary["node_count"], "(KG-linked:", summary["kg_linked_nodes"], ")")
	print("Edges:", summary["edge_count"])
	print("Omic counts:", summary["omic_counts"])
	print("Predicates (top 5):", Counter(summary["predicate_counts"]).most_common(5))


if __name__ == "__main__":
	# Skip argparse when running inside notebooks/VS Code interactive to avoid errors.
	if "ipykernel" not in sys.modules:
		main()

