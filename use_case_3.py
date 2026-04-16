"""

### **UC3 — Phenotype-anchored trans-omic subgraph extraction** *(il più adatto alle "visualizzazioni" del supervisor)*

**Idea:** Parti da un nodo fenotipo nel KG (HP ontology) biologicamente rilevante, traversa il grafo per costruire un gene set semanticamente grounded, poi estrai e visualizza il profilo multi-omico di quel gene set su campioni TCGA.

**Esempio concreto con i tuoi dati:**
- Nodo seed: un fenotipo HP legato a BRCA (es. "HP:0003002" = Breast carcinoma, o fenotipi di risposta cellulare)
- Traversal 2-hop: fenotipo → geni (via "has phenotype" / "causally influenced by") → proteine (via "gene product of")
- Gene set risultante: 50-200 geni semanticamente grounded al fenotipo
- Estrazione quantitativa: espressione TPM + CNV + beta metilazione + abbondanza proteica per tutti i campioni BRCA
- Visualizzazione: una "trans-omic map" del fenotipo

**Cosa mostra che è nuovo:** lo stesso gene set definito semanticamente mostra pattern multi-omici coerenti e interpretabili. Un gene set definito puramente statisticamente (es. top varianza) non ha questa interpretabilità.

**Visualizzazioni (questo è il cuore per il supervisor):**
- **Subgraph semantico** (Cytoscape/D3): fenotipo → geni → proteine → GO terms, con nodi colorati per espressione media in BRCA
- **Heatmap multi-omica**: geni (righe) × sample BRCA (colonne), con barre di annotazione per layer (TPM / CNV / beta / RPPA)
- **Radar/spider chart**: profilo multi-omico medio per il gene set, per ogni layer
- **Sankey diagram**: fenotipo → geni → pathway → layer omici, con spessore proporzionale alla forza del segnale quantitativo
"""


"""
Source di codice:

codice che costruisce la rete trasomica per paziente:
scripts/build_transomic_network.py

analyzer:
scripts/analyzie_kg_transomics.py

utils:
scripts/omics_query_utils.py
scripts/query_utils.py
scripts/kg_transomics_aql_query.md


Folder:
- transomic-networks: contiene i JSON dei grafi trasomici per ogni paziente

Referece:
- for the structure of the script UC3, you can refer to use_case_1.py, which has a similar structure but focuses on a different biological question. 


Studia source di codice e scrivi codice per realizzare UC3, con visualizzazioni.
"""
