
"""
### **UC1 — Semantic proximity predicts cross-layer quantitative coherence** *(il più originale)*

**Idea:** Se due geni sono connessi nel KG da una relazione semantica specifica (es. "causally influences", "molecularly interacts with"), i loro profili quantitativi attraverso i layer omici sono più coerenti rispetto a coppie di geni non connesse?

**Come si fa con il tuo framework:**
```
1. Estrai coppie di geni connesse nel KG per tipo di relazione 
   (1-hop: "causally influences" / "molecularly interacts with" / "has participant")
2. Per ogni coppia, calcola correlazione mRNA-mRNA, protein-protein, 
   e cross-layer mRNA-protein su campioni BRCA
3. Confronta: coppie semanticamente connesse vs coppie random
4. Stratifica per tipo di relazione: "causally influences" 
   dovrebbe dare correlazioni più alte di "located_in"
```

**Perché è pubblicabile:** dimostra quantitativamente che il layer semantico KG è un prior biologicamente valido per l'analisi multi-omica. Non è mai stato fatto con questa architettura. Nessun paper TCGA multi-omics ha mai usato le relazioni ontologiche come variabile indipendente per predire la coerenza quantitativa cross-layer.

**Visualizzazioni:**
- Boxplot: distribuzione delle correlazioni per tipo di relazione semantica vs random
- Network: subgrafo KG con edges colorati per forza di correlazione mRNA-protein
- Scatter: distanza semantica (numero di hop nel KG) vs correlazione cross-layer
"""


"""
Source di codice:

codice che costruscie la rete trasomica per paziente:
scripts/build_transomic_network.py

analyzer:
scripts/analyzie_kg_transomics.py

utils:
scripts/omics_query_utils.py
scripts/query_utils.py
scripts/kg_transomics_aql_query.md



"""