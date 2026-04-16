"""### **UC2 — CNV-expression-protein discordance con spiegazione semantica** *(il più biologicamente ricco)*

La discordanza tra amplificazione del DNA, abbondanza di mRNA e abbondanza proteica è un fenomeno documentato: nel caso di SLPI in BRCA si osserva forte amplificazione CNV ma mRNA e proteina molto bassi, suggerendo dosage compensation. Questo fenomeno è noto, ma non è mai stato esplorato con contesto semantico.

**Idea:** per i geni con dati RPPA disponibili in BRCA, identifica quelli "concordanti" (CNV→mRNA→proteina tutti correlati) e quelli "discordanti" (CNV alto, proteina bassa). Poi interroga il KG per capire se i geni discordanti condividono proprietà semantiche specifiche — GO terms, pathway, tipo di interazione — che spiegano la compensazione.

**Come si fa:**
```
1. Per ogni gene nell'indice RPPA (487 proteine):
   - Calcola correlazione CNV→mRNA su campioni BRCA
   - Calcola correlazione mRNA→protein su campioni BRCA
2. Classifica: concordant (entrambe alte) / discordant (CNV-mRNA alta, mRNA-prot bassa)
3. Per i geni discordanti: query AQL → nodi KG → 
   quali GO terms, pathway, predicati li caratterizzano?
4. Risultato: i geni discordanti sono arricchiti in GO terms 
   legati a traduzione? regolazione post-traduzionale? 
   Sono più spesso "molecularly interacts with" ubiquitin ligases?
```

**Perché è pubblicabile:** la domanda biologica è genuinamente interessante (dosage compensation + semantic context) e la risposta non è banale. Il framework è l'unico modo per fare questa analisi in modo sistematico su tutti i geni con coverage multi-layer.

**Visualizzazioni:**
- Scatter plot 2D: correlazione CNV-mRNA (asse x) vs correlazione mRNA-protein (asse y), quadranti colorati
- Network: geni discordanti nel KG, colorati per tipo di GO process
- Heatmap: geni × layer omici per il subset discordante vs concordante"""