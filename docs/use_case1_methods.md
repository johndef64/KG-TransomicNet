***

# UC1 Architecture

## Pipeline (6 steps)

1. **Find samples** con complete omics (expression + protein) from TCGA-BRCA Primary Tumor
2. **Extract gene pairs by predicate**  
   Queries ArangoDB per EntrezID gene pairs connessi da 7 semantic relationship types:  
   `causally influences`, `molecularly interacts with`, `interacts with`,  
   `genetically interacts with`, `has participant`, `has phenotype`, `has gene product`
3. **Generate random baseline**  
   Size-matched random gene pairs from the expression index
4. **Compute correlations** per ogni pair across all BRCA samples:  
   - mRNA-mRNA (Pearson r di TPM vectors)  
   - Protein-Protein (Pearson r di RPPA abundance)  
   - Cross-layer mRNA-Protein (mRNA gene A vs Protein gene B)
5. **Statistical testing**  
   Mann-Whitney U test (one-sided) di ogni predicate group vs random, con effect size
6. **Hop-distance analysis**  
   Correlations per gene pairs a 1, 2, 3 hops nel KG (testa l'ipotesi del gradient)

## 5 Visualizzazioni

| # | Tipo      | File                        | Descrizione                              |
|---|-----------|-----------------------------|------------------------------------------|
| 1 | Boxplot   | `uc1_boxplot_mrna.png`      | mRNA correlation distributions by predicate vs random |
| 2 | Boxplot   | `uc1_boxplot_protein.png`   | Protein correlation distributions by predicate vs random |
| 3 | Boxplot   | `uc1_boxplot_cross_layer.png` | Cross-layer mRNA-Protein correlations |
| 4 | Heatmap   | `uc1_heatmap_summary.png`   | Median correlations: predicate type x omic metric |
| 5 | Network   | `uc1_network_*.png`         | KG subgraph con edges colorati by correlation strength |
| 6 | Scatter   | `uc1_scatter_hop_distance.png` | Hop distance vs correlation (gradient plot) |

## How to run

### Opzione 1: CLI
```bash
conda activate gnn
python use_case_1.py
```

### Opzione 2: VS Code (interattivo)
- Usa cell mode con `#%%`
- Esegui solo l'ultima cell che chiama `run_uc1()`

**Risultati salvati in**: `results/uc1/` (CSVs + PNGs)

***

Cosa stiamo valutando
L'idea di UC1 è: le relazioni semantiche nel KG predicono la coerenza quantitativa tra geni?

Il test principale e sufficiente per dimostrarlo è:

mRNA-mRNA correlation cross-sample: prendo geneA e geneB connessi nel KG, e calcolo se quando geneA è sovraespresso in un paziente, anche geneB tende ad esserlo. Faccio questo su 878 pazienti BRCA. Se la correlazione è più alta delle coppie random, il KG "predice" la co-espressione.

Le correlazioni protein e cross-layer servono?
Sono un plus per il paper ma non sono il core:

Protein-protein: verifica che la coerenza si mantiene anche a livello proteomico (rinforza il risultato)
Cross-layer mRNA-protein: verifica che la relazione KG predice coerenza tra layer diversi (l'aspetto più originale del paper)
Ma il problema pratico è che RPPA (proteomica) copre pochi geni (~200 antibodies), quindi pochissime coppie hanno dati protein validi. I risultati lo confermano: n_valid_protein = 0 per quasi tutto.

Cosa funziona già
Il tuo primo run ha già il risultato chiave:

Gruppo	Mediana corr mRNA	p-value
genetically interacts with	0.144	2.6e-25
random	0.013	baseline
Questo è pubblicabile: le coppie di geni connesse nel KG da "genetically interacts with" hanno correlazione mRNA 10x più alta del random, con significatività altissima.



==========================
==================
==========================


***

# Valutazione dei risultati

## Conferma metodologica: ✅ SI

I dati confermano che l'approccio è **metodologicamente corretto**:

- **Random baseline**: mediana ~0.013 (quasi zero) — coerente con geni casuali non co-espressi
- **Campione robusto**: 878 pazienti, >350 coppie valide per gruppo — sufficiente per inferenza statistica
- **Coerenza biologica interna**: ordine dei predicati biologicamente sensato  
  `co-pathway (0.196) > molecularly interacts with (0.154) > genetically interacts with (0.144) > co-disease (0.067) > random (0.013)`

> **Spiegazione biologica**: I geni co-pathway hanno la coerenza più alta (co-regolazione trascrizionale). La co-malattia ha segnale più debole (meccanismi diversi).

## Significatività statistica: 🌟 ECCELLENTE

Tutti i 4 predicati sono **significativi vs random** con **p < 1e-10** (tutti `***`):

| Predicato                | Mediana r | p-value   | Pubblicabile?  |
|--------------------------|-----------|-----------|----------------|
| `co-pathway`             | 0.196     | 2e-38     | **molto forte** |
| `molecularly interacts with` | 0.154  | 1.4e-34   | **molto forte** |
| `genetically interacts with` | 0.144  | 2.6e-25   | **molto forte** |
| `co-disease`             | 0.067     | 8.7e-11   | **forte**      |
| `random (baseline)`      | 0.013     | —         | —              |

## Verdetto complessivo

**✅ Risultati significativi e pubblicabili**

**Finding principale**:  
Le relazioni semantiche nel Knowledge Graph predicono **significativamente** la co-espressione genica cross-paziente (tutti i predicati **p < 1e-10** vs random), con gradiente biologicamente coerente:  
`co-pathway > interazione molecolare > interazione genetica > co-malattia > random`

**Note tecniche**:
- Bug effect size **corretto** (ora usa rank-biserial correlation dall'U statistic)
- Limiti protein/cross-layer: da discutere come **limitazione copertura RPPA** nel paper
- **Rilancia lo script** per output finali con effect size corretto

***



======================================================================
  UC1: Semantic Proximity Predicts Cross-Layer Quantitative Coherence
======================================================================
🔌 Setting up ArangoDB connection to http://localhost:8529...
Database PKT_test10000 already exists
✓ Connected to database: PKT_test10000

[1/6] Finding samples with complete omics...
  Found 910 samples with both expression + protein data
  Of which 878 are Primary Tumor samples

[2/6] Extracting gene pairs from KG by predicate type...

  --- genetically interacts with (strategy: direct) ---
  Found 1174 unique gene-gene pairs
  Sampled down to 500 pairs
  Fetching expression matrix for 243 genes across 878 samples...
  Expression matrix: 878 samples x 203 genes
  Fetching protein abundance for 243 genes...
  Protein matrix: 878 samples x 1 genes
  Valid mRNA correlations: 354/500

  --- molecularly interacts with (strategy: via_protein) ---
    Building PR -> gene lookup map...
    PR->gene map: 19091 proteins mapped
    Fetched 4500 edges with PR endpoints
  Found 4137 unique gene-gene pairs
  Sampled down to 500 pairs
  Fetching expression matrix for 890 genes across 878 samples...
  Expression matrix: 878 samples x 885 genes
  Fetching protein abundance for 890 genes...
  Protein matrix: 878 samples x 45 genes
  Valid mRNA correlations: 490/500

  --- co-pathway (has participant) (strategy: co_pathway) ---
  Found 1497 unique gene-gene pairs
  Sampled down to 500 pairs
  Fetching expression matrix for 92 genes across 878 samples...
  Expression matrix: 878 samples x 92 genes
  Fetching protein abundance for 92 genes...
  Protein matrix: 878 samples x 39 genes
  Valid mRNA correlations: 484/500

  --- co-disease (strategy: co_disease) ---
  Found 1387 unique gene-gene pairs
  Sampled down to 500 pairs
  Fetching expression matrix for 249 genes across 878 samples...
  Expression matrix: 878 samples x 246 genes
  Fetching protein abundance for 249 genes...
  Protein matrix: 878 samples x 36 genes
  Valid mRNA correlations: 481/500

[3/6] Generating random gene pair baseline...
  Total genes in expression index: 41604
  Fetching expression matrix for 987 genes across 878 samples...
  Expression matrix: 878 samples x 987 genes
  Fetching protein abundance for 987 genes...
  Protein matrix: 878 samples x 5 genes
  Valid random mRNA correlations: 472/500

  Saved all results to G:\Altri computer\Horizon\horizon_workspace\projects\ActiveProjects\KG-TransomicNet\results\uc1/uc1_correlation_results.csv

[4/6] Computing summary statistics...
             predicate_group  n_pairs  n_valid_mrna  n_valid_protein  n_valid_cross  median_corr_mrna  mean_corr_mrna  std_corr_mrna  median_corr_protein  mean_corr_protein  median_corr_cross_ab  mean_corr_cross_ab
  genetically interacts with      500           354                0              9          0.143816        0.214139       0.252196                  NaN                NaN              0.036502            0.028939
  molecularly interacts with      500           490                2             22          0.153982        0.195436       0.233333             0.011994           0.011994             -0.012995            0.033922
co-pathway (has participant)      500           484               92            227          0.195583        0.179055       0.201222             0.030057           0.024959             -0.004036           -0.003148
                  co-disease      500           481               21             85          0.066533        0.106929       0.197649             0.041808           0.016756             -0.008862            0.001601
                      random      500           472                0              4          0.013383        0.030897       0.106437                  NaN                NaN              0.021939           -0.009473

======================================================================
  Mann-Whitney U test: each predicate vs random (corr_mrna)
======================================================================
  Predicate                               n    Median      U stat     p-value  Effect size r
  ----------------------------------- -----  --------  ----------  ----------  ------------
  random (baseline)                     472    0.0134
  genetically interacts with            354    0.1438      118592    2.63e-25      0.4195  ***
  molecularly interacts with            490    0.1540      168231    1.42e-34      0.4548  ***
  co-pathway (has participant)          484    0.1956      169316    2.05e-38      0.4823  ***
  co-disease                            481    0.0665      140633    8.70e-11      0.2389  ***
  Not enough random pairs for statistical testing.
  Not enough random pairs for statistical testing.