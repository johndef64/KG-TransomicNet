


KG di riferimento utilizzati includono:
- **PheKnowLator**: per associazioni tra geni, malattie e fenotipi, con oltre 10 milioni di triple, basato su ontologie OWL/RDF e URI semantici, formato come RDF/OWL. Supporta rappresentazioni flessibili (simple, hybrid, complex graphs) e reasoning automatico;
- **Monarch**: per relazioni cross-species tra geni, malattie e fenotipi, integra 33 fonti di dati biomedici, utilizza il modello Biolink e PHENIO come semantic layer, formato come KGX/RDF;
- **Clinical Knowledge Graph (CKG)**: per associazioni tra geni, proteine, malattie e farmaci; incentrato sulla proteomica clinica; formato Property Graph (Neo4J) con 20 milioni di nodi e relazioni;
- **Metabolomics Knowledge Graph (MetaKG)**: framework comprehensivo per l'integrazione di dati metabolomici da HMDB, SMPDB e KEGG in una rappresentazione unificata. Supporta machine learning con KG embeddings (TransE, RotatE, ComplEx, ConvE), analisi statistiche, ricerca di pathway e predizione di link. Include API per estrazione dati e strumenti di visualizzazione avanzati; 



Benchmark tra KG di riferimento:
| Knowledge Graph            | Entity Types                        | Relation Types         | Triple Count           |
|---------------------------|-------------------------------------|-----------------------|------------------------|
| PheKnowLator              | Gene, Disease, Phenotype, Protein, Pathway, Variant, Chemical, Anatomy, Tissue, Cell | 100+                  | 1,000,000+ (fino a 10 milioni+ in alcune versioni) |
| Monarch                   | Gene, Disease, Phenotype, Pathway, Anatomy, Chemical | 50+                    | 2,000,000+             |
| Clinical Knowledge Graph  | 36                                  | 47                    | ~20 milioni (triples)  |
| Metabolomics Knowledge Graph | Metabolite, Enzyme, Pathway, Disease, Protein | 50+                   | 1,000,000+ (integrazione HMDB, SMPDB, KEGG) |


<!-- riesci ad estrami anche le ontologie per nodi ed archi che sono rappresentate nel CKG? -->

# Ontolgies Usage

## Clinical Knowledge Graph (CKG)

Some of the nodes in our graph data model could not be described using ontologies or existing terminologies, and they needed to be standardized using identifiers from the selected biomedical databases

**Per i nodi (entità):**
- UniProt (proteine)
- HMDB (metaboliti)  
- DrugBank (farmaci)
- SNOMED-CT (termini clinici)
- Gene Ontology (GO) per annotazioni funzionali
- Altre ontologie specifiche per tipo di entità (es. malattie, tessuti)

**Per gli archi (relazioni):**
- Relation Ontology (RO) per relazioni generali tra entità biologiche (es. part_of, regulates, derives_from)
- Altre ontologie e terminologie standard per tipi di relazione specifici (es. HAS_PARENT, HAS_QUANTIFIED_PROTEIN)
- SNOMED-CT per relazioni cliniche

## PheKnowLator

**Per i nodi (entità):**
- Gene Ontology (GO) per annotazioni funzionali di geni e proteine
- Human Phenotype Ontology (HPO) per fenotipi umani
- Mondo Disease Ontology (MONDO) per malattie
- Chemical Entities of Biological Interest (CHEBI) per entità chimiche
- Protein Ontology (PRO) per proteine e loro forme
- Sequence Ontology (SO) per sequenze genomiche
- Cell Ontology (CL) per tipi cellulari
- Uberon per anatomia comparativa
- OBO foundational ontologies

**Per gli archi (relazioni):**
- Relation Ontology (RO) per relazioni fondamentali (part_of, regulates, etc.)
- Web Ontology Language (OWL) expressions per reasoning complesso
- RDF/RDFS per triple semantiche
- Biolink Model compatibile per standardizzazione

## Monarch

**Per i nodi (entità):**
- Biolink Model come schema principale per entità biologiche
- PHENIO (Phenomics Integrated Ontology) come semantic layer unificante
- Gene Ontology (GO) per funzioni geniche (BP, MF, CC)
- Human Phenotype Ontology (HPO) per fenotipi umani
- Mondo Disease Ontology (MONDO) per malattie
- Chemical Entities of Biological Interest (CHEBI) per entità chimiche
- Anatomical ontologies (UBERON, CL) per anatomia e cellule
- Species-specific ontologies per organismi modello

**Per gli archi (relazioni):**
- Biolink predicates standardizzati
- RO (Relation Ontology) per relazioni biologiche fondamentali  
- Cross-species orthology relations
- Phenotype similarity measures
- Disease-gene-phenotype association predicates




## MetaKG (Metabolomics Knowledge Graph)

**Per i nodi (entità):**
- HMDB (Human Metabolome Database) identifiers per metaboliti umani
- SMPDB (Small Molecule Pathway Database) identifiers per pathway e proteine
- KEGG (Kyoto Encyclopedia of Genes and Genomes) identifiers per pathway e enzimi
- Chemical Entities of Biological Interest (CHEBI) per classificazione chimica
- Custom entity schemas per malattie, proteine ed enzimi
- Cross-database entity alignment e mapping automatico

**Per gli archi (relazioni):**
- has_pathway (metabolita → pathway)
- has_disease (metabolita → malattia)  
- has_protein (pathway → proteina)
- involved_in (metabolita → processo biologico)
- interacts_with (metabolita ↔ metabolita)
- Custom relation types per integrazione multi-database
- Machine learning derived relationships via KG embeddings

**Caratteristiche tecniche:**
- Schema unificato per integrazione HMDB, SMPDB, KEGG
- Supporto KG embedding models (TransE/D/H/R, RotatE, ComplEx, ConvE)
- Pipeline ML per link prediction e scoperta di nuove associazioni
- API per estrazione dati e query flessibili
- Strumenti di visualizzazione (Sankey diagrams, network plots)
- Sistema di analisi statistiche e enrichment analysis

## Sintesi Comparativa

| KG | Ontologie Nodi | Ontologie Archi | Modello |
|---|---|---|---|
| CKG | UniProt, HMDB, DrugBank, SNOMED-CT, GO | RO, SNOMED-CT, custom relations | Property Graph |
| PheKnowLator | GO, HPO, MONDO, CHEBI, PRO, SO, CL, Uberon | RO, OWL expressions, RDF/RDFS | RDF/OWL |
| Monarch | Biolink, PHENIO, GO, HPO, MONDO, CHEBI, Uberon | Biolink predicates, RO, orthology | KGX/RDF |
| MetaKG | HMDB, SMPDB, KEGG, CHEBI, custom schemas | Custom relations, ML-derived, pathway relations | Multi-format (XML, JSON, CSV) |

