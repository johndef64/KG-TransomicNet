#!/usr/bin/env python3
"""
Script per creare dataset di mapping tra ENSG (Gene ID) e ENST (Transcript ID)
usando BioMart API tramite pybiomart.

Autore: Generated script
Data: 2025-12-09
"""

import pandas as pd
from pybiomart import Dataset
import time
from pathlib import Path
from typing import Optional, List
import logging

# Configurazione logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ENSGtoENSTMapper:
    """Classe per gestire il mapping ENSG -> ENST da Ensembl BioMart."""

    def __init__(self, host: str = 'http://www.ensembl.org'):
        """
        Inizializza la connessione a BioMart.

        Args:
            host: URL del server Ensembl BioMart
        """
        logger.info(f"Connessione a BioMart: {host}")
        self.dataset = Dataset(
            name='hsapiens_gene_ensembl',
            host=host
        )
        logger.info("Connessione stabilita")

    def query_chromosome(
        self, 
        chromosome: str,
        attributes: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Effettua query per un singolo cromosoma.

        Args:
            chromosome: Nome del cromosoma (1-22, X, Y, MT)
            attributes: Lista attributi da recuperare (default: attributi base)

        Returns:
            DataFrame con mappings per il cromosoma specificato
        """
        if attributes is None:
            attributes = self._get_default_attributes()

        logger.info(f"Query cromosoma {chromosome}...")

        try:
            df = self.dataset.query(
                attributes=attributes,
                filters={'chromosome_name': [chromosome]}
            )

            logger.info(f"Cromosoma {chromosome}: {len(df)} righe recuperate")
            return df

        except Exception as e:
            logger.error(f"Errore query cromosoma {chromosome}: {e}")
            raise

    def _get_default_attributes(self) -> List[str]:
        """
        Restituisce la lista di attributi di default per il mapping ENSG->ENST.

        Returns:
            Lista di attributi BioMart
        """
        return [
            # Gene identifiers
            'ensembl_gene_id',
            'ensembl_gene_id_version',
            'external_gene_name',
            'gene_biotype',

            # Transcript identifiers
            'ensembl_transcript_id',
            'ensembl_transcript_id_version',
            'external_transcript_name',
            'transcript_biotype',
            'transcript_is_canonical',
            'transcript_tsl',  # Transcript Support Level
            'transcript_appris',  # Principal isoform annotation
            'transcript_mane_select',  # MANE Select

            # Genomic coordinates
            'chromosome_name',
            'start_position',
            'end_position',
            'strand',
            'transcript_start',
            'transcript_end',
            'transcript_length',

            # Additional info
            'description',
        ]

    def create_full_mapping(
        self,
        chromosomes: Optional[List[str]] = None,
        output_path: Optional[str] = None,
        save_temp: bool = True,
        temp_dir: str = './temp'
    ) -> pd.DataFrame:
        """
        Crea il mapping completo per tutti i cromosomi.

        Args:
            chromosomes: Lista cromosomi da processare (default: tutti)
            output_path: Path file output (default: ensg_enst_mapping.tsv)
            save_temp: Se True, salva file temporanei per cromosoma
            temp_dir: Directory per file temporanei

        Returns:
            DataFrame con mapping completo
        """
        if chromosomes is None:
            chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']

        if output_path is None:
            output_path = 'ensg_enst_mapping.tsv'

        # Crea directory temporanea se necessaria
        if save_temp:
            Path(temp_dir).mkdir(parents=True, exist_ok=True)

        master_df = pd.DataFrame()

        for chrom in chromosomes:
            logger.info(f"\n{'='*60}")
            logger.info(f"Processamento cromosoma: {chrom}")
            logger.info(f"{'='*60}")

            try:
                # Query cromosoma
                df_chrom = self.query_chromosome(chrom)

                # Salva file temporaneo
                if save_temp:
                    temp_file = Path(temp_dir) / f"chr{chrom}_mapping.tsv"
                    df_chrom.to_csv(temp_file, sep='\t', index=False)
                    logger.info(f"File temporaneo salvato: {temp_file}")

                # Aggiungi al master dataframe
                master_df = pd.concat([master_df, df_chrom], ignore_index=True)

                # Pausa per non sovraccaricare il server
                time.sleep(1)

            except Exception as e:
                logger.error(f"Errore processamento cromosoma {chrom}: {e}")
                continue

        # Rinomina colonne per maggiore chiarezza
        master_df = self._rename_columns(master_df)

        # Salva risultato finale
        logger.info(f"\n{'='*60}")
        logger.info(f"Salvataggio dataset finale: {output_path}")
        master_df.to_csv(output_path, sep='\t', index=False)

        # Statistiche finali
        self._print_statistics(master_df)

        return master_df

    def _rename_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Rinomina le colonne con nomi più chiari.

        Args:
            df: DataFrame originale

        Returns:
            DataFrame con colonne rinominate
        """
        column_mapping = {
            'Gene stable ID': 'gene_id',
            'Gene stable ID version': 'gene_id_version',
            'Gene name': 'gene_name',
            'Gene type': 'gene_biotype',
            'Transcript stable ID': 'transcript_id',
            'Transcript stable ID version': 'transcript_id_version',
            'Transcript name': 'transcript_name',
            'Transcript type': 'transcript_biotype',
            'Transcript is canonical': 'is_canonical',
            'Transcript support level (TSL)': 'tsl',
            'APPRIS annotation': 'appris',
            'MANE Select transcript': 'mane_select',
            'Chromosome/scaffold name': 'chromosome',
            'Gene start (bp)': 'gene_start',
            'Gene end (bp)': 'gene_end',
            'Strand': 'strand',
            'Transcript start (bp)': 'transcript_start',
            'Transcript end (bp)': 'transcript_end',
            'Transcript length (including UTRs and CDS)': 'transcript_length',
            'Gene description': 'description',
        }

        return df.rename(columns=column_mapping)

    def _print_statistics(self, df: pd.DataFrame):
        """
        Stampa statistiche sul dataset creato.

        Args:
            df: DataFrame con mappings
        """
        logger.info(f"{'='*60}")
        logger.info("STATISTICHE DATASET")
        logger.info(f"{'='*60}")
        logger.info(f"Righe totali: {len(df):,}")
        logger.info(f"Geni unici (ENSG): {df['gene_id'].nunique():,}")
        logger.info(f"Trascritti unici (ENST): {df['transcript_id'].nunique():,}")
        logger.info(f"Media trascritti per gene: {len(df) / df['gene_id'].nunique():.2f}")

        logger.info(f"\nDistribuzione gene_biotype:")
        for biotype, count in df['gene_biotype'].value_counts().head(10).items():
            logger.info(f"  {biotype}: {count:,}")

        if 'is_canonical' in df.columns:
            canonical_count = df['is_canonical'].sum()
            logger.info(f"\nTrascritti canonici: {canonical_count:,}")

        if 'tsl' in df.columns:
            logger.info(f"\nDistribuzione TSL:")
            for tsl, count in df['tsl'].value_counts().items():
                logger.info(f"  TSL {tsl}: {count:,}")


def main():
    """Main entry point: query BioMart for the full ENSG -> ENST mapping."""
    from pathlib import Path

    # Anchor outputs to <repo_root>/data/mappings/ regardless of the cwd.
    REPO_ROOT = Path(__file__).resolve().parents[2]
    OUTPUT_FILE = str(REPO_ROOT / "data" / "mappings" / "biomart_ensg_enst_mapping.tsv")
    TEMP_DIR = str(REPO_ROOT / "data" / "mappings" / "temp" / "ensg_enst_chr")

    mapper = ENSGtoENSTMapper()

    # Esegui mapping completo
    # Per test rapido, usa subset di cromosomi: chromosomes=['21', '22']
    df_mapping = mapper.create_full_mapping(
        chromosomes=None,  # None = tutti i cromosomi
        output_path=OUTPUT_FILE,
        save_temp=True,
        temp_dir=TEMP_DIR
    )

    logger.info(f"\n{'='*60}")
    logger.info("COMPLETATO!")
    logger.info(f"File output: {OUTPUT_FILE}")
    logger.info(f"{'='*60}")

    return df_mapping


if __name__ == "__main__":
    main()
