import requests
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path

path = "../data/pkt/builds/v3.0.2/PKT_NodeLabels_with_metadata_v3.0.2.csv"
pkt_df = pd.read_csv(path)
pkt_df.head()

# get all "entity" for class_code == "dbSNP"
rsids = pkt_df[pkt_df['class_code'] == 'dbSNP']['entity'].unique().tolist()
print(f"Total rsids in PKT: {len(rsids)}")

def fetch_rsid_data(rsid, session):
    """Fetch data for a single rsID using a shared session."""
    url = f"https://clinicaltables.nlm.nih.gov/api/snps/v3/search?terms={rsid}&q=rsNum:{rsid}&ef=38.chr,38.pos,38.alleles,38.gene"
    
    try:
        response = session.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        # Convert to dataframe
        if data[0] > 0:  # Check if results exist
            df_rs = pd.DataFrame(data[2])
            df_rs["snp_id"] = rsid
            # Move snp_id to first column
            cols = ['snp_id'] + [col for col in df_rs.columns if col != 'snp_id']
            df_rs = df_rs[cols]
            return df_rs
        else:
            # Return empty dataframe with snp_id if no results
            return pd.DataFrame({'snp_id': [rsid]})
    except Exception as e:
        print(f"Error fetching {rsid}: {e}")
        return pd.DataFrame({'snp_id': [rsid], 'error': [str(e)]})

def load_checkpoint(checkpoint_file):
    """Load checkpoint file and return processed rsIDs."""
    if os.path.exists(checkpoint_file):
        print(f"Loading checkpoint from {checkpoint_file}")
        df = pd.read_csv(checkpoint_file, sep='\t')
        processed_rsids = set(df['snp_id'].unique())
        print(f"Found {len(processed_rsids)} already processed rsIDs")
        return df, processed_rsids
    return pd.DataFrame(), set()

def save_checkpoint(df, checkpoint_file):
    """Save checkpoint to TSV file."""
    df.to_csv(checkpoint_file, sep='\t', index=False)
    print(f"Checkpoint saved: {len(df)} total records")

def fetch_multiple_rsids_with_checkpoint(rsids, checkpoint_file='rsid_checkpoint.tsv', 
                                         checkpoint_interval=500, max_workers=10):
    """
    Fetch data for multiple rsIDs using multithreading with checkpoint support.
    
    Parameters:
    -----------
    rsids : list
        List of rsIDs to fetch
    checkpoint_file : str
        Path to checkpoint TSV file
    checkpoint_interval : int
        Save checkpoint every N results
    max_workers : int
        Number of concurrent threads
    """
    
    # Load existing checkpoint if available
    master_rs_df, processed_rsids = load_checkpoint(checkpoint_file)
    
    # Filter out already processed rsIDs
    remaining_rsids = [rsid for rsid in rsids if rsid not in processed_rsids]
    
    if not remaining_rsids:
        print("All rsIDs already processed!")
        return master_rs_df
    
    print(f"Processing {len(remaining_rsids)} remaining rsIDs out of {len(rsids)} total")
    
    # Create a session with connection pooling
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(
        pool_connections=max_workers,
        pool_maxsize=max_workers
    )
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    
    results = []
    completed_count = 0
    
    # Use ThreadPoolExecutor for concurrent requests
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_rsid = {
            executor.submit(fetch_rsid_data, rsid, session): rsid 
            for rsid in remaining_rsids
        }
        
        # Process completed tasks with progress bar
        for future in tqdm(as_completed(future_to_rsid), total=len(remaining_rsids)):
            rsid = future_to_rsid[future]
            try:
                df_result = future.result()
                results.append(df_result)
                completed_count += 1
                
                # Save checkpoint every checkpoint_interval results
                if completed_count % checkpoint_interval == 0:
                    temp_df = pd.concat(results, ignore_index=True)
                    current_df = pd.concat([master_rs_df, temp_df], ignore_index=True)
                    save_checkpoint(current_df, checkpoint_file)
                    
            except Exception as e:
                print(f"Exception for {rsid}: {e}")
                results.append(pd.DataFrame({'snp_id': [rsid], 'error': [str(e)]}))
    
    # Final concatenation and save
    if results:
        new_results_df = pd.concat(results, ignore_index=True)
        master_rs_df = pd.concat([master_rs_df, new_results_df], ignore_index=True)
    
    # Save final checkpoint
    save_checkpoint(master_rs_df, checkpoint_file)
    
    return master_rs_df

# Usage
# rsids = ['rs121913529', 'rs121913530', 'rs121913531', ...]
master_rs_df = fetch_multiple_rsids_with_checkpoint(
    rsids, 
    checkpoint_file='../data/omics/maps/rsid_checkpoint.tsv',
    checkpoint_interval=1000,
    max_workers=10
)
master_rs_df
