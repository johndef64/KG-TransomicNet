#%%
"""
This script creates a sample of 1000 rows from PKT dataframes.
Processes two compressed files:
- PheKnowLator NT triples file (tar.gz)
- PheKnowLator CSV file (gz)
Outputs simple CSV samples.
"""
from  pkt_utils import *
set_working_directory()

file_root = "PheKnowLator_v3.0.2_full_instance_inverseRelations_OWLNETS_INSTANCE_purified"
file_root = "PKT"

kg_files = {
# other files
1:f"{file_root}_NetworkxMultiDiGraph.gpickle.tar.gz",
2:f"{file_root}_decoding_dict.pkl.tar.gz",

# ntriples
3:f"{file_root}.nt.tar.gz",

# tables
5:f"{file_root}_NodeLabels.txt.tar.gz",
6:f"{file_root}_Triples_Identifiers.txt.tar.gz",
7:f"{file_root}_Triples_Integers.txt.tar.gz",

}

import pandas as pd
import gzip
import tarfile
import os


def sample_nt_file(input_file, output_file, sample_size=1000):
    """
    Sample NT triples file and save as CSV.
    NT format: subject predicate object .
    """
    print(f"Processing NT file: {input_file}")
    
    triples = []
    with tarfile.open(input_file, 'r:gz') as tar:
        # Get the first member (should be the .nt file)
        members = tar.getmembers()
        if not members:
            print("No files found in tar.gz")
            return
        
        nt_member = members[0]
        print(f"Extracting: {nt_member.name}")
        
        # Read the NT file
        f = tar.extractfile(nt_member)
        for i, line in enumerate(f):
            if i >= sample_size:
                break
            
            line = line.decode('utf-8').strip()
            if not line or line.startswith('#'):
                continue
            
            # Parse NT triple: subject predicate object .
            parts = line.rstrip(' .').split(None, 2)
            if len(parts) == 3:
                # Remove < and > from URIs
                subject = parts[0].strip('<>')
                predicate = parts[1].strip('<>')
                obj = parts[2].strip('<>')
                
                triples.append({
                    'subject': subject,
                    'predicate': predicate,
                    'object': obj
                })
    
    # Create DataFrame and save
    df = pd.DataFrame(triples)
    df.to_csv(output_file, index=False)
    print(f"Saved {len(df)} rows to {output_file}")
    print(f"Sample preview:\n{df.head()}\n")

def sample_csv_file(input_file, output_file, sample_size=1000):
    """
    Sample CSV file and save.
    """
    print(f"Processing CSV file: {input_file}")
    
    # Read compressed CSV
    df = pd.read_csv(input_file, nrows=sample_size, sep='\t')
    
    # Save sample
    df.to_csv(output_file, index=False)
    print(f"Saved {len(df)} rows to {output_file}")
    print(f"Sample preview:\n{df.head()}\n")

if __name__ == "__main__":
    print("Starting PKT data sampling...\n")

    # Define file paths
    for n in [5,6,7]:  # CSV files
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        NT_FILE = os.path.join(SCRIPT_DIR, kg_files[3])
        CSV_FILE = os.path.join(SCRIPT_DIR, kg_files[n])

        # Output files
        NT_SAMPLE_OUTPUT = os.path.join(f"{NT_FILE}_SAMPLE_1000.csv")
        CSV_SAMPLE_OUTPUT = os.path.join(f"{CSV_FILE}_SAMPLE_1000.csv")


        
        # Sample NT file
        if not os.path.exists(NT_SAMPLE_OUTPUT):
            if os.path.exists(NT_FILE):
                sample_nt_file(NT_FILE, NT_SAMPLE_OUTPUT, sample_size=1000)
            else:
                print(f"NT file not found: {NT_FILE}")
        else:
            print(f"NT sample file already exists: {NT_SAMPLE_OUTPUT}")
            
        # Sample CSV file
        if not os.path.exists(CSV_SAMPLE_OUTPUT):
            if os.path.exists(CSV_FILE):
                sample_csv_file(CSV_FILE, CSV_SAMPLE_OUTPUT, sample_size=1000)
            else:
                print(f"CSV file not found: {CSV_FILE}")
        else:
            print(f"CSV sample file already exists: {CSV_SAMPLE_OUTPUT}")


#%%
