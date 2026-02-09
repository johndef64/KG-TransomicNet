#%%
import json
import os
from pathlib import Path
from  pkt_utils import *
set_working_directory()

# Define the file path
json_file = "PheKnowLator_v3.0.2_full_instance_inverseRelations_OWLNETS_Triples_Integer_Identifier_Map.json"
file_path = Path(__file__).parent / json_file

print(f"Loading JSON file: {file_path}")
print(f"File exists: {file_path.exists()}")
print(f"File size: {os.path.getsize(file_path) / (1024*1024):.2f} MB\n")

# Load the JSON file
try:
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    print("JSON file loaded successfully!")
    print(f"Type: {type(data)}")
    
    # Display basic information about the data
    if isinstance(data, dict):
        print(f"Number of keys: {len(data)}")
        print(f"\nFirst 10 keys:")
        for i, key in enumerate(list(data.keys())[:10]):
            print(f"  {i+1}. {key}: {data[key]}")
        
        if len(data) > 10:
            print(f"\n... and {len(data) - 10} more entries")
    
    elif isinstance(data, list):
        print(f"Number of items: {len(data)}")
        print(f"\nFirst 5 items:")
        for i, item in enumerate(data[:5]):
            print(f"  {i+1}. {item}")
        
        if len(data) > 5:
            print(f"\n... and {len(data) - 5} more items")
    
    # You can now work with the 'data' variable
    # Example: print specific entries, analyze structure, etc.
    
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
except json.JSONDecodeError as e:
    print(f"Error decoding JSON: {e}")
except Exception as e:
    print(f"Error loading file: {e}")


#%%
# load json in pandas
import pandas as pd

# Since the JSON is a dictionary mapping, we need to convert it properly
# Option 1: Create a DataFrame from the dictionary with orient='index'
df = pd.DataFrame.from_dict(data, orient='index', columns=['value'])
print(f"DataFrame shape: {df.shape}")
print(f"\nFirst 10 rows:")
print(df.head(10))
df
# Option 2: If you want key-value pairs as columns
# df = pd.DataFrame(list(data.items()), columns=['key', 'value'])

#%%


# Load NodeLabels_with_metadata_v3.0.2.csv and in clumn bioentity_type, replace rna with sequence
NodeLabels_meta = pd.read_csv("PKT_NodeLabels_with_metadata_v3.0.2.csv")
NodeLabels_meta['bioentity_type'] = NodeLabels_meta['bioentity_type'].replace('rna', 'sequence')
#%%
NodeLabels_meta.head()