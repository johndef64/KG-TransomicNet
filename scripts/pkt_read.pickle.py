#%%
# read  pickle the file PheKnowLator_v3.0.2_full_instance_inverseRelations_OWLNETS_INSTANCE_purified_decoding_dict.pkl.tar

import pickle
import tarfile
import io
import os
from  pkt_utils import *
set_working_directory()

# Open the tar.gz file and extract the pickle
tar_file_path = "PheKnowLator_v3.0.2_full_instance_inverseRelations_OWLNETS_INSTANCE_purified_decoding_dict.pkl.tar.gz"

with tarfile.open(tar_file_path, 'r:gz') as tar:
    # List all files in the archive
    print("Files in archive:")
    for member in tar.getmembers():
        print(f"  - {member.name}")
    
    # Extract the first file (should be the .pkl file)
    tar.extractall()
    
    # Get the name of the extracted file
    pkl_file = tar.getmembers()[0].name
    print(f"\nReading pickle file: {pkl_file}")

# Load the pickle file
with open(pkl_file, 'rb') as f:
    decoding_dict = pickle.load(f)

print(f"\nLoaded decoding dictionary with {len(decoding_dict)} entries")
print(f"Type: {type(decoding_dict)}")

# # Show a sample of the data
# print("\nFirst 10 entries:")
# for i, (key, value) in enumerate(list(decoding_dict.items())[:10]):
#     print(f"  {key}: {value}")

# # Show basic statistics
# if isinstance(decoding_dict, dict):
#     print(f"\nKeys type: {type(list(decoding_dict.keys())[0]) if decoding_dict else 'N/A'}")
#     print(f"Values type: {type(list(decoding_dict.values())[0]) if decoding_dict else 'N/A'}")

# %%
decoding_dict.keys()
# pretti print(list(decoding_dict.keys())[:10])
# %%
import pandas as pd
list(decoding_dict.values())[6]
# pd.Series(list(decoding_dict.values())).value_counts().head(10)
# %%
# save decodict as json with indentation
# Need to convert sets to lists for JSON serialization
import json

def convert_sets_to_lists(obj):
    """Recursively convert sets to lists for JSON serialization"""
    if isinstance(obj, set):
        return list(obj)
    elif isinstance(obj, dict):
        return {key: convert_sets_to_lists(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_sets_to_lists(item) for item in obj]
    else:
        return obj

# Convert sets to lists
decoding_dict_serializable = convert_sets_to_lists(decoding_dict)

with open("decoding_dict_v3.0.2.json", "w") as f:
    json.dump(decoding_dict_serializable, f, indent=4)
    
print("Successfully saved decoding_dict to JSON!")