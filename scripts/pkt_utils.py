#%%

# set working dir
import os

def set_working_directory():
    print(f"Current working directory: {os.getcwd()}")
    if os.getcwd().endswith("scripts"):
        print("Changing working directory to ..\\data\\pkt\\builds\\v3.0.2")
        os.chdir("..\\data\\pkt\\builds\\v3.0.2")
        print(f"Current working directory: {os.getcwd()}")

def print_file_contents(file_path=os.getcwd(), num_lines=10):
    """Print the files contained in the working directory sorted for file size and extension.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(os.path.join(file_path, f))]
    file_info = []
    for f in files:
        full_path = os.path.join(file_path, f)
        size = os.path.getsize(full_path)
        ext = os.path.splitext(f)[1]
        file_info.append((f, size, ext))
    
    # Sort by size (descending) and then by extension
    file_info.sort(key=lambda x: (-x[1], x[2]))
    
    print(f"Files in directory: {file_path}\n")
    for f, size, ext in file_info:
        print(f"File: {f} -  Size: {size / (1024*1024):.2f} MB -  Extension: {ext}\n")
        # Print first num_lines of the file
        # try:

        #     with open(os.path.join(file_path, f), 'r', encoding='utf-8') as file:
        #         print(f"  First {num_lines} lines:")
        #         for i in range(num_lines):
        #             line = file.readline()
        #             if not line:
        #                 break
        #             print(f"    {line.strip()}")
        # except Exception as e:
        #     print(f"  Could not read file: {e}")
        print("-" * 40)


import pandas as pd
import tarfile

def read_tar_rdf_csv(file_path):
    print(f"Reading file ID  {file_path}")    
    with tarfile.open(file_path, 'r:gz') as tar:
        members = tar.getmembers()
        extracted_file = tar.extractfile(members[0])
        df = pd.read_csv(extracted_file, sep=' ', header=None, names=['subject', 'predicate', 'object', '.'])
        return df[['subject', 'predicate', 'object']]



if __name__ == "__main__":
    set_working_directory()
    # Example usage of print_file_contents
    # print_file_contents("example_file.txt", 5)
    print_file_contents()