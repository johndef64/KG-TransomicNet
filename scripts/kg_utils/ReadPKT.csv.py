#%%
import pandas as pd

pkl_file = r'..\temp_dir\PheKnowLator_v2.0.0_REDUX_NO-GO.csv.gz'
pkf_df = pd.read_csv(pkl_file)
pkf_df.head()
#%%
print(pkf_df.head().to_csv())
# pkf_df data structure (Tabular URI):
#  ,subject,predicate,object
# 0,http://purl.obolibrary.org/obo/PR_Q7Z4G1,http://purl.obolibrary.org/obo/RO_0002436,http://purl.obolibrary.org/obo/PR_Q96L50
# 1,http://purl.obolibrary.org/obo/CHEBI_7553,http://purl.obolibrary.org/obo/RO_0002434,https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000586592
# 2,http://purl.obolibrary.org/obo/CHEBI_46024,http://purl.obolibrary.org/obo/RO_0002434,https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000468105
# 3,http://purl.obolibrary.org/obo/PR_Q8N2Z9,http://purl.obolibrary.org/obo/RO_0002436,http://purl.obolibrary.org/obo/PR_Q9NVI1-3
# 4,http://www.ncbi.nlm.nih.gov/gene/145773,http://purl.obolibrary.org/obo/RO_0002511,https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000557914

# Extract subset: only lines that contain "rdf" in any column
redf_pkt = pkf_df[pkf_df.apply(lambda row: row.astype(str).str.contains('www.w3.org').any(), axis=1)]
print("Rows containing 'rdf' in any column:")
print(redf_pkt)

# Compare number of rows in native dataset vs subset
print(f"Total rows in dataset: {len(pkf_df)}")
print(f"Rows with 'rdf': {len(redf_pkt)}")

#%%

#%%
