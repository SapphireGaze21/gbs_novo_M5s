import pandas as pd
import numpy as np

# File names
nfhf = 'NFHS-5-States.csv'
data_csv = 'data.csv'
who_xlsx = 'who-data.xlsx'
pubmed_csv = 'cleaned_pubmed_data.csv'

# Step 1: Load datasets

df_nfhs = pd.read_csv(nfhf, low_memory=False)
df_data = pd.read_csv(data_csv, low_memory=False)
df_pubmed = pd.read_csv(pubmed_csv, low_memory=False)
df_who = pd.read_excel(who_xlsx, sheet_name=None)

# Step 2: Extract all meaningful numeric columns (drop ID-like, or nearly all NaN)
def is_sufficient_numeric(col):
    if not np.issubdtype(col.dtype, np.number):
        return False
    # At least 20% non-NA, not all same value
    return col.count() > 0.2 * len(col) and col.nunique() > 1

def extract_numeric(df):
    numerics = df.loc[:, df.apply(is_sufficient_numeric, axis=0)]
    return numerics

nfhs_numeric = extract_numeric(df_nfhs)
data_numeric = extract_numeric(df_data)
pubmed_numeric = extract_numeric(df_pubmed)

# For WHO XLSX: scan all sheets and extract numeric columns
df_who_numeric_list = []
for name, sheet in df_who.items():
    numerics = extract_numeric(sheet)
    if not numerics.empty:
        numerics['sheet'] = name
        df_who_numeric_list.append(numerics)
if df_who_numeric_list:
    who_numeric = pd.concat(df_who_numeric_list, ignore_index=True)
else:
    who_numeric = pd.DataFrame()

# Step 3: Merge all extracted numeric columns (align on available shared keys if possible)
# For illustration, simply concatenate along columns (as features may not align)
combined_numeric = pd.concat([
    nfhs_numeric.add_prefix('NFHS5_'),
    data_numeric.add_prefix('DATA_'),
    pubmed_numeric.add_prefix('PUBMED_'),
    who_numeric.add_prefix('WHO_')
], axis=1)

# Step 4: Save an ML-suitable CSV
combined_numeric.to_csv('ml_ready_numeric_data.csv', index=False)

# Save column descriptions for user review
with open('ml_ready_numeric_columns.txt','w',encoding='utf-8') as f:
    for col in combined_numeric.columns:
        f.write(col+'\n')

# Print shape for confirmation
combined_numeric.shape