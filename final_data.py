import pandas as pd
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split

# Load datasets
nfhs_df = pd.read_csv('NFHS-5-States.csv', low_memory=False)
data_df = pd.read_csv('data.csv')
who_excel = pd.ExcelFile('who-data.xlsx')

# Load first sheet from WHO Excel file
who_df = who_excel.parse(who_excel.sheet_names[0])

# Filter WHO data for India only
if 'Country' in who_df.columns:
    who_df_india = who_df[who_df['Country'].str.lower() == 'india']
else:
    who_df_india = who_df.copy()

# Define relevant columns to look for
key_columns_all = ['Year', 'State', 'District', 'Gender', 'Age_Group', 'Urban_Rural']
key_columns_nfhs = [col for col in key_columns_all if col in nfhs_df.columns]
key_columns_data = [col for col in key_columns_all if col in data_df.columns]
key_columns_who = [col for col in key_columns_all if col in who_df_india.columns]

# Identify obesity-related columns in datasets
obesity_keywords = ['obesity', 'bmi', 'overweight', 'prevalence']
obesity_cols_nfhs = [col for col in nfhs_df.columns if any(kw in col.lower() for kw in obesity_keywords)]
obesity_cols_data = [col for col in data_df.columns if any(kw in col.lower() for kw in obesity_keywords)]
obesity_cols_who = [col for col in who_df_india.columns if any(kw in col.lower() for kw in obesity_keywords)]

# Filter datasets
nfhs_filtered = nfhs_df[key_columns_nfhs + obesity_cols_nfhs]
data_filtered = data_df[key_columns_data + obesity_cols_data]
who_filtered = who_df_india[key_columns_who + obesity_cols_who]

# Merge NFHS and data.csv on intersecting key columns
common_keys_nfhs_data = list(set(key_columns_nfhs).intersection(set(key_columns_data)))
if common_keys_nfhs_data:
    merged1 = pd.merge(nfhs_filtered, data_filtered, on=common_keys_nfhs_data, how='outer')
else:
    merged1 = pd.concat([nfhs_filtered, data_filtered], axis=1)

# Merge with WHO data on intersecting keys treating population independently
common_keys_with_who = list(set(common_keys_nfhs_data).intersection(set(key_columns_who)))
if common_keys_with_who:
    merged = pd.merge(merged1, who_filtered.drop(columns=['Population'], errors='ignore'), on=common_keys_with_who, how='outer')
else:
    merged = pd.concat([merged1, who_filtered], axis=1)

# Handle missing numeric data by median imputation
numeric_cols = merged.select_dtypes(include=['float64', 'int64']).columns.tolist()
for col in numeric_cols:
    merged[col] = merged[col].fillna(merged[col].median())

# Encode categorical variables
categorical_cols = merged.select_dtypes(include=['object']).columns.tolist()
le_dict = {}
for col in categorical_cols:
    le = LabelEncoder()
    merged[col] = le.fit_transform(merged[col].astype(str))
    le_dict[col] = le

# Normalize numeric columns if present
if len(numeric_cols) > 0:
    scaler = StandardScaler()
    merged[numeric_cols] = scaler.fit_transform(merged[numeric_cols])

# Identify a label column for obesity, pick first matching column if any
possible_labels = [col for col in merged.columns if any(kw in col.lower() for kw in obesity_keywords)]
label_col = possible_labels[0] if possible_labels else None

if label_col:
    X = merged.drop(columns=[label_col])
    y = merged[label_col]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
else:
    X_train = X_test = y_train = y_test = None

# Save final ML ready dataset and splits
merged.to_csv('final_ml_ready_dataset.csv', index=False)
if X_train is not None:
    X_train.to_csv('X_train.csv', index=False)
    X_test.to_csv('X_test.csv', index=False)
    y_train.to_csv('y_train.csv', index=False)
    y_test.to_csv('y_test.csv', index=False)

# Summary report
summary = {
    'final_dataset_shape': merged.shape,
    'label_column': label_col,
    'train_set_shape': X_train.shape if X_train is not None else None,
    'test_set_shape': X_test.shape if X_test is not None else None
}

summary
