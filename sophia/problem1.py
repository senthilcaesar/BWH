import pandas as pd

# Load Dataset_A.csv into a DataFrame
df_A = pd.read_csv('/Users/sq566/Desktop/sophia/Dataset_A.csv')

# Load Dataset_B.csv into a DataFrame
df_B = pd.read_csv('/Users/sq566/Desktop/sophia/Dataset_B.csv')

# Display the first few rows of both DataFrames
df_A.head(), df_B.head()

# Renaming the columns to align similar data for merging
df_B_renamed = df_B.rename(columns={ "PatientID": "SampleID"})  # Assuming "PatientID" and "SampleID" are the same

# Replacing 'P' with 'S' in SampleID column of Dataset_B to match the format of SampleID in Dataset_A
df_B_renamed['SampleID'] = df_B_renamed['SampleID'].str.replace('P', 'S')

# Merging the datasets on the common column "SampleID"
merged_df = pd.merge(df_A, df_B_renamed, on="SampleID", how="outer")

# Checking for duplicate rows and removing them
merged_df_no_duplicates = merged_df.drop_duplicates()

# Identifying missing data
missing_data = merged_df_no_duplicates.isnull().sum()

print(missing_data)

# Get all rows that contain NaN values
rows_with_nan = merged_df_no_duplicates[merged_df_no_duplicates.isnull().any(axis=1)]

# Print the rows with NaN values
print(rows_with_nan)

merged_df_filled = merged_df_no_duplicates.fillna("Unknown")

# Dropping rows with 'Unknown' values
df_cleaned = merged_df_filled.replace("Unknown", pd.NA).dropna()

# Checking for anomalies in columns
categorical_columns = ['Mutation', 'Alteration', 'Chromosome', 'Impact' ,'Effect']

# Getting the unique values in each column to inspect for anomalies
anomalies_categorical = {col: sorted(df_cleaned[col].unique()) for col in categorical_columns}

# Checking the range of the 'Location' column for outliers (assuming valid locations are non-negative)
location_outliers = df_cleaned[df_cleaned['Location'] < 0]

# Displaying results: unique values in categorical columns and any outliers in 'Location'
anomalies_categorical, location_outliers

merged_data = merged_df_filled
merged_data.head()