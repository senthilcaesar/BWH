import pandas as pd

# Read the CSV file into a DataFrame
df = pd.read_csv('patient_medication_data.csv')

# Count the number of duplicate rows
duplicates_all = df.duplicated().sum()

# Print the results
print("The number of duplicate rows in the dataset is:", duplicates_all)

# Drop duplicate rows
df_cleaned = df.drop_duplicates()

# Count the number of duplicate rows in the cleaned dataset
duplicates_all_remaining = df_cleaned.duplicated().sum()

# Print the results
print("The number of duplicate rows in the dataset is:", duplicates_all_remaining)