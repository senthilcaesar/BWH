import pandas as pd

# Load the CSV file into a pandas dataframe

file_path = '/Users/sq566/Desktop/sophia/mutation_data.csv'
mutation_data = pd.read_csv(file_path)

# Display the first few rows of the dataframe to understand its structure
mutation_data.head()

