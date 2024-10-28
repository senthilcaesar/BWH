import pandas as pd

# Load the provided CSV files
commercial_panel_path = 'commercial_panel.csv'
in_house_panel_path = 'in_house_panel.csv'

# Read the data into DataFrames
commercial_panel_df = pd.read_csv(commercial_panel_path)
in_house_panel_df = pd.read_csv(in_house_panel_path)

# Display the first few rows and the column names of each DataFrame for inspection
commercial_panel_df_info = commercial_panel_df.head(), commercial_panel_df.columns
in_house_panel_df_info = in_house_panel_df.head(), in_house_panel_df.columns

print(commercial_panel_df_info), print(in_house_panel_df_info)

merged_df = pd.merge(commercial_panel_df, in_house_panel_df, on=['Gene', 'MutationType'], how='outer', indicator=True)

# Split the merged data into shared mutations, unique to commercial panel, and unique to in-house panel
shared_mutations = merged_df[merged_df['_merge'] == 'both']
unique_commercial = merged_df[merged_df['_merge'] == 'left_only']
unique_in_house = merged_df[merged_df['_merge'] == 'right_only']

# Remove the merge indicator for clarity in results
shared_mutations_clean = shared_mutations.drop(columns=['_merge'])
unique_shared_mutations = shared_mutations_clean.drop_duplicates(subset=['Gene', 'MutationType', 'Annotation'])

commercial_clean = unique_commercial.drop(columns=['_merge'])
unique_commercial_clean = commercial_clean.drop_duplicates(subset=['Gene', 'MutationType', 'Annotation'])

in_house_clean = unique_in_house.drop(columns=['_merge', 'Annotation'])
unique_in_house_clean = in_house_clean.drop_duplicates(subset=['Gene', 'MutationType'])

discrepancies = shared_mutations_clean[shared_mutations_clean['Annotation'].isnull() | (shared_mutations_clean['Annotation'] != shared_mutations_clean['ResolvedAnnotation'])]

# Suggest resolution: Fill missing annotations from commercial panel
# Since in-house panel lacks annotations, we will use the commercial panel's annotation wherever available.
discrepancies['SuggestedResolution'] = discrepancies['ResolvedAnnotation'].fillna(discrepancies['Annotation'])


len(discrepancies)