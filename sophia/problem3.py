import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

# Load the CSV file into a DataFrame
df = pd.read_csv('/Users/sq566/Desktop/sophia/ML_dataset.csv')

# Step 1: Encode categorical features
le_gene = LabelEncoder()
le_mutation = LabelEncoder()
le_effect = LabelEncoder()

df['Gene_encoded'] = le_gene.fit_transform(df['Gene'])
df['Mutation_encoded'] = le_mutation.fit_transform(df['Mutation'])
df['Effect_encoded'] = le_effect.fit_transform(df['Effect'])

# Drop original categorical columns
df_processed = df.drop(['Gene', 'Mutation', 'Effect'], axis=1)

# Step 2: Define features (X) and target (y)
X = df_processed.drop(['ClinicalOutcome'], axis=1)
y = df_processed['ClinicalOutcome']

# Step 3: Train-test split (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Step 4: Train a Random Forest model
rf_model = RandomForestClassifier(random_state=42)
rf_model.fit(X_train, y_train)

# Step 5: Predictions and evaluation
y_pred = rf_model.predict(X_test)
report = classification_report(y_test, y_pred)

# Output classification report
print(report)


