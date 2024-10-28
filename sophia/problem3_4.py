import pandas as pd
from pytorch_tabular import TabularModel
from pytorch_tabular.config import DataConfig, ModelConfig
from pytorch_tabular.models import TabTransformerConfig
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

file_path = '/Users/sq566/Desktop/sophia/ML_dataset.csv'
df = pd.read_csv(file_path)

# Preparing the data (assuming df is the dataset)
label_encoder = LabelEncoder()
df['Gene_encoded'] = label_encoder.fit_transform(df['Gene'])
df['Mutation_encoded'] = label_encoder.fit_transform(df['Mutation'])
df['Effect_encoded'] = label_encoder.fit_transform(df['Effect'])

X = df[['Gene_encoded', 'Mutation_encoded', 'Location', 'Effect_encoded']]
y = df['ClinicalOutcome']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Prepare the dataset for TabTransformer
train_data = pd.concat([X_train, y_train], axis=1)
test_data = pd.concat([X_test, y_test], axis=1)

# Step 1: Define DataConfig
data_config = DataConfig(
    target=['ClinicalOutcome'],  # specify the target column name
    continuous_cols=['Location'],
    categorical_cols=['Gene_encoded', 'Mutation_encoded', 'Effect_encoded']
)

# Step 2: Define the TabTransformer ModelConfig
model_config = TabTransformerConfig(
    task="classification",
    #n_blocks=6,
    #n_heads=4,
    #attention_dropout=0.2,
    #mlp_hidden_dims=[128, 64]
)

# Step 3: Initialize the TabularModel with both configs
tabular_model = TabularModel(
    data_config=data_config,
    model_config=model_config
)

# Train the model
tabular_model.fit(train=train_data, test=test_data)

# Evaluate the model
result = tabular_model.evaluate(test_data)
print("Accuracy on the test set:", result)
