# Re-importing the necessary preprocessing and model training tools after disconnection
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score

# Reload the dataset
file_path = '/Users/sq566/Desktop/sophia/ML_dataset.csv'
df = pd.read_csv(file_path)

# Encode categorical variables
label_encoder = LabelEncoder()

# Encoding 'Gene', 'Mutation', 'Effect'
df['Gene_encoded'] = label_encoder.fit_transform(df['Gene'])
df['Mutation_encoded'] = label_encoder.fit_transform(df['Mutation'])
df['Effect_encoded'] = label_encoder.fit_transform(df['Effect'])

# Selecting features and target
X = df[['Gene_encoded', 'Mutation_encoded', 'Location', 'Effect_encoded']]
y = df['ClinicalOutcome']

# Splitting the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Define the parameter grid for Gradient Boosting
param_grid = {
    'n_estimators': [100, 200, 300],
    'learning_rate': [0.01, 0.1, 0.2],
    'max_depth': [3, 5, 7],
    'subsample': [0.8, 1.0],
    'min_samples_split': [2, 5, 10]
}

# Initialize the Gradient Boosting classifier
gb_model = GradientBoostingClassifier(random_state=42)

# Perform GridSearchCV for hyperparameter tuning
grid_search = GridSearchCV(estimator=gb_model, param_grid=param_grid, cv=5, n_jobs=-1, verbose=1)
grid_search.fit(X_train, y_train)

# Get the best parameters from the grid search
best_params = grid_search.best_params_

# Train the Gradient Boosting model with the best parameters
best_gb_model = GradientBoostingClassifier(**best_params, random_state=42)
best_gb_model.fit(X_train, y_train)

# Make predictions on the test set
y_pred_best_gb = best_gb_model.predict(X_test)

# Evaluate the model's accuracy
accuracy_best_gb = accuracy_score(y_test, y_pred_best_gb)

accuracy_best_gb, best_params
