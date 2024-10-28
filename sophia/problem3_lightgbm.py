import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import lightgbm as lgb
from lightgbm import early_stopping
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import accuracy_score, precision_score, roc_auc_score, recall_score, f1_score, classification_report


# Step 1: Data Preprocessing
# Load the CSV file into a DataFrame
df = pd.read_csv('/Users/sq566/Desktop/sophia/ML_dataset.csv')

# Encode the categorical features in a separate variable
encoded_features = df[['Gene', 'Mutation', 'Effect']].copy()

# Apply label encoding to each categorical feature
for col in encoded_features.columns:
    le = LabelEncoder()
    encoded_features[col] = le.fit_transform(df[col])

# Combine the encoded features with the rest of the data (e.g., 'Location')
X = pd.concat([encoded_features, df[['Location']]], axis=1)


y = df['ClinicalOutcome']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# Initialize LightGBM model
lgb_train = lgb.Dataset(X_train, y_train)
lgb_test = lgb.Dataset(X_test, y_test, reference=lgb_train)

# Set parameters for LightGBM
params = {
    'objective': 'binary',
    'metric': 'auc',
    'boosting_type': 'gbdt',
    'num_leaves': 31,
    'learning_rate': 0.05,
    'verbosity': -1
}

# Train the model
lgb_model = lgb.train(params, lgb_train, valid_sets=lgb_test, num_boost_round=100, 
                      callbacks=[early_stopping(stopping_rounds=10)])


# Make predictions
y_pred_proba = lgb_model.predict(X_test, num_iteration=lgb_model.best_iteration)
y_pred = (y_pred_proba > 0.5).astype(int)  # converting probabilities to binary predictions

# Calculate AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

# Calculate other metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)

# Generate classification report
classification_rep = classification_report(y_test, y_pred)

# Display the metrics
print(f"Accuracy: {accuracy}")
print(f"Precision: {precision}")
print(f"Recall: {recall}")
print(f"F1-score: {f1}")
print(f"AUC score: {auc_score}")
print("Classification Report:\n", classification_rep)


# Calculate the ROC curve values
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)

# Plot the AUC curve
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='blue', label=f'AUC = {auc_score:.2f}')
plt.plot([0, 1], [0, 1], color='red', linestyle='--', label='Random Classifier')
plt.title('ROC Curve')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()

# Plot feature importance
lgb.plot_importance(lgb_model, max_num_features=10)
plt.show()
