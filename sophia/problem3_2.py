# Preparing the data for modeling
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, accuracy_score

df = pd.read_csv('/Users/sq566/Desktop/sophia/ML_dataset.csv')

# Encode categorical variables
label_encoder = LabelEncoder()

# Encoding 'Gene', 'Mutation', 'Effect'
df['Gene_encoded'] = label_encoder.fit_transform(df['Gene'])
df['Mutation_encoded'] = label_encoder.fit_transform(df['Mutation'])
df['Effect_encoded'] = label_encoder.fit_transform(df['Effect'])

# Selecting features and target
X = df[['Gene_encoded', 'Mutation_encoded', 'Effect_encoded']]
y = df['ClinicalOutcome']

# Splitting the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Initializing the Random Forest classifier
rf_model = RandomForestClassifier(n_estimators=20)

# Training the model
rf_model.fit(X_train, y_train)

# Making predictions
y_pred = rf_model.predict(X_test)

# Evaluating the model
accuracy = accuracy_score(y_test, y_pred)
classification_rep = classification_report(y_test, y_pred)

print(accuracy)
print(classification_rep)
