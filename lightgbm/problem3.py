import pandas as pd
import lightgbm as lgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.preprocessing import LabelEncoder

file_path = '/Users/senthilp/Desktop/ML_Dataset.csv'
data = pd.read_csv(file_path)

label_encoder_gene = LabelEncoder()
label_encoder_mutation = LabelEncoder()
label_encoder_effect = LabelEncoder()

data['Gene_encoded'] = label_encoder_gene.fit_transform(data['Gene'])
data['Mutation_encoded'] = label_encoder_mutation.fit_transform(data['Mutation'])
data['Effect_encoded'] = label_encoder_effect.fit_transform(data['Effect'])

X_lgb = data[['Gene_encoded', 'Mutation_encoded', 'Effect_encoded', 'Location']]
y_lgb = data['ClinicalOutcome']

# Train-test split
X_train_lgb, X_test_lgb, y_train_lgb, y_test_lgb = train_test_split(X_lgb, y_lgb, test_size=0.2, random_state=42)

# Create LightGBM dataset
lgb_train = lgb.Dataset(X_train_lgb, y_train_lgb)
lgb_test = lgb.Dataset(X_test_lgb, y_test_lgb, reference=lgb_train)

# Set up LightGBM parameters
params = {
    'objective': 'binary',
    'metric': 'auc',
    'boosting_type': 'gbdt',
    'num_leaves': 31,
    'learning_rate': 0.05,
}

callbacks = [lgb.early_stopping(stopping_rounds=10)]

# Train LightGBM model
lgbm_model = lgb.train(params, lgb_train, num_boost_round=100, valid_sets=[lgb_train, lgb_test], callbacks=callbacks)

# Make predictions
y_pred_lgb = (lgbm_model.predict(X_test_lgb) > 0.5).astype(int)

# Evaluate the model
print(classification_report(y_test_lgb, y_pred_lgb))
