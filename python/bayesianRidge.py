import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import BayesianRidge

# Generate synthetic data
np.random.seed(42)
X = np.sort(5 * np.random.rand(80, 1), axis=0)
y_true = np.sin(X).ravel()
y_obs = y_true + 0.1 * np.random.randn(80)

# Fit Bayesian Ridge Regression
model = BayesianRidge(compute_score=True)
model.fit(X, y_obs)

# Predictions
X_test = np.linspace(0, 5, 1000)[:, np.newaxis]
y_pred, y_std = model.predict(X_test, return_std=True)

# Plot the true function, observed data, and Bayesian Ridge Regression predictions
plt.figure(figsize=(12, 6))
plt.errorbar(X, y_obs, 0.1, fmt='r.', markersize=10, label='Observations with Noise')
plt.plot(X_test, y_true, 'k', label='True Function')
plt.plot(X_test, y_pred, 'b', label='Bayesian Ridge Regression Prediction')
plt.fill_between(X_test.ravel(), y_pred - y_std, y_pred + y_std, alpha=0.2, color='blue')

plt.title('Bayesian Ridge Regression Example')
plt.xlabel('X')
plt.ylabel('y')
plt.legend()
plt.show()
