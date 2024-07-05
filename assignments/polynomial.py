import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import Ridge


# Generate synthetic data
np.random.seed(0)  # For reproducibility
n_points = 1000
x = np.linspace(0, 1, n_points)
y = np.sin(2 * np.pi * x) + np.random.normal(0, 0.1, n_points)  # Adding some noise

# Generate polynomial features
degree = 9
poly = PolynomialFeatures(degree)
X_poly = poly.fit_transform(x.reshape(-1, 1))

# Fit a ridge regression model
ridge_reg = Ridge(alpha=1e-3)  # Regularization term
ridge_reg.fit(X_poly, y)

# Generate fitted values
y_ridge_fit = ridge_reg.predict(X_poly)

# Plot the original data and the ridge regression polynomial fit
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='blue', label='Synthetic Data')
plt.plot(x, y_ridge_fit, color='red', label=f'Ridge Polynomial Fit (degree={degree})')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Ridge Polynomial Fit to Sinusoidal Data')
plt.legend()
plt.grid(True)
plt.show()

# Display the coefficients of the polynomial
ridge_coefficients = ridge_reg.coef_

# Since the first coefficient is the intercept, we remove it for clarity
ridge_coefficients = ridge_coefficients[1:]

print(ridge_coefficients)