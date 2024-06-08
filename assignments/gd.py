import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

def pred_sklearn(X, lr_sklearn):
    X_2D = X[:, np.newaxis]
    Y = lr_sklearn.predict(X_2D)    
    return Y

adv = pd.read_csv("/Users/sq566/Desktop/tvmarketing.csv")
adv.plot(x='TV', y='Sales', kind='scatter', c='black')

lr_sklearn = LinearRegression()

X = np.array(adv['TV'])
Y = np.array(adv['Sales'])

X_sklearn = X[:, np.newaxis]
Y_sklearn = Y[:, np.newaxis]

print(f"Shape of new X array: {X_sklearn.shape}")
print(f"Shape of new Y array: {Y_sklearn.shape}")

lr_sklearn.fit(X_sklearn, Y_sklearn)

m_sklearn = lr_sklearn.coef_[0][0]  # slope
b_sklearn = lr_sklearn.intercept_[0] # intercept

print(f"Linear regression using Scikit-Learn. Slope: {m_sklearn}. Intercept: {b_sklearn}")

X_pred = np.array([10, 50, 120, 300])
Y_pred = pred_sklearn(X_pred, lr_sklearn)
print(f"TV marketing expenses:\n{X_pred}")
print(f"Predictions of sales using Scikit_Learn linear regression:\n{Y_pred.T}")


plt.figure(figsize=(10, 6))
plt.scatter(X, Y, color='blue', label='Original data')
plt.scatter(X_pred, Y_pred, color='green', label='Predicted data')
plt.plot(X_pred, Y_pred, color='red', label=f'Fitted line (y = {m_sklearn:.2f}x + {b_sklearn:.2f})', linewidth=2)
plt.xlabel('X')
plt.ylabel('y')
plt.title('Linear Regression Fit')
plt.legend()
plt.grid(True)
plt.show()



# Linear Regression Gradient Descent
X_norm = (X - np.mean(X))/np.std(X)
Y_norm = (Y - np.mean(Y))/np.std(Y)

def E(m, b, X, Y):
    return 1/(2*len(Y))*np.sum((m*X + b - Y)**2)


def dEdm(m, b, X, Y):
    res = 1/len(X)*np.sum(np.dot(m * X + b - Y, X))    
    return res
    
def dEdb(m, b, X, Y):
    res = 1/len(X)*np.sum((m * X + b - Y))    
    return res


def gradient_descent(dEdm, dEdb, m, b, X, Y, learning_rate = 0.001, num_iterations = 1000, print_cost=False):
    for iteration in range(num_iterations):
        m_new = m - learning_rate * dEdm(m, b, X, Y)
        b_new = b - learning_rate * dEdb(m, b, X, Y)
        m = m_new
        b = b_new
        if print_cost:
            print (f"Cost after iteration {iteration}: {E(m, b, X, Y)}")
    return m, b



m_initial = 0; b_initial = 0; num_iterations = 30; learning_rate = 1.2
m_gd, b_gd = gradient_descent(dEdm, dEdb, m_initial, b_initial, 
                              X_norm, Y_norm, learning_rate, num_iterations, print_cost=True)

print(f"Gradient descent result: m_min, b_min = {m_gd}, {b_gd}") 

X_pred = np.array([10, 50, 120, 300])
X_pred_norm = (X_pred - np.mean(X))/np.std(X)
Y_pred_gd_norm = m_gd * X_pred_norm + b_gd
Y_pred_gd = Y_pred_gd_norm * np.std(Y) + np.mean(Y)

print(f"TV marketing expenses:\n{X_pred}")
print(f"Predictions of sales using Scikit_Learn linear regression:\n{Y_pred.T}")
print(f"Predictions of sales using Gradient Descent:\n{Y_pred_gd}")















