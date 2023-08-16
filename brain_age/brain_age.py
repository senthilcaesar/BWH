
import pickle
import numpy as np

# Input Features
X = [8.037083723,7.080192415,2.94300045,0.843729859,0.871010226, 
      0.903848936,5.188447503,22.65975292,0.918476141,276.5,2.67833698,
      21.93443638,0.025051651]

  
X = np.array(X)
X = X.reshape(1, 13)

f_model = '/Users/sq566/Desktop/Haoqi/SleepEEGBasedBrainAge/brain_age_model_c/stable_brain_age_model.pickle'

with open(f_model, 'rb') as f1:
     model = pickle.load(f1)

print(X)

coef = model['predictor'].coef_
intercept = model['predictor'].intercept_

age = 55
BA = model.predict(X, y=age)[0]
BAI = BA-age

print("\n")
print(BA, BAI)
print("\n")

X = model._validate_data(X, accept_sparse=["csr", "csc", "coo"], reset=False)


from sklearn import linear_model
clf = linear_model.BayesianRidge()

def _decision_function(model, X):

    coef = model['predictor'].coef_
    intercept = model['predictor'].intercept_
    X = model._validate_data(X, accept_sparse=["csr", "csc", "coo"], reset=False)
    
    return safe_sparse_dot(X, coef_.T, dense_output=True) + intercept




yp = 50.13990535479082
yp2 = yp - (model['bias_correction'].correction_slope * age + model['bias_correction'].correction_intercept)
yp2 = np.log1p(np.exp(-np.abs(X2))) + np.maximum(X2,0)



Original data
[[8.03708372e+00 7.08019241e+00 2.94300045e+00 8.43729859e-01
  8.71010226e-01 9.03848936e-01 5.18844750e+00 2.26597529e+01
  9.18476141e-01 2.76500000e+02 2.67833698e+00 2.19344364e+01
  2.50516507e-02]]

Transformaed data
[[ 0.22501175 -0.09673319  0.06812488 -0.83194508 -0.86245684 -0.87465325
  -0.08597223  1.3859261  -0.29070937 -0.46841362 -0.96005772  1.41872477
  -0.90946905]]




