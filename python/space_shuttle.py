from ucimlrepo import fetch_ucirepo 
  
# fetch dataset 
challenger_usa_space_shuttle_o_ring = fetch_ucirepo(id=92) 
  
# data (as pandas dataframes) 
dataframe = challenger_usa_space_shuttle_o_ring.data.features 
y = challenger_usa_space_shuttle_o_ring.data.targets 
  
# metadata 
print(challenger_usa_space_shuttle_o_ring.metadata) 
  
# variable information 
print(challenger_usa_space_shuttle_o_ring.variables) 



import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.scatter(dataframe['launch_temp'], dataframe['num_thermal_distress'], color='blue', label='Data Points')
plt.xlabel('Launch Temperature (°F)')
plt.ylabel('Number of O-rings in Thermal Distress')
plt.xlim(20, 100)  # Set x-axis limits
plt.ylim(0, 6)    # Set y-axis limits
plt.grid(True)
plt.legend()
plt.savefig('/Users/sq566/Desktop/Coursera/thermal_distress_plot1.png', dpi=300)
plt.show()




from sklearn.linear_model import LogisticRegression
import numpy as np

# Define the independent variable (launch_temp) and the dependent variable (num_thermal_distress)
X = dataframe[['launch_temp']]
y = dataframe['num_thermal_distress']

# Fit the logistic regression model
logistic_model = LogisticRegression(max_iter=1000)
logistic_model.fit(X, y)

# Prepare new data for prediction
new_data = np.array([[20], [25], [30], [35], [40], [45], [50], [55], [60], [65], [70], [75],
                     [80], [85], [90],[95], [100]])

# Predict the probability of thermal distress
predicted_probs = logistic_model.predict_proba(new_data)[:, 1]

# Calculate the expected number of O-rings experiencing thermal distress
predicted_num_thermal_distress = predicted_probs * 6

# Round the predictions to the nearest whole number
rounded_predictions = np.round(predicted_num_thermal_distress).astype(int)

rounded_predictions


temp_range = np.linspace(20, 100, 17).reshape(-1, 1)
plt.figure(figsize=(10, 6))
plt.scatter(dataframe['launch_temp'], dataframe['num_thermal_distress'], color='blue', label='Data Points')
plt.plot(temp_range, predicted_probs * 6, color='red', linewidth=2, label='Fitted Model')
plt.scatter([20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100], rounded_predictions, color='green', marker='x', s=100, label='Predicted Data Points')
plt.xlabel('Launch Temperature (F)')
plt.ylabel('Number of O-rings in Thermal Distress')
plt.title('Fitted Binomial Logistic Regression Model')
plt.xlim(20, 100)
plt.ylim(0, 6)
plt.legend()
plt.grid(True)
plt.savefig('/Users/sq566/Desktop/Coursera/thermal_distress_plot2.png', dpi=300)
plt.show()


'''
The probabilities are multiplied by 6 to estimate the expected number of O-rings experiencing thermal distress for each launch temperature. 
Here's the reasoning behind this step:

Logistic Regression Prediction: Logistic regression provides a probability for each data point that an event (thermal distress) will occur. 
In this case, the model predicts the probability that thermal distress will occur for a given launch temperature.

O-rings in the Shuttle: The context suggests that there are 6 O-rings in the shuttle's solid rocket boosters. 
Therefore, for each launch temperature, the model predicts the probability that any one O-ring will experience thermal distress.

Expected Number of Distressed O-rings: To estimate the expected number of O-rings in distress, 
you multiply the probability of distress (from the logistic regression model) by the total number of O-rings (which is 6). 
This gives an expected value for the number of O-rings in distress at a given launch temperature.

For example:

If the model predicts a 0.5 probability of distress at a certain temperature, multiplying this by 6 gives an expected number of 3 O-rings in distress.
If the probability is 0.2, the expected number of O-rings in distress is 0.2 × 6 = 1.2
This approach provides a way to interpret the probabilistic output of the logistic regression model in terms of the actual 
number of O-rings that might be affected, making it more practical and understandable for decision-making.
'''



