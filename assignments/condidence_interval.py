import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
mu = 0
sigma = 1
confidence_level = 0.95
z_value = 1.96  # Critical value for 95% confidence
sample_sizes = [10, 50, 100, 200, 300, 400, 500]
num_intervals = 100

# Function to calculate confidence intervals
def calculate_confidence_intervals(sample_size):
    sample_means = []
    intervals = []

    for _ in range(num_intervals):
        sample = np.random.normal(mu, sigma, sample_size)
        sample_mean = np.mean(sample)
        sem = sigma / np.sqrt(sample_size)
        margin_of_error = z_value * sem
        interval = (sample_mean - margin_of_error, sample_mean + margin_of_error)
        sample_means.append(sample_mean)
        intervals.append(interval)

    return sample_means, intervals

# Calculate confidence intervals for each sample size
results = {}
for size in sample_sizes:
    sample_means, intervals = calculate_confidence_intervals(size)
    results[size] = intervals

# Create a DataFrame to visualize the intervals
df_intervals = pd.DataFrame({
    "Sample Size": np.repeat(sample_sizes, num_intervals),
    "Lower Bound": [interval[0] for size in sample_sizes for interval in results[size]],
    "Upper Bound": [interval[1] for size in sample_sizes for interval in results[size]]
})


# Plotting the confidence intervals
plt.figure(figsize=(12, 8))

for size in sample_sizes:
    intervals = results[size]
    plt.plot([size] * num_intervals, [interval[0] for interval in intervals], 'r_', markersize=10)
    plt.plot([size] * num_intervals, [interval[1] for interval in intervals], 'b_', markersize=10)

plt.xlabel('Sample Size')
plt.ylabel('Confidence Interval')
plt.title('95% Confidence Intervals for Different Sample Sizes')
plt.xticks(sample_sizes)
plt.grid(True)
plt.savefig('/Users/sq566/Desktop/confidence_inetrval_plot.png', dpi=300)
plt.show()

