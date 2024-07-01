import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm

mu = 10
sigma = 5

gaussian_population = np.random.normal(mu, sigma, 100_000)

sns.histplot(gaussian_population, stat="density")
plt.show()

gaussian_pop_mean = np.mean(gaussian_population)
gaussian_pop_std = np.std(gaussian_population)

print(f"Gaussian population has mean: {gaussian_pop_mean:.1f} and std: {gaussian_pop_std:.1f}")


def sample_means(data, sample_size):
    # Save all the means in a list
    means = []

    # For a big number of samples
    # This value does not impact the theorem but how nicely the histograms will look (more samples = better looking)
    for _ in range(10_000):
        # Get a sample of the data WITH replacement
        sample = np.random.choice(data, size=sample_size)

        # Save the mean of the sample
        means.append(np.mean(sample))

    # Return the means within a numpy array
    return np.array(means)



# Compute the sample means
gaussian_sample_means = sample_means(gaussian_population, sample_size=5)

# Plot a histogram of the sample means
sns.histplot(gaussian_sample_means, stat="density")
plt.show()



# Compute estimated mu
mu_sample_means = mu

# Compute estimated sigma
# 5 is being used because you used a sample size of 5
sigma_sample_means = sigma / np.sqrt(5)

# Define the x-range for the Gaussian curve (this is just for plotting purposes)
x_range = np.linspace(min(gaussian_sample_means), max(gaussian_sample_means), 100)

# Plot everything together
sns.histplot(gaussian_sample_means, stat="density")
plt.plot(
    x_range,
    norm.pdf(x_range, loc=mu_sample_means, scale=sigma_sample_means),
    color="black",
)
plt.show()



# Histogram of sample means (blue)
sns.histplot(gaussian_sample_means, stat="density", label="hist")

# Estimated PDF of sample means (red)
sns.kdeplot(
    data=gaussian_sample_means,
    color="crimson",
    label="kde",
    linestyle="dashed",
    fill=True,
)

# Gaussian curve with estimated mu and sigma (black)
plt.plot(
    x_range,
    norm.pdf(x_range, loc=mu_sample_means, scale=sigma_sample_means),
    color="black",
    label="gaussian",
)

plt.legend()
plt.show()



# Create the QQ plot
fig, ax = plt.subplots(figsize=(6, 6))
res = stats.probplot(gaussian_sample_means, plot=ax, fit=True)
plt.show()









