import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# Sample data
data = [8.7, 14.2, 18.3, 18.4, 23.2, 25.9, 29.7, 35.2, 51.2, 54.7, 65.9, 75]

# Calculate Q1, Q2 (median), and Q3
Q1 = np.percentile(data, 25)
Q2 = np.median(data)
Q3 = np.percentile(data, 75)

iqr = Q3 - Q1
print("Interquartile range (IQR):", iqr)

middle_50_percent = [x for x in data if Q1 <= x <= Q3]
middle_50_percent

print(f'The points within the middle 50% of the data, which lie between Q1 ({Q1}) and Q3 ({Q3}), are: {middle_50_percent}')


plt.figure(figsize=(10, 6))
plt.boxplot(data, vert=False)

# Adding annotations for Q1, Q2, and Q3
plt.text(Q1, 1.05, 'Q1', horizontalalignment='center', fontsize=12, color='blue')
plt.text(np.median(data), 1.05, 'Q2 (Median)', horizontalalignment='center', fontsize=12, color='blue')
plt.text(Q3, 1.05, 'Q3', horizontalalignment='center', fontsize=12, color='blue')

# Marking the points in the middle 50%
for point in data:
    if Q1 <= point <= Q3:
        plt.plot(point, 1, 'ro')

plt.title('Box Plot')
plt.xlabel('Value')
plt.grid(True)
plt.savefig('/Users/senthilp/Desktop/box.png', dpi=300)
plt.show()


# Create a violin plot
plt.figure(figsize=(10, 7))
sns.violinplot(data=data, inner=None)

# Overlay the box plot with Q1, Median, Q3, and Whiskers
sns.boxplot(data=data, whis=1.5, width=0.2, showcaps=True, showbox=True, showfliers=False,
            whiskerprops={'linewidth': 2}, medianprops={'color': 'yellow', 'linewidth': 2})

# Annotate Q1, Median, Q3
plt.text(Q1, 0.15, 'Q1', horizontalalignment='center',
         fontsize=12, color='blue')
plt.text(np.median(data), 0.15, 'Median',
         horizontalalignment='center', fontsize=12, color='blue')
plt.text(Q3, 0.15, 'Q3', horizontalalignment='center',
         fontsize=12, color='blue')

plt.title('Violin Plot with Q1, Median, Q3, Whiskers, and KDE Curve')
plt.xlabel('Value')
plt.ylabel('Density')
plt.grid(True)

# Save the plot
plt.savefig('/Users/senthilp/Desktop/violin.png', dpi=300)
plt.show()



import scipy.stats as stats

mean = np.mean(data)
median = np.median(data)
std_dev = np.std(data)
variance = np.var(data)
percentiles = np.percentile(data, [25, 50, 75])

# Plotting Histogram
plt.figure(figsize=(12, 6))
sns.histplot(data, bins=10, kde=True, color='blue')
plt.title('Histogram with KDE')
plt.xlabel('Value')
plt.ylabel('Frequency')

# Plotting Normal Distribution for Comparison
plt.figure(figsize=(12, 6))
norm_data = np.random.normal(mean, std_dev, 1000)
sns.histplot(norm_data, bins=10, kde=True, color='green')
plt.title('Histogram of Normal Distribution with Same Mean and Std Dev')
plt.xlabel('Value')
plt.ylabel('Frequency')

# Q-Q Plot
plt.figure(figsize=(12, 6))
stats.probplot(data, dist="norm", plot=plt)
plt.title('Q-Q Plot')

plt.show()



