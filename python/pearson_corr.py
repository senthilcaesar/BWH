import numpy as np
from scipy.stats import pearsonr

# Given data for salt intake (in mEq) and systolic blood pressure (in mmHg)
salt_intake = np.array([106.0960, 194.7779, 275.2025, 397.4523, 497.3065, 574.1339, 705.6480,
                        801.5520, 881.2873, 999.4862])

blood_pressure = np.array([100.99, 105.58, 114.04, 114.79, 115.99, 117.13, 122.20, 124.84, 126.01, 129.70])

# Calculate the Pearson correlation coefficient
correlation_coefficient, p_value = pearsonr(salt_intake, blood_pressure)

print("\n")
print(f"The Pearson correlation coefficient is approximately {correlation_coefficient}, which indicates\n"
      "a strong positive correlation between salt intake and systolic blood pressure.")
print("\n")
print(f"The p-value is extremely small {p_value}, suggesting that the correlation is statistically significant.")
print("\n")
