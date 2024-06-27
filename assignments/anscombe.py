import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt


df_anscombe = pd.read_csv('/Users/senthilp/Desktop/df_anscombe.csv')
df_anscombe.head()
df_anscombe.group.nunique()
df_anscombe.groupby('group').describe()
df_anscombe.groupby('group').corr()



fig, axs = plt.subplots(2, 2, figsize=(16, 8))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

fig.suptitle("Anscombe's quartet", fontsize=16)

# Plot each group in its own subplot
for i, ax in enumerate(axs.flatten(), 1):
    group_data = df_anscombe[df_anscombe['group'] == i]
    ax.scatter(group_data['x'], group_data['y'])
    ax.set_title(f'Group {i}')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 15)

# Show the plot
plt.show()