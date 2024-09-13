import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

df = pd.read_csv('/Users/sq566/Desktop/mesa-5/res/slppos.csv', sep="\t")
df['N_norm'] = df['N']/32

grouped_sum = df.groupby('VALUE')['N_norm'].sum()


fig, ax = plt.subplots(figsize=(10, 6))
grouped_sum.plot(kind='bar', ax=ax)
ax.set_title('Total seconds (during sleep)')
plt.xticks(rotation=0)  # Keeps the labels horizontal

ax.set_xlabel('')
ax.set_ylabel('Seconds')

positions = ["0-Right", "1-Back", "2-Left", "3-Prone", "4-Upright"]
ax.set_xticks(range(len(positions)))  # Set custom labels
ax.set_xticklabels(positions, rotation=0)
ax.set_ylim(0, 20000000)
ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x))))


plt.grid(False)
plt.savefig('/Users/sq566/Desktop/mesa-5/plots/Position.png', dpi=300)
plt.show()
