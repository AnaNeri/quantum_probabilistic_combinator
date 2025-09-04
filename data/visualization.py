import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df1 = pd.read_csv('out_test2_1500_003.csv')
df2 = pd.read_csv('out_test2_1000_002.csv')

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

axes[0].plot(df1['prob_error'], df1['0_corr'], label="0 with correction", color='blue')
axes[0].plot(df1['prob_error'], df1['0_no_corr'], label="0 without correction", color='red')
axes[0].set_title("Test2 Results (1500 trials)")
axes[0].set_xlabel('Error Probability')
axes[0].set_ylabel('Probability of 0')
axes[0].legend()
axes[0].grid(True)

axes[1].plot(df2['prob_error'], df2['0_corr'], label="0 with correction", color='blue')
axes[1].plot(df2['prob_error'], df2['0_no_corr'], label="0 without correction", color='red')
axes[1].set_title("Test2 Results (1000 trials)")
axes[1].set_xlabel('Error Probability')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()