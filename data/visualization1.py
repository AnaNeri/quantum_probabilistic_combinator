import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Load the data
df1 = pd.read_csv('out_test2_1000_002.csv')
df4 = pd.read_csv('out_test2_1000_correction_less_50percent.csv')  
df2 = pd.read_csv('out_test2_1000_correction_less_75percent.csv')
df3 = pd.read_csv('out_test2_1000_no_noise_in_correction.csv')

# Create a 2x2 grid for the plots
fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharey=True)

# Plot for df1
axes[0, 0].plot(df1['prob_error'], df1['0_corr'], label="0 with correction", color='blue')
axes[0, 0].plot(df1['prob_error'], df1['0_no_corr'], label="0 without correction", color='red')
axes[0, 0].set_title("Test2 Results (1000 trials)")
axes[0, 0].set_xlabel('Error Probability')
axes[0, 0].set_ylabel('Probability of 0')
axes[0, 0].legend()
axes[0, 0].grid(True)

# Plot for df4 (new graph)
axes[0, 1].plot(df4['prob_error'], df4['0_corr'], label="0 with correction", color='blue')
axes[0, 1].plot(df4['prob_error'], df4['0_no_corr'], label="0 without correction", color='red')
axes[0, 1].set_title("Test2 Results - Noise is 50% less in correction")
axes[0, 1].set_xlabel('Error Probability')
axes[0, 1].legend()
axes[0, 1].grid(True)

# Plot for df2
axes[1, 0].plot(df2['prob_error'], df2['0_corr'], label="0 with correction", color='blue')
axes[1, 0].plot(df2['prob_error'], df2['0_no_corr'], label="0 without correction", color='red')
axes[1, 0].set_title("Test2 Results - Noise is 75% less in correction")
axes[1, 0].set_xlabel('Error Probability')
axes[1, 0].set_ylabel('Probability of 0')
axes[1, 0].legend()
axes[1, 0].grid(True)

# Plot for df3
axes[1, 1].plot(df3['prob_error'], df3['0_corr'], label="0 with correction", color='blue')
axes[1, 1].plot(df3['prob_error'], df3['0_no_corr'], label="0 without correction", color='red')
axes[1, 1].set_title("Test2 Results - No Noise in Correction")
axes[1, 1].set_xlabel('Error Probability')
axes[1, 1].legend()
axes[1, 1].grid(True)

# Adjust layout and save the figure
plt.tight_layout()
now = datetime.now().strftime("%Y%m%d%H%M")
plt.savefig(f"{now}_comparison.png")
plt.show()