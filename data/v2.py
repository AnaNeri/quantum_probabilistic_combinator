import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

df = pd.read_csv('data/out_test2_10000_bitFlip_0.0.csv')


plt.plot(df['prob_error'], df['0_corr'], label="0 with correction", color='blue')
plt.plot(df['prob_error'], df['0_no_corr'], label="0 without correction", color='red')

plt.title("Test Results - Noisy correction")
plt.xlabel('Error Probability')
plt.ylabel('Probability of 0')
plt.legend()
plt.grid(True)

plt.tight_layout()
now = datetime.now().strftime("%Y%m%d%H%M")
plt.savefig(f"{now}_comparison.png")
plt.show()