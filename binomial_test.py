import pandas as pd
from scipy.stats import binom_test
import numpy as np

file_path = r"PATH" #Path to file containing superpathway counts and empirical probabilities for input and reference gene sets
df = pd.read_excel(file_path) 

# choose significance lvl
alpha = 0.05

results = []

for index, row in df.iterrows():
    successes = row['D6 Count'] #row showing number of counts assopciated with each pathway detected in input gene set
    observed_prob = row['CSR Probability'] #row with emperical probabilities of the refernce gene set

    if pd.notna(observed_prob):
        total_trials = 1597 #number of genes in the input gene set 

        p_value = binom_test(successes, total_trials, observed_prob)

        log_p_value = -np.log10(p_value)

        results.append({
            'Pathway': row['DisplayName'],
            'Successes': successes,
            'Observed Probability': observed_prob,
            'P-value': p_value,
            'Log of P-value': log_p_value,
            'Significant': p_value < alpha
        })

results_df = pd.DataFrame(results)

print(results_df)
results_df.to_csv(r"PATH", sep="\t", index=False) #specify path to output



