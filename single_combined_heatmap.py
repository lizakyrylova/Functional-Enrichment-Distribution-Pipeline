import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def add_significance_hatches(ax, significance_data, hatch='///'):
    for y in range(len(significance_data)):
        for x in range(len(significance_data.columns)):
            if significance_data.iat[y, x] == 0: 
                ax.add_patch(Rectangle((x, y), 1, 1, fill=False, hatch=hatch, edgecolor='grey', lw=0))

# Load pointwise relative entropy data and binomial significance data
df_PRE = pd.read_excel(r"C:\Users\Лиза\MSci Project\Level 2 Annotations\A_metabolism_nucleo_dig\Root to Sub\heatmapdata.xlsx")
df_binomial = pd.read_csv(r"C:\Users\Лиза\MSci Project\Level 2 Annotations\A_metabolism_nucleo_dig\Root to Sub\binomial.tsv", sep='\t')

df_binomial_hatching = df_binomial.set_index('Pathway').filter(regex='Significant').astype(int)

# Find overall min and max values for color scaling
overall_min = df_PRE.iloc[:, 1:].min().min()
overall_max = df_PRE.iloc[:, 1:].max().max()

cmap = sns.color_palette("RdBu_r", as_cmap=True)
fig, ax = plt.subplots(figsize=(12, 8))

# Plotting the heatmap
sns.heatmap(data=df_PRE.set_index('Pathway'), cmap=cmap, annot=True, fmt=".4f",
            linewidths=.5, vmin=overall_min, vmax=overall_max, center=0, cbar=True, cbar_kws={'label': 'Values'}, ax=ax)

# Add hatches for statistical significance on the CSR heatmap
add_significance_hatches(ax, df_binomial_hatching)

plt.xticks(rotation=0)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

ax.set_xlabel('Timepoint', fontsize=13)
ax.set_ylabel('Pathway', fontsize=13)
# Adjust layout and title
plt.suptitle('Heatmap of ...')
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save the figure
output_file = r"PATH"
plt.savefig(output_file, bbox_inches='tight')

# Show the plot
plt.show()