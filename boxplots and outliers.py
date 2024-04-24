import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import Normalize, LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.cm import ScalarMappable

file_path = r"C:\Users\Лиза\MSci Project\Expression Data\A_D6_metabolipathway_associated.tsv" #specify path to document with genes associated to a functional annotations, log fold difference of expression and pointwise relative entropy values
data = pd.read_csv(file_path, sep='\t')

# Specify the superpathways to exclude based on binomial results
exclude_superpathways = [
"Metabolism of nucleotides", "Peptide hormone metabolism", "Metabolism of non-coding RNA"
]

# Filter out the specified superpathways
filtered_data = data[~data['Pathway'].isin(exclude_superpathways)]

# Calculate median Pointwise Relative Entropy for ordering
median_entropy = filtered_data.groupby('Pathway')['Pointwise Relative Entropy'].median().reset_index()

# Normalize the median Pointwise Relative Entropy values for color mapping
norm = norm = TwoSlopeNorm(vmin=median_entropy['Pointwise Relative Entropy'].min(), 
                    vmax=median_entropy['Pointwise Relative Entropy'].max(), 
                    vcenter=0.0)

cmap = sns.color_palette("RdBu_r", as_cmap=True)

colors = [cmap(norm(value)) for value in median_entropy['Pointwise Relative Entropy']]
color_map = dict(zip(median_entropy['Pathway'], colors))

fig, ax1 = plt.subplots(figsize=(18, 12))
boxplot = sns.boxplot(x='Log fold difference', y='Pathway', data=filtered_data, palette=color_map, ax=ax1)

ax1.set_title('Boxplot of ...')
ax1.set_xlabel('Log fold difference')
ax1.set_ylabel('Pathway')

for line in boxplot.lines:
    if line.get_marker() == 'o':  # Outliers are marked with 'o' by default
        line.set_marker('x')  # Change to x


sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
colorbar = plt.colorbar(sm, ax=ax1)
min_entropy = median_entropy['Pointwise Relative Entropy'].min()
max_entropy = median_entropy['Pointwise Relative Entropy'].max()

# Set the colorbar ticks to show minimum, middle (0.0), and maximum values
colorbar.set_ticks([min_entropy, 0, max_entropy])

# Optional: Format the tick labels, for example, to show limited decimal places
colorbar.set_ticklabels([f'{min_entropy:.3f}', '0.0', f'{max_entropy:.3f}'])
colorbar.set_label('Pointwise Relative Entropy')
output_file = r"PATH" #specify path for the outputted image

plt.savefig(output_file, bbox_inches='tight')

plt.show()

#Looking at outliers in the box plots

all_outliers = pd.DataFrame()

def find_outliers(data):
    Q1 = data.quantile(0.25)
    Q3 = data.quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return (data < lower_bound) | (data > upper_bound)

for superpathway in filtered_data['Sub-Pathway'].unique():
    superpathway_data = filtered_data[filtered_data['Sub-Pathway'] == superpathway]
    outliers_mask = find_outliers(superpathway_data['Log fold difference'])
    outliers = superpathway_data[outliers_mask]
    all_outliers = pd.concat([all_outliers, outliers])

    if not outliers.empty:
        print(f"Outliers for {superpathway}:")
        print(outliers[['Gene Identifier', 'Log fold difference']])
        print("\n")


all_outliers = all_outliers.drop_duplicates()
total_outliers = len(all_outliers)
unique_gene_identifiers = all_outliers['Gene Identifier'].nunique()

print(f"Total detected outliers across all root superpathways: {total_outliers}")
print(f"Total non-redundant Gene Identifiers across all outliers: {unique_gene_identifiers}")

unique_gene_identifiers_list = all_outliers['Gene Identifier'].unique()
unique_gene_identifiers_list = list(unique_gene_identifiers_list)

print("List of all non-redundant Gene Identifiers across all outliers:")
print(unique_gene_identifiers_list)