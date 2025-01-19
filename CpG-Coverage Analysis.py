import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv(r'C:\Users\iucbr\Downloads\Pupil-Bio-PMP-Dataset\sample-data.csv')

coverage_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
print(coverage_columns)

grouped = data.groupby('Tissue')[coverage_columns]

# Calculate median
medians = grouped.median()

# Calculate coefficient of variation (CV) = (std / mean) * 100
means = grouped.mean()
std_devs = grouped.std()

cv = (std_devs / means) * 100

# Set the aesthetic style of the plots
sns.set(style="whitegrid")

# Define colors for each coverage category (this can be customized)
coverage_colors = sns.color_palette("Set2", len(coverage_columns))

# --- Plot for Median CpG Coverage ---
plt.figure(figsize=(12, 7))

# Loop through the coverage categories and plot each one with a unique color
bar_width = 0.1  # Width of the bars

for i, coverage_column in enumerate(coverage_columns):
    # Slightly offset the x-position of each bar to avoid overlap
    x_positions = np.arange(len(medians)) + i * bar_width
    plt.bar(x_positions, medians[coverage_column],
            color=coverage_colors[i],
            label=coverage_column,
            width=bar_width,  # Narrow bars to fit multiple
            edgecolor='black')

plt.title('Median CpG Coverage for Each Tissue', fontsize=16, weight='bold')
plt.xlabel('Tissue', fontsize=12)
plt.ylabel('Median Coverage', fontsize=12)
plt.xticks(np.arange(len(medians)) + bar_width * (len(coverage_columns) // 2),  # Adjust x-ticks to center on tissue
           medians.index, rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)
plt.grid(True, axis='y', linestyle='--', alpha=0.7)

# Add a legend to identify each coverage category
plt.legend(title="Coverage Categories", fontsize=10)

# Add value labels on the bars
for i, coverage_column in enumerate(coverage_columns):
    for j, v in enumerate(medians[coverage_column]):
        plt.text(j + i * bar_width, v + 0.5, f'{v:.0f}', ha='center', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.show()

sns.set(style="whitegrid")

# Define colors for each coverage category
coverage_colors = sns.color_palette("Set2", len(coverage_columns))

# --- Plot for Median CpG Coverage ---
plt.figure(figsize=(12, 7))

bar_width = 0.1  # Width of the bars

# Loop through the coverage categories and plot each one with a unique color
for i, coverage_column in enumerate(coverage_columns):
    # Slightly offset the x-position of each bar to avoid overlap
    x_positions = np.arange(len(medians)) + i * bar_width
    plt.bar(x_positions, medians[coverage_column],
            color=coverage_colors[i],
            label=coverage_column,
            width=bar_width,  # Narrow bars to fit multiple
            edgecolor='black')

# Title and Labels
plt.title('Median CpG Coverage for Each Tissue', fontsize=16, weight='bold')
plt.xlabel('Tissue', fontsize=12)
plt.ylabel('Median Coverage', fontsize=12)
plt.xticks(np.arange(len(medians)) + bar_width * (len(coverage_columns) // 2),  # Adjust x-ticks to center on tissue
           medians.index, rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)
plt.grid(True, axis='y', linestyle='--', alpha=0.7)

# Add a legend to identify each coverage category
plt.legend(title="Coverage Categories", fontsize=10)

# Add value labels on the bars
for i, coverage_column in enumerate(coverage_columns):
    for j, v in enumerate(medians[coverage_column]):
        # Adjust y-position slightly above the bar height
        plt.text(j + i * bar_width, v + 0.5, f'{v:.0f}', ha='center', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.show()