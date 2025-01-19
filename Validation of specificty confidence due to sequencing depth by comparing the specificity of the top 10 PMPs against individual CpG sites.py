import pandas as pd
import numpy as np

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(r'F:\DEEPTHI\SELF STUDY\sample-data.csv')

# Define a threshold for confident methylation call (e.g., 30 reads)
threshold = 30

# Filter the data for Tissue #2 (cfDNA), assuming 'Tissue' column indicates the tissue type
tissue_df = df[df['Tissue'] == 'cfDNA']

# Define methylation columns
methylation_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110',
                       '`111']  # Modify if needed based on the actual column names


# Function to calculate specificity for individual CpG sites
def calculate_specificity_for_sites(df, methylation_columns):
    # For each CpG site, calculate the distribution of methylation states
    specificity_results = []

    for _, row in df.iterrows():
        methylation_counts = row[methylation_columns].values
        total_reads = sum(methylation_counts)

        # Specificity: proportion of methylation states with sufficient reads
        if total_reads >= threshold:
            specificity = np.sum(methylation_counts > 0) / len(
                methylation_columns)  # Proportion of methylation states observed
            specificity_results.append({
                'CpG_Coordinates': row['CpG_Coordinates'],
                'Specificity': specificity,
                'Total_Reads': total_reads,
                'Methylation_Counts': methylation_counts
            })

    return pd.DataFrame(specificity_results)


# Function to calculate specificity for PMPs (combination of three CpG sites)
def calculate_specificity_for_PMPs(df, methylation_columns):
    # For each PMP (combination of CpG sites), calculate the distribution of methylation states
    pmp_specificity_results = []

    for _, row in df.iterrows():
        methylation_counts = row[methylation_columns].values
        total_reads = sum(methylation_counts)

        # Specificity: Proportion of observed methylation patterns for the combined sites
        if total_reads >= threshold:
            specificity = np.sum(methylation_counts > 0) / len(
                methylation_columns)  # Proportion of methylation states observed
            pmp_specificity_results.append({
                'CpG_Coordinates': row['CpG_Coordinates'],
                'Specificity': specificity,
                'Total_Reads': total_reads,
                'Methylation_Counts': methylation_counts
            })

    return pd.DataFrame(pmp_specificity_results)


# Calculate specificity for individual CpG sites
individual_specificity_df = calculate_specificity_for_sites(tissue_df, methylation_columns)

# Calculate specificity for the top 10 PMPs (here, we use the first 10 rows as an example, but you can sort based on any criteria like methylation counts)
top_10_pmp_df = tissue_df.head(10)
pmp_specificity_df = calculate_specificity_for_PMPs(top_10_pmp_df, methylation_columns)

# Compare Specificity between PMPs and individual CpG sites
# Merge the two DataFrames on CpG_Coordinates (if the data allows, or just compare based on order)
comparison_df = pd.merge(individual_specificity_df, pmp_specificity_df, on='CpG_Coordinates',
                         suffixes=('_individual', '_pmp'))

# Output the comparison for analysis
print("Specificity Comparison between Individual CpG Sites and PMPs:")
print(comparison_df)

# Save the results to CSV files if needed
individual_specificity_df.to_csv('individual_cpg_specificity.csv', index=False)
pmp_specificity_df.to_csv('pmp_specificity.csv', index=False)
comparison_df.to_csv('specificity_comparison.csv', index=False)
