import pandas as pd

df = pd.read_csv(r'F:\DEEPTHI\SELF STUDY\sample-data.csv')

# Define a threshold for confident methylation call (e.g., 30 reads)
threshold = 30

# Filter the data for Tissue #2 (cfDNA), assuming 'Tissue' column indicates the tissue type
tissue_df = df[df['Tissue'] == 'cfDNA']


# Create a function to estimate the threshold of reads for each PMP
def estimate_read_threshold(df, threshold):
    # For each unique set of coordinates (PMP), check if the methylation state has sufficient reads
    results = []

    for _, row in df.iterrows():
        methylation_counts = row[['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']].values
        total_reads = sum(methylation_counts)

        # If total reads are below the threshold, add it to the result for further evaluation
        if total_reads < threshold:
            results.append({
                'CpG_Coordinates': row['CpG_Coordinates'],
                'Total_Reads': total_reads,
                'Methylation_Counts': methylation_counts,
                'Required_Reads': threshold - total_reads
            })

    # Return results of PMPs that need more reads
    return pd.DataFrame(results)


# Estimate the threshold of reads required for confident methylation calling
insufficient_reads_df = estimate_read_threshold(tissue_df, threshold)

# Display the results
print(f"PMPs with insufficient reads (below {threshold} reads):")
print(insufficient_reads_df)

# To calculate and display the results for all PMPs (not just below threshold)
all_pmp_df = pd.DataFrame({
    'CpG_Coordinates': tissue_df['CpG_Coordinates'],
    'Total_Reads': tissue_df[['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']].sum(axis=1),
    'Methylation_Counts': tissue_df[['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']].values.tolist()
})

print("\nAll PMPs with total reads:")
print(all_pmp_df)

# Save the results to a CSV
insufficient_reads_df.to_csv('insufficient_reads.csv', index=False)
all_pmp_df.to_csv('all_pmp_reads.csv', index=False)
