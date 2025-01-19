import pandas as pd

# Load your data (replace 'your_file.tsv' with the actual file path)
# Assuming the file is tab-delimited; if it's CSV, use read_csv instead of read_table
df = pd.read_csv(r'F:\DEEPTHI\SELF STUDY\sample-data.csv')

# List of columns that contain the variant read counts (the binary columns)
variant_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']

# Calculate the total number of reads for each row
df['Total_Reads'] = df[variant_columns].sum(axis=1)

# Calculate the number of variant reads (sum of the variant columns)
df['Variant_Reads'] = df[variant_columns].sum(axis=1)

# Calculate the Variant Read Fraction (VRF) for each row
df['VRF'] = df['Variant_Reads'] / df['Total_Reads']

# Group by the CpG_Coordinates and calculate the mean VRF for each unique CpG
mean_vrf_by_cpg = df.groupby('CpG_Coordinates')['VRF'].mean()

# Print the result
print(mean_vrf_by_cpg)

# Optionally, save the result to a file
mean_vrf_by_cpg.to_csv('mean_vrf_by_cpg.csv', header=True)
