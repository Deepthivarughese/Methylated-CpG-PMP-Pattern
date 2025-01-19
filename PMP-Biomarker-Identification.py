import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.preprocessing import LabelEncoder

# Step 1: Import the CSV file
# Replace 'your_data.csv' with the path to your CSV file
df = pd.read_csv(r'F:\DEEPTHI\SELF STUDY\sample-data.csv')
df.columns = df.columns.str.strip()  # Removes extra spaces from column names
#print("Column names:", df.columns)
#
# # Step 2: Data Preprocessing
# # Convert methylation columns ('000', '001', '010', etc.) into integer format
methylation_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
#print (methylation_columns[2])
df[methylation_columns] = df[methylation_columns].astype(int)

# # Convert 'strand' to numerical values: 'r' -> 1, 'f' -> 0
df['strand'] = df['strand'].map({'r': 1, 'f': 0})
#
# # Encode 'Tissue' (cfDNA, or any other tissue types) using LabelEncoder
label_encoder = LabelEncoder()
df['Tissue'] = label_encoder.fit_transform(df['Tissue'])
#
# # Step 3: Statistical Testing (Chi-square or Fisher’s Exact Test)
# # For each methylation pattern column, test its association with Tissue type
#
p_values = {}
for methyl_col in methylation_columns:
    # Create a contingency table
    contingency_table = pd.crosstab(df[methyl_col], df['Tissue'])

    # Check the shape of the contingency table
    print(f"\nContingency table for {methyl_col}:")
    print(contingency_table)

    if contingency_table.shape == (2, 2):
        # If the table is 2x2, use Fisher's Exact Test
        _, p_value = fisher_exact(contingency_table)
    else:
        # Otherwise, use Chi-square test
        _, p_value, _, _ = chi2_contingency(contingency_table)

    # Store the p-values
    p_values[methyl_col] = p_value

    # Output p-values for each methylation column
    print("\nP-values for each methylation pattern:")
    for methyl_col, p_value in p_values.items():
        print(f'{methyl_col}: p-value = {p_value}')

    # Step 3: Interpretation and Thresholding
    significant_pmp = {key: value for key, value in p_values.items() if value < 0.05}
    print("\nSignificant PMPs with p-value < 0.05:")
    for methyl_col, p_value in significant_pmp.items():
        print(f'{methyl_col}: p-value = {p_value}')



