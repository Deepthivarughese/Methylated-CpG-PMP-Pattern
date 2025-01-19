import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import cross_val_score

# Step 1: Import the CSV file
# Replace 'your_data.csv' with the path to your CSV file
df = pd.read_csv(r'F:\DEEPTHI\SELF STUDY\sample-data.csv')
df.columns = df.columns.str.strip()  # Removes extra spaces from column names

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


# Step 4: Logistic Regression (to assess the odds of tissue classification based on methylation patterns)
X = df[['strand'] + methylation_columns]  # Features: strand and methylation status columns
y = df['Tissue']  # Target: Tissue type

# Split data into training and testing sets (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Fit a logistic regression model
logreg = LogisticRegression(max_iter=1000)
logreg.fit(X_train, y_train)


# Get the p-values for each feature using the model's coefficients
import statsmodels.api as sm
X_train_sm = sm.add_constant(X_train)  # Add constant term for intercept
logreg_sm = sm.Logit(y_train, X_train_sm)
result = logreg_sm.fit()

# Print p-values for each feature
print("\nP-values for logistic regression coefficients:")
print(result.summary())

# Step 5: Predict on test data and evaluate model
y_pred = logreg.predict(X_test)

# Print classification report and confusion matrix
print("\nClassification Report:")
print(classification_report(y_test, y_pred))

print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))