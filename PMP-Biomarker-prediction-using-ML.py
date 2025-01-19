import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
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

# Step 3: Define features (X) and target variable (y)
X = df[['strand'] + methylation_columns]  # Features: strand and methylation status columns
y = df['Tissue']  # Target: Tissue type


# Step 4: Split data into training and test sets (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Step 5: Train a Random Forest classifier
clf = RandomForestClassifier(random_state=42)
clf.fit(X_train, y_train)

# Step 6: Predict on test data
y_pred = clf.predict(X_test)

# Step 7: Evaluate the model
print("Classification Report:")
print(classification_report(y_test, y_pred))
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))

# Step 8: Cross-validation to assess model performance
cv_scores = cross_val_score(clf, X, y, cv=5, scoring='precision_macro')
print(f"Cross-validation precision scores: {cv_scores}")
print(f"Mean cross-validation precision: {cv_scores.mean()}")

# Step 9: Prediction Probabilities (Confidence Scores)
y_pred_prob = clf.predict_proba(X_test)
print("Prediction Probabilities for Test Set:")
print(y_pred_prob)
