# =============================================================================
# Train-Test Dataset Split
# Purpose: Load methylation data, merge with phenotype information, and split
#          into training and test sets for machine learning model development
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interp
from scipy.stats import pearsonr

from sklearn import tree, svm, naive_bayes, neighbors, metrics, linear_model
from sklearn.ensemble import (BaggingClassifier, AdaBoostClassifier, 
                             RandomForestClassifier, GradientBoostingClassifier,
                             ExtraTreesClassifier)
from sklearn.model_selection import (train_test_split, cross_val_score, ShuffleSplit, 
                                    GridSearchCV, KFold, StratifiedKFold)
from sklearn.metrics import (average_precision_score, precision_score, roc_curve, 
                            auc, f1_score, classification_report, roc_auc_score)
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.svm import SVR, LinearSVC, SVC
from sklearn.base import TransformerMixin, BaseEstimator
from sklearn.feature_selection import (VarianceThreshold, SelectFromModel, SelectKBest,
                                      chi2, mutual_info_classif, RFECV, RFE)
from sklearn.linear_model import LogisticRegression, LassoCV, Ridge
from sklearn.datasets import make_classification
# =============================================================================
# Configuration
# =============================================================================

# File paths
RESPONSE_KEY_FILE = "responseKey.txt"      # Format: sample_id, Type/condition
METHYLATION_DATA_FILE = "data.txt"         # Features as rows, samples as columns

# Output file paths
TRAIN_SAMPLES_OUTPUT = "Train_samples.csv"
TEST_SAMPLES_OUTPUT = "Test_samples.csv"

# Train-test split parameters
TEST_SIZE = 0.3                             # 30% for test, 70% for training
RANDOM_STATE = 42                           # For reproducibility
STRATIFY = True                             # Stratified split to maintain class distribution
# =============================================================================
# Load and Prepare Data
# =============================================================================

print("Loading response key (phenotype data)...")
response_data = pd.read_csv(RESPONSE_KEY_FILE, sep=r'\s+')
print(f"Response key shape: {response_data.shape}")

print("\nLoading methylation data...")
methylation_data = pd.read_csv(METHYLATION_DATA_FILE, sep=r'\s')
print(f"Methylation data shape: {methylation_data.shape}")
print("Note: Features as rows, samples as columns")

# Transpose methylation data (features → columns, samples → rows)
methylation_data_T = methylation_data.T
print(f"Transposed data shape: {methylation_data_T.shape}")

# Merge methylation data with phenotype information
print("\nMerging methylation data with phenotype information...")
merged_data = methylation_data_T.merge(
    response_data, 
    left_index=True, 
    right_index=True, 
    how='inner'
)
print(f"Merged data shape: {merged_data.shape}")
print("\nFirst few rows of merged data:")
print(merged_data.head())

# =============================================================================
# Prepare Features and Target Variable
# =============================================================================

print("\n" + "="*60)
print("Preparing features and target variable")
print("="*60)

# Identify target column (phenotype/class label)
TARGET_COLUMN = "Type"

# Check if target column exists
if TARGET_COLUMN not in merged_data.columns:
    raise ValueError(f"Target column '{TARGET_COLUMN}' not found in data. "
                    f"Available columns: {list(merged_data.columns)}")

# Separate features (X) and target (y)
X = merged_data.drop([TARGET_COLUMN], axis=1)
y = merged_data[TARGET_COLUMN]

print(f"\nFeatures (X) shape: {X.shape}")
print(f"Target (y) shape: {y.shape}")
print(f"\nClass distribution:")
print(y.value_counts())
print(f"\nClass proportions:")
print(y.value_counts(normalize=True))

# =============================================================================
# Train-Test Split
# =============================================================================

print("\n" + "="*60)
print("Splitting data into training and test sets")
print("="*60)
print(f"Test size: {TEST_SIZE * 100}%")
print(f"Training size: {(1 - TEST_SIZE) * 100}%")

X_train, X_test, y_train, y_test = train_test_split(
    X, y, 
    test_size=TEST_SIZE, 
    stratify=y if STRATIFY else None,
    random_state=RANDOM_STATE
)

print(f"\nTraining set size: {X_train.shape[0]} samples")
print(f"Test set size: {X_test.shape[0]} samples")
print(f"\nTraining set class distribution:")
print(y_train.value_counts())
print(f"\nTest set class distribution:")
print(y_test.value_counts())

# =============================================================================
# Export Sample Information
# =============================================================================

print("\n" + "="*60)
print("Exporting sample assignments")
print("="*60)

# Get training sample information
training_samples = response_data.loc[X_train.index]
print(f"\nTraining samples: {training_samples.shape[0]}")
print(training_samples.head())

# Export training samples
print(f"\nExporting training samples to: {TRAIN_SAMPLES_OUTPUT}")
training_samples.to_csv(TRAIN_SAMPLES_OUTPUT, sep="\t")
print("✓ Training samples exported successfully")

# Get test sample information
test_samples = response_data.loc[X_test.index]
print(f"\nTest samples: {test_samples.shape[0]}")
print(test_samples.head())

# Export test samples
print(f"\nExporting test samples to: {TEST_SAMPLES_OUTPUT}")
test_samples.to_csv(TEST_SAMPLES_OUTPUT, sep="\t")
print("✓ Test samples exported successfully")

print("\n" + "="*60)
print("Data split and export completed!")
print("="*60) 
