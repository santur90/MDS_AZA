1. # =============================================================================
# Machine Learning Pipeline: Multi-omics Data Integration and Classification
# Purpose: Integrate methylation, gene expression, and mutation data to build
#          a predictive machine learning model using Random Forest classifier
# =============================================================================

import mglearn
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as backend
from scipy import interp
from scipy.stats import pearsonr

from sklearn import preprocessing, tree, svm, naive_bayes, neighbors, metrics, linear_model
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score, ShuffleSplit, GridSearchCV, KFold, StratifiedKFold
from sklearn.metrics import (average_precision_score, precision_score, roc_curve, auc, 
                            f1_score, classification_report, roc_auc_score, 
                            mean_squared_error, r2_score, make_scorer, accuracy_score)
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.svm import SVR, LinearSVC, SVC
from sklearn.base import TransformerMixin, BaseEstimator
from sklearn.feature_selection import SelectKBest, f_classif, f_regression, SelectFpr, SelectFdr, SelectFwe, GenericUnivariateSelect
from sklearn.datasets import make_classification
from sklearn.tree import plot_tree
from sklearn.decomposition import PCA
import os

# =============================================================================
# Configuration
# =============================================================================

# Working directory
os.chdir("code/")

# File paths
PHENOTYPE_FILE = "pData.csv"               # Sample phenotype/class labels
METHYLATION_FILE = "DNAme_percentage.txt"  # DNA methylation data (features × samples)
EXPRESSION_FILE = "GE.txt"                 # Gene expression data (features × samples)
MUTATION_FILE = "Mutation_Cli_infor.txt"   # Mutation and clinical information
TRAIN_SAMPLES_FILE = "Train_samples.csv"   # Training sample indices
TEST_SAMPLES_FILE = "Test_samples.csv"     # Test sample indices

# Model hyperparameters (to be configured)
N_ESTIMATORS_RANGE = (50, 200, 10)        # Range for number of trees: (start, stop, step)
MAX_DEPTH_RANGE = (5, 20, 2)              # Range for tree depth
FEATURE_K_RANGE = (10, 100, 10)           # Range for number of selected features
CV_FOLDS = 5                               # Number of cross-validation folds
RANDOM_STATE = 43                          # For reproducibility
SCORING_METRIC = "accuracy"                # Evaluation metric

print("="*70)
print("MULTI-OMICS MACHINE LEARNING PIPELINE")
print("="*70)
print("="*70)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

print("\n[STEP 1] Loading multi-omics data...")
print("-"*70)

# Load phenotype data
print("Loading phenotype data...")
phenotype_data = pd.read_csv(PHENOTYPE_FILE, sep=r'\s+')
print(f"✓ Phenotype data shape: {phenotype_data.shape}")

# Load methylation data (features × samples)
print("Loading DNA methylation data...")
methylation_data = pd.read_csv(METHYLATION_FILE, sep=r'\s', engine='python')
print(f"✓ Methylation data shape: {methylation_data.shape}")

# Load gene expression data (features × samples)
print("Loading gene expression data...")
expression_data = pd.read_csv(EXPRESSION_FILE, sep=r'\s', engine='python')
print(f"✓ Expression data shape: {expression_data.shape}")

# Load mutation and clinical information
print("Loading mutation and clinical data...")
mutation_data = pd.read_csv(MUTATION_FILE, sep=r'\s', engine='python')
print(f"✓ Mutation/Clinical data shape: {mutation_data.shape}")

# Load training and test sample information
print("Loading sample assignments...")
train_samples = pd.read_csv(TRAIN_SAMPLES_FILE, sep=r'\s', engine='python')
test_samples = pd.read_csv(TEST_SAMPLES_FILE, sep=r'\s', engine='python')
print(f"✓ Training samples: {train_samples.shape[0]}")
print(f"✓ Test samples: {test_samples.shape[0]}")

print(f"✓ Test samples: {test_samples.shape[0]}")

# =============================================================================
# 2. PREPARE DATA MATRICES
# =============================================================================

print("\n[STEP 2] Preparing data matrices...")
print("-"*70)

# Extract training and test samples for methylation
print("Filtering methylation data...")
train_methylation = methylation_data.loc[:, train_samples.index]
test_methylation = methylation_data.loc[:, test_samples.index]
print(f"✓ Train methylation: {train_methylation.shape}")
print(f"✓ Test methylation: {test_methylation.shape}")

# Extract training and test samples for expression
print("Filtering expression data...")
train_expression = expression_data.loc[:, train_samples.index]
test_expression = expression_data.loc[:, test_samples.index]
print(f"✓ Train expression: {train_expression.shape}")
print(f"✓ Test expression: {test_expression.shape}")

# Process mutation/clinical data
print("Processing mutation/clinical data...")
mutation_processed = mutation_data.copy()
mutation_processed = mutation_processed.fillna(-1)
mutation_processed = mutation_processed.astype(int)
mutation_processed = mutation_processed.astype(str)
mutation_processed = mutation_processed.replace('-1', np.nan)

# Extract training and test samples for mutation
train_mutation = mutation_processed.loc[train_samples.index, :]
test_mutation = mutation_processed.loc[test_samples.index, :]
print(f"✓ Train mutation/clinical: {train_mutation.shape}")
print(f"✓ Test mutation/clinical: {test_mutation.shape}")

# =============================================================================
# 3. PREPARE TRAINING DATASET
# =============================================================================

print("\n[STEP 3] Preparing training dataset...")
print("-"*70)

# Transpose data (features × samples → samples × features)
train_methylation_T = train_methylation.T
train_expression_T = train_expression.T

# Merge all omics data
train_merged = train_methylation_T.join(train_expression_T, how='inner')
train_merged = train_merged.join(train_mutation, how='inner')
print(f"✓ Merged training data shape: {train_merged.shape}")

# Merge with phenotype information
train_final = train_merged.merge(
    phenotype_data, 
    left_index=True, 
    right_index=True, 
    how='inner',
    sort=True
)
print(f"✓ Final training data shape: {train_final.shape}")

# Separate features and target
X_train = train_final.drop(['Type'], axis=1)
y_train = train_final["Type"]
print(f"X_train shape: {X_train.shape}")
print(f"Training class distribution:")
print(y_train.value_counts())

# =============================================================================
# 4. PREPARE TEST DATASET
# =============================================================================

print("\n[STEP 4] Preparing test dataset...")
print("-"*70)

# Transpose data
test_methylation_T = test_methylation.T
test_expression_T = test_expression.T

# Merge all omics data
test_merged = test_methylation_T.join(test_expression_T, how='inner')
test_merged = test_merged.join(test_mutation, how='inner')
print(f"✓ Merged test data shape: {test_merged.shape}")

# Merge with phenotype information
test_final = test_merged.merge(
    phenotype_data, 
    left_index=True, 
    right_index=True, 
    how='inner',
    sort=True
)
print(f"✓ Final test data shape: {test_final.shape}")

# Separate features and target
X_test = test_final.drop(['Type'], axis=1)
y_test = test_final["Type"]
print(f"X_test shape: {X_test.shape}")
print(f"Test class distribution:")
print(y_test.value_counts())
print(y_test.value_counts())

# =============================================================================
# 5. BUILD MACHINE LEARNING PIPELINE
# =============================================================================

print("\n[STEP 5] Building ML pipeline with hyperparameter tuning...")
print("-"*70)

# Define feature selector
print("Configuring feature selection method...")
feature_selector = SelectKBest()

# Define hyperparameter search space
print("Defining hyperparameter search space...")
hyperparameters = {
    "clf__n_estimators": list(range(*N_ESTIMATORS_RANGE)),
    "clf__max_depth": list(range(*MAX_DEPTH_RANGE)),
    "clf__max_features": ['log2', 'sqrt'],
    "clf__bootstrap": [True, False],
    "clf__criterion": ['gini', 'entropy'],
    "feature__k": list(range(*FEATURE_K_RANGE)),
    "feature__score_func": [f_classif],
}

# Build pipeline: Scaling → Feature Selection → Classification
print("Constructing pipeline...")
pipeline = Pipeline(steps=[
    ('scale', preprocessing.StandardScaler()),      # Data normalization
    ('feature', feature_selector),                  # Feature selection
    ('clf', RandomForestClassifier(random_state=RANDOM_STATE))  # Classifier
])

# Configure GridSearchCV for hyperparameter optimization
print("Setting up GridSearchCV for hyperparameter tuning...")
grid_search = GridSearchCV(
    pipeline,
    hyperparameters,
    cv=StratifiedKFold(CV_FOLDS),
    scoring=SCORING_METRIC,
    error_score='raise',
    verbose=1,
    n_jobs=-1
)

# =============================================================================
# 6. TRAIN MODEL AND FIND BEST PARAMETERS
# =============================================================================

print("\n[STEP 6] Training model with hyperparameter tuning...")
print("-"*70)
print(f"Training on {X_train.shape[0]} samples with {X_train.shape[1]} features...")

grid_search.fit(X_train, y_train)

# Extract best model and parameters
best_model = grid_search.best_estimator_
best_params = grid_search.best_params_
feature_selector_best = best_model.named_steps['feature']
classifier_best = best_model.named_steps['clf']

print(f"\n✓ Best parameters found:")
for param, value in best_params.items():
    print(f"  {param}: {value}")
print(f"\nBest CV score ({SCORING_METRIC}): {grid_search.best_score_:.4f}")

# =============================================================================
# 7. FEATURE IMPORTANCE AND SELECTION
# =============================================================================

print("\n[STEP 7] Analyzing feature importance...")
print("-"*70)

# Get feature importance from the best classifier
feature_importances = classifier_best.feature_importances_

# Get selected feature indices
selected_features_mask = feature_selector_best.get_support()
selected_features = X_train.columns[selected_features_mask]

print(f"Number of selected features: {len(selected_features)}")
print(f"\nTop 20 important features:")
feature_importance_df = pd.DataFrame({
    'feature': X_train.columns,
    'importance': feature_importances
}).sort_values('importance', ascending=False)
print(feature_importance_df.head(20).to_string())

# Apply feature selection to training data
X_train_selected = feature_selector_best.transform(X_train)
X_train_selected_df = pd.DataFrame(X_train_selected, columns=selected_features, index=X_train.index)
print(f"\nTraining data after feature selection: {X_train_selected_df.shape}")

# Apply feature selection to test data
X_test_selected = feature_selector_best.transform(X_test)
X_test_selected_df = pd.DataFrame(X_test_selected, columns=selected_features, index=X_test.index)
print(f"Test data after feature selection: {X_test_selected_df.shape}")

# =============================================================================
# 8. MODEL EVALUATION ON TEST SET
# =============================================================================

print("\n[STEP 8] Evaluating model on test set...")
print("-"*70)

# Make predictions on test set
y_pred_proba = best_model.predict_proba(X_test)
y_pred = best_model.predict(X_test)

# Calculate ROC curve and AUC
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba[:, 1])
roc_auc = auc(fpr, tpr)

print(f"✓ ROC-AUC Score: {roc_auc:.4f}")
print(f"✓ Accuracy: {accuracy_score(y_test, y_pred):.4f}")
print(f"✓ F1-Score: {f1_score(y_test, y_pred, average='weighted'):.4f}")

print(f"\nClassification Report:")
print(classification_report(y_test, y_pred))

print("\n" + "="*70)
print("PIPELINE COMPLETED SUCCESSFULLY")
print("="*70)
