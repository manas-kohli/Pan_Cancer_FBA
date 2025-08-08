from pandas.core.common import random_state
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sksurv.metrics import concordance_index_ipcw
from sklearn.model_selection import GridSearchCV
from sksurv.ensemble import RandomSurvivalForest
import pandas as pd
import numpy as np


TCGA_flux= pd.read_csv('flux_input/TCGA_flux_combined_input.csv')

TCGA_flux['OS_binary'] = TCGA_flux['OS'].map({1: True, 0: False})
TCGA_flux_clean_OS = TCGA_flux.dropna(subset=['OS_binary', 'OS.time'])
OS_array = np.array([(event, time) for event, time in zip(TCGA_flux_clean_OS['OS_binary'], TCGA_flux_clean_OS['OS.time'])],
                            dtype=[('event', bool), ('time', float)])
TCGA_flux_input = TCGA_flux_clean_OS.iloc[:, 10:]
columns_to_drop = ['OS_binary']
TCGA_flux_input = TCGA_flux_input.drop(columns_to_drop, axis = 1)
TCGA_flux_input = TCGA_flux_input.dropna(axis = 1)
TCGA_flux_input = TCGA_flux_input.loc[:, (TCGA_flux_input != 0).any(axis = 0)]

random_state = 20
# Step 1: Bin survival time for stratification
num_bins = 10
y_binned = pd.qcut(OS_array["time"], q=num_bins, labels=False, duplicates="drop")

# Step 2: Convert tissue type into categorical labels
tissue_labels = TCGA_flux["type"].astype("category").cat.codes

# Step 3: Create combined stratification label (time_bin + tissue_type)
stratify_labels = tissue_labels.astype(str) + "_" + y_binned.astype(str)
# Count group sizes
group_counts = stratify_labels.value_counts()

# Filter out groups with fewer than 2 samples
valid_groups = group_counts[group_counts >= 2].index
filtered_idx = stratify_labels.isin(valid_groups)

# Filter Data and Labels
X_filtered = TCGA_flux_input[filtered_idx]
y_filtered = OS_array[filtered_idx]
stratify_filtered = stratify_labels[filtered_idx]

# ✅ Step 4: Train/Test Split with Stratification
X_train, X_test, y_train, y_test = train_test_split(
    X_filtered,
    y_filtered,
    stratify=stratify_filtered,
    test_size=0.2,
    random_state=random_state,
)

strat_train = stratify_filtered.loc[X_train.index]

# ✅ Step 5: Stratified Shuffle Split Matching Training Data
cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=random_state)

# IPCW Scorer
def ipcw_scorer(estimator, X, y):
    """Custom scorer using IPCW C-index."""
    surv_probs = estimator.predict(X)
    c_index_ipcw = concordance_index_ipcw(y, y, surv_probs)[0]
    return c_index_ipcw

# ✅ Step 6: Random Survival Forest with Grid Search
param_grid = {
    "n_estimators": [100, 300, 500],
    "min_samples_split": [5, 10, 20],
    "max_depth": [None, 10, 20],
}

rsf = RandomSurvivalForest(random_state=random_state)

grid_search = GridSearchCV(
    estimator=rsf,
    param_grid=param_grid,
    scoring=ipcw_scorer,
    cv=cv.split(X_train, strat_train),  # Use stratified training folds
    n_jobs=-1,
)

# ✅ Step 7: Fit the Model
grid_search.fit(X_train, y_train)
cv_results_rf_train = pd.DataFrame(grid_search.cv_results_)

# ✅ Step 8: Evaluate on Test Data
best_model = grid_search.best_estimator_

# ✅ Step 9: Predict Risk Scores for Test Data
risk_scores = grid_search.best_estimator_.predict(X_test)

# ✅ Step 10: Create DataFrame of Risk Scores
risk_scores_df = pd.DataFrame({
    "sample_index": X_test.index,
    "risk_score": risk_scores
}).set_index("sample_index")

# ✅ Step 11: Merge with Original TCGA Data
risk_output = TCGA_flux.loc[X_test.index].copy()
risk_output["risk_score"] = risk_scores_df["risk_score"]

# ✅ Step 12: Save Risk Score Data
risk_output.to_csv("flux_input/TCGA_flux_risk_scores.csv")