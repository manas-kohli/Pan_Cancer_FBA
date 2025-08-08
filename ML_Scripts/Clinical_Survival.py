import pandas as pd
import numpy as np
import matplotlib
from sklearn import set_config
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OrdinalEncoder

from sksurv.datasets import load_gbsg2
from sksurv.ensemble import RandomSurvivalForest
from sksurv.preprocessing import OneHotEncoder
from sklearn.model_selection import GridSearchCV
from sksurv.metrics import as_concordance_index_ipcw_scorer
from sklearn.model_selection import StratifiedShuffleSplit

#### Initial Common Preprocessing

matplotlib.use('TkAgg')
TCGA_clinical= pd.read_csv('clinical_input/TCGA_clinical_input.csv')

TCGA_clinical['OS_binary'] = TCGA_clinical['OS'].map({1: True, 0: False})
TCGA_clinical_clean_OS = TCGA_clinical.dropna(subset=['OS_binary', 'OS.time'])
OS_array = np.array([(event, time) for event, time in zip(TCGA_clinical_clean_OS['OS_binary'], TCGA_clinical_clean_OS['OS.time'])],
                    dtype=[('event', bool), ('time', float)])
TCGA_clinical_input = TCGA_clinical_clean_OS.iloc[:, 1:]
columns_to_drop = ['OS_binary', 'OS', 'OS.time']
TCGA_clinical_input = TCGA_clinical_input.drop(columns_to_drop, axis = 1)
TCGA_clinical_input = TCGA_clinical_input.dropna(axis = 1)
TCGA_clinical_input = TCGA_clinical_input.loc[:, (TCGA_clinical_input != 0).any(axis = 0)]

columns_to_encode = [0] + list(range(2, 7))  # Adjusted for your dataset
categorical_data = TCGA_clinical_input.iloc[:, columns_to_encode].astype("category")  # Convert to string
non_categorical_data = TCGA_clinical_input.drop(TCGA_clinical_input.columns[columns_to_encode], axis=1)  # Drop categorical columns
# Apply OneHotEncoder
ohe = OneHotEncoder()  # drop="first" to avoid redundancy
encoded_array = ohe.fit_transform(categorical_data)

# Convert encoded data to DataFrame
encoded_df = pd.DataFrame(encoded_array, columns=ohe.get_feature_names_out(TCGA_clinical_input.columns[columns_to_encode]))

# Concatenate the processed categorical data with the original DataFrame
TCGA_clinical_input = pd.concat([non_categorical_data, encoded_df], axis=1)

set_config(display="text")  # displays text representation of estimators
random_state = 20

# Step 2: Bin survival time for stratification
num_bins = 10  # Adjust based on data distribution
y_binned = pd.qcut(OS_array["time"], q=num_bins, labels=False, duplicates="drop")

# Step 3: Convert tissue type into categorical labels
tissue_labels = TCGA_clinical["type"].astype("category").cat.codes

# Step 4: Create combined stratification label (time_bin + tissue_type)
stratify_labels = tissue_labels.astype(str)

# Define Stratified Shuffle Split
cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=random_state)

from sksurv.metrics import concordance_index_ipcw

### Define Scorers here

def ipcw_scorer(estimator, X, y):
    """Custom scorer using IPCW C-index."""
    surv_probs = estimator.predict(X)
    c_index_ipcw = concordance_index_ipcw(y, y, surv_probs)[0]  # IPCW-adjusted C-index
    return c_index_ipcw

### Random Survival Forest

param_grid = {
    "n_estimators": [100, 300, 500],  # Number of trees
    "min_samples_split": [5, 10, 20],  # Minimum samples to split
    "max_depth": [None, 10, 20],  # Maximum depth of trees
}

rsf = RandomSurvivalForest(random_state=random_state)

grid_search = GridSearchCV(
    estimator=rsf,
    param_grid=param_grid,
    scoring=ipcw_scorer,  # Use IPCW scorer
    cv=cv.split(TCGA_clinical_input, tissue_labels),  # Use Stratified Shuffle Split
    n_jobs=-1  # Use all available cores
)

grid_search.fit(TCGA_clinical_input, OS_array)


#SAVE THE PREVIOUS MODEL AND THE PREVIOUS CV RESULTS BEFORE YOU DO ANYTHING ELSE!!!!
cv_results_rf = pd.DataFrame(grid_search.cv_results_)
cv_results_rf.to_csv('clinical_input/cv_results_clinical_rxns_RF.csv', index = True)

best_model_rf = grid_search.best_estimator_
import pickle

model_pkl_file = "clinical_input/RF_clinical_rxns.pkl"

with open(model_pkl_file, 'wb') as file:
    pickle.dump(best_model_rf, file)

from sklearn.inspection import permutation_importance

result = permutation_importance(best_model_rf, TCGA_clinical_input, OS_array, n_repeats=3, random_state=random_state)
fimportance = pd.DataFrame(
    {
        k: result[k]
        for k in (
            "importances_mean",
            "importances_std",
        )
    },
    index=TCGA_clinical_input.columns,
).sort_values(by="importances_mean", ascending=False)

fimportance.to_csv('clinical_input/feature_importance_clinical_RF.csv', index = True)

### Cox Elastic Net Regression

from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis
import warnings
from sklearn.exceptions import FitFailedWarning
from sklearn.model_selection import GridSearchCV, KFold

coxnet_pipe = CoxnetSurvivalAnalysis(l1_ratio=0.9, max_iter=100)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", FitFailedWarning)
coxnet_pipe.fit(TCGA_clinical_input, OS_array)

estimated_alphas = coxnet_pipe.alphas_
gcv = GridSearchCV(
    CoxnetSurvivalAnalysis(l1_ratio=0.9, max_iter=100),
    param_grid={"alphas": [[v] for v in estimated_alphas]},
    scoring=ipcw_scorer,
    cv=cv.split(TCGA_clinical_input, tissue_labels),
    error_score=0.5,
    n_jobs=-1,
).fit(TCGA_clinical_input, OS_array)

cv_results_cox = pd.DataFrame(gcv.cv_results_)

best_model_Cox = gcv.best_estimator_
best_coefs_Cox = pd.DataFrame(best_model_Cox.coef_, index=TCGA_clinical_input.columns, columns=["coefficient"])

non_zero = np.sum(best_coefs_Cox.iloc[:, 0] != 0)
print(f"Number of non-zero coefficients: {non_zero}")

non_zero_coefs = best_coefs_Cox.query("coefficient != 0")
coef_order = non_zero_coefs.abs().sort_values("coefficient").index

cv_results_cox.to_csv('clinical_input/cv_results_clinical_rxns_COX.csv', index = True)
import pickle
model_pkl_file = "clinical_input/COX_clinical_rxns.pkl"
with open(model_pkl_file, 'wb') as file:
    pickle.dump(best_model_Cox, file)

best_coefs_Cox.to_csv('clinical_input/feature_importance_clinical_rxns_COX.csv', index = True)

### Gradient Boosting

from sksurv.ensemble import ComponentwiseGradientBoostingSurvivalAnalysis, GradientBoostingSurvivalAnalysis

param_grid_gbsa = {
    "learning_rate": [0.1, 0.25, 0.5, 0.9, 1.0],  # Step size for boosting
    "max_depth": [2, 3, 5, 7],
    "n_estimators": [100, 300, 500] # Maximum tree depth
}

gbsa = GradientBoostingSurvivalAnalysis(random_state=random_state)

gcv = GridSearchCV(
    estimator=gbsa,
    param_grid=param_grid_gbsa,
    scoring=ipcw_scorer,  # Use IPCW-adjusted C-index
    cv=cv.split(TCGA_clinical_input, tissue_labels),
    n_jobs=-1,  # Use all cores
    verbose=1  # Show progress
)

gcv.fit(TCGA_clinical_input, OS_array)

cv_results_gbsa = pd.DataFrame(gcv.cv_results_)

best_model_gbsa = gcv.best_estimator_
best_coefs_gbsa = pd.DataFrame(best_model_gbsa.feature_importances_, best_model_gbsa.feature_names_in_)

cv_results_gbsa.to_csv('clinical_input/cv_results_clinical_rxns_GBSA.csv', index = True)
import pickle
model_pkl_file = "clinical_input/GBSA_clinical_rxns.pkl"
with open(model_pkl_file, 'wb') as file:
    pickle.dump(best_model_gbsa, file)

best_coefs_gbsa.to_csv('clinical_input/feature_importance_clinical_rxns_GBSA.csv', index = True)

### Support Vector Machine
from sksurv.svm import FastSurvivalSVM

svm = FastSurvivalSVM(max_iter=1000, tol=1e-5, random_state=0)
param_grid_svm = {"alpha": 2.0 ** np.arange(-12, 13, 2)}
gcv = GridSearchCV(
    estimator=svm,
    param_grid=param_grid_svm,
    scoring=ipcw_scorer,  # Use IPCW-adjusted C-index
    cv=cv.split(TCGA_clinical_input, tissue_labels),
    n_jobs=-1,  # Use all cores
    verbose=1  # Show progress
)

gcv.fit(TCGA_clinical_input, OS_array)
cv_results_svm = pd.DataFrame(gcv.cv_results_)

best_model_svm = gcv.best_estimator_
best_coefs_svm = pd.DataFrame(best_model_svm.coef_, best_model_svm.feature_names_in_)

cv_results_svm.to_csv('clinical_input/cv_results_clinical_rxns_SVM.csv', index = True)
import pickle
model_pkl_file = "clinical_input/SVM_clinical_rxns.pkl"
with open(model_pkl_file, 'wb') as file:
    pickle.dump(best_model_svm, file)

best_coefs_svm.to_csv('clinical_input/feature_importance_clinical_rxns_SVM.csv', index = True)
### Kernel Support Vector Machine
from sksurv.kernels import clinical_kernel
from sksurv.svm import FastKernelSurvivalSVM

kernel_matrix = clinical_kernel(TCGA_clinical_input)
kssvm = FastKernelSurvivalSVM(optimizer="rbtree", kernel="precomputed", random_state=0)
gcv = GridSearchCV(
    estimator=kssvm,
    param_grid=param_grid_svm,
    scoring=ipcw_scorer,  # Use IPCW-adjusted C-index
    cv=cv.split(TCGA_clinical_input, tissue_labels),
    refit=True,
    n_jobs=-1,  # Use all cores
    verbose=1  # Show progress
)

gcv.fit(kernel_matrix, OS_array)

cv_results_ksvm = pd.DataFrame(gcv.cv_results_)

best_model_ksvm = gcv.best_estimator_
result = permutation_importance(best_model_ksvm, TCGA_clinical_input, OS_array, n_repeats=3, random_state=random_state)
fimportance = pd.DataFrame(
    {
        k: result[k]
        for k in (
            "importances_mean",
            "importances_std",
        )
    },
    index=TCGA_clinical_input.columns,
).sort_values(by="importances_mean", ascending=False)

fimportance.to_csv('clinical_input/feature_importance_clinical_indivrxns_KSVM.csv', index = True)

cv_results_ksvm.to_csv('clinical_input/cv_results_clinical_indiv_rxns_KSVM.csv', index = True)
import pickle
model_pkl_file = "clinical_input/KSVM_clinical_indiv_rxns.pkl"
with open(model_pkl_file, 'wb') as file:
    pickle.dump(best_model_ksvm, file)

