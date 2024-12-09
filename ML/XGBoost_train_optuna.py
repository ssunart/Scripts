### XGBoost regressor training & hyperparameter optimization with Optuna

#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
from tqdm import tqdm

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sklearn
from tqdm import tqdm # to show a progress bar.


#load_descriptor - AD_GPU
df_score = pd.read_csv('total_scoreDF.csv')
df_score


score_df=df_score["AD_GPU"]



#load_descriptor - ECFP
#fp_df = pd.read_csv('total_ecfp4_1024.csv')
#fp_df


#load_descriptor - MACCS key
#maccs_df = pd.read_csv('total_maccsKey.csv')
#maccs_df


#load_descriptor - RDKit key
#desc_df = pd.read_csv('total_RDkit195.csv')
#desc_df


#merge_descriptor - ecfp,RDKit,maccs key
###dataframe 합치기
merged_df = pd.read_csv('mergedTotal_desc.csv')


##model 학습##

X = merged_df #maccs_df #desc_df  #fp_df #new_mordred_df
y = score_df #AD_GPU Top score / best_score



# In[7]:


import sklearn.model_selection
import sklearn.ensemble
import optuna
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_regression
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_squared_error


X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)


# In[ ]:


##Parameter optimization
def objective(trial):
    n_estimators = trial.suggest_int('n_estimators', 10, 3000)
    learning_rate = trial.suggest_float('learning_rate', 0.0001, 0.05)
    max_depth = trial.suggest_int('max_depth', 1, 100)
    max_leaves = trial.suggest_int('max_leaves', 0, 100) 
    min_child_weight = trial.suggest_int('min_child_weight', 1, 100)

    clf = xgb.XGBRegressor(
        n_estimators=n_estimators,
        learning_rate=learning_rate,
        max_depth=max_depth,
        max_leaves=max_leaves,
        min_child_weight=min_child_weight,
        tree_method='gpu_hist',
        gpu_id=0,
        n_jobs=3
    )
    clf.fit(X_train, y_train)
    y_pred=clf.predict(X_test)
    return mean_squared_error(y_test, y_pred)

study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=100)
trial = study.best_trial
print('MSE: {}'.format(trial.value))
print("Best hyperparameters: {}".format(trial.params))

x = optuna.visualization.matplotlib.plot_param_importances(study)
x.figure.savefig('param_import.pdf')

x = optuna.visualization.matplotlib.plot_optimization_history(study)
x.figure.savefig('opt_his.pdf')

x = optuna.visualization.matplotlib.plot_slice(study)
x.figure.savefig('slice.pdf')
