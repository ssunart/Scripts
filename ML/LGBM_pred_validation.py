###
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
fp_df = pd.read_csv('total_ecfp4_1024.csv')
fp_df

#load_descriptor - MACCS key
maccs_df = pd.read_csv('total_maccsKey.csv')
maccs_df

#load_descriptor - RDKit key
desc_df = pd.read_csv('total_RDkit195.csv')
desc_df

#merge_descriptor - ecfp,RDKit,maccs key
###dataframe 합치기
merged_df = pd.concat([fp_df, maccs_df, desc_df], axis=1)
merged_df.to_csv('merged100K_desc.csv', index=False)
merged_df

##model 학습##

X = desc_df #maccs_df #desc_df  #fp_df #new_mordred_df
y = score_df #AD_GPU Top score / best_score

import sklearn.model_selection
import sklearn.ensemble
import optuna
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

def objective(trial):
    import lightgbm as ltb
    X = desc_df #maccs_df #desc_df  #fp_df #new_mordred_df
    y = score_df #AD_GPU Top score / best_score
    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)
    
    n_estimators = trial.suggest_int('n_estimators', 100, 3000)
    max_depth = int(trial.suggest_float('max_depth', 1, 100, log=True))
    learning_rate = float(trial.suggest_float('learning_rate', 0.01, 0.1))
    min_child_samples = trial.suggest_int('min_child_samples', 1, 20)
    clf = ltb.LGBMRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        min_child_samples=min_child_samples,
        n_jobs=3
    )
    return cross_val_score(clf, X_train, y_train,
                           n_jobs=3, cv=3,
                           scoring='neg_mean_squared_error').mean()
study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=100)
trial = study.best_trial
print('MSE: {}'.format(trial.value))
print("Best hyperparameters: {}".format(trial.params))

import sklearn.model_selection
import lightgbm as ltb
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
# Optuna를 통해 얻은 최적의 하이퍼파라미터
X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

best_params = trial.params
# 모델 생성
model = ltb.LGBMRegressor(
    n_estimators=best_params['n_estimators'],
    max_depth=26,
    learning_rate=best_params['learning_rate'],
    min_child_samples=best_params['min_child_samples'],
    n_jobs=3,
    random_state=42)
# 모델 학습 - LGBM
model.fit(X_train, y_train)
# 예측 및 성능 평가
y_pred = model.predict(X_test)

##no optimization RF

my_model = RandomForestRegressor(n_jobs=3, random_state=42)
my_model.fit(X_train, y_train)

y_pred = my_model.predict(X_test)

##index 중복
cols=pd.Series(merged_df.columns)
for dup in cols[cols.duplicated()].unique(): 
    cols[cols[cols == dup].index.values.tolist()] = [dup + '_' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]
merged_df.columns=cols

X = maccs_df #maccs_df #desc_df  #fp_df #new_mordred_df
y = score_df #AD_GPU Top score / best_score

##LGBM model##
import sklearn.model_selection
import lightgbm as ltb
X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

import sklearn.ensemble
from sklearn.ensemble import RandomForestRegressor
my_model = ltb.LGBMRegressor(n_estimators=100, n_jobs=2, random_state=42)
my_model.fit(X_train, y_train)

y_pred = my_model.predict(X_test)

##결과보기##
from sklearn.metrics import mean_squared_error
mse1 = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error of this model is {mse1:.3f} (docking score unit)")

np.corrcoef(y_test, y_pred)

print(f"Pearson's correlation coefficient of our model is {np.corrcoef(y_test,y_pred)[0,1]:.3f}")

corr = np.corrcoef(y_test,y_pred)[0,1].round(3)
rmse = mean_squared_error(y_test, y_pred) ** 0.5
rmse = rmse.round(3)

plt.scatter(y_test, y_pred, s=3, alpha=0.2)
plt.xlabel("AutoDock-GPU docking score", fontsize='xx-large')
plt.ylabel("Predicted docking score", fontsize='xx-large')
plt.grid()
plt.plot(range(-15, 0), range(-15, 0), "r--", linewidth=0.8, )
plt.xlim(-15, 0)
plt.ylim(-15, 0)
plt.text(-15.4, -1.1, f"Pearson's corr:{corr}", fontsize=12)
plt.text(-15.4, -2.0, f"RMSE:{rmse}", fontsize=12)
plt.axis('square')
plt.savefig("TotalRplot_rdkit195_LGBM_optim.pdf")
