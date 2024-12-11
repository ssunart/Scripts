###여러 feature를 이용한 물질특성 예측모델 제작 플로우
import pandas as pd
raw_data=pd.read_csv('cleaned_results_Total.csv')

feat_01=pd.read_csv('maccsKey.csv')

###mol2vec 해당###
# mol2vec-000부터의 열만 선택
selected_columns = [col for col in feat_01.columns if col.startswith("mol2vec-")]

# 선택된 열만으로 새로운 DataFrame 생성
feat_01_selected = feat_01[selected_columns]
feat_01 = feat_01[selected_columns]

feat_02=pd.read_csv('mordred_desc.csv')

feat_tot=pd.concat([feat_01, feat_02],axis=1)
cols=pd.Series(feat_tot.columns)
for dup in cols[cols.duplicated()].unique():
    cols[cols[cols == dup].index.values.tolist()] = [dup + '_' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]

feat_tot.columns=cols
feat_tot

target=raw_data['MLM']
data=pd.concat([feat_tot, target],axis=1)

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

##model 학습##

X = feat_tot #maccs_df #desc_df  #fp_df #merged_df
y = target #MLM HLM

import sklearn.model_selection
import sklearn.ensemble

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)
y_train = y_train.values.ravel()
y_test = y_test.values.ravel()

import sklearn.model_selection
import lightgbm as ltb
X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

y_train = y_train.values.ravel()
y_test = y_test.values.ravel()

import sklearn.ensemble
from sklearn.ensemble import RandomForestRegressor
my_model = ltb.LGBMRegressor(device_type='cpu', random_state=42)
my_model.fit(X_train, y_train)

y_pred = my_model.predict(X_test)
booster = model.booster_
booster.save_model('lgbmDEF_MLM_dual_03.txt')

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
plt.xlabel("HLM", fontsize='xx-large')
plt.ylabel("Predicted HLM", fontsize='xx-large')
plt.grid()
plt.plot(range(0, 100), range(0, 100), "r--", linewidth=0.8, )
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.text(0, 95, f"Pearson's corr:{corr}", fontsize=12)
plt.text(0, 85, f"RMSE:{rmse}", fontsize=12)
plt.axis('square')
plt.savefig("LGBM__MLM_dual_03.pdf")

###최적화
import optuna

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

y_train = y_train.values.ravel()
y_test = y_test.values.ravel()

import lightgbm as ltb

def objective(trial):

    X = feat_tot #maccs_df #desc_df  #fp_df #new_mordred_df
    y = target #MLM

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

    y_train = y_train.values.ravel()
    y_test = y_test.values.ravel()

    n_estimators = trial.suggest_int('n_estimators', 100, 3000)
    #max_depth = int(trial.suggest_float('max_depth', 1, 100, log=True))
    max_depth = int(trial.suggest_float('max_depth', 1, 32, log=True))
    #num_leaves = trial.suggest_int('num_leaves', 2, 2**max_depth)
    learning_rate = float(trial.suggest_float('learning_rate', 0.01, 0.1))
    min_child_samples = trial.suggest_int('min_child_samples', 1, 20)
    clf = ltb.LGBMRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        #num_leaves=num_leaves,
        min_child_samples=min_child_samples,
        device_type='cpu',
        n_jobs=3, verbose=-1
    )
    return cross_val_score(clf, X_train, y_train,
                           cv=3,
                           scoring='neg_mean_squared_error').mean()
study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=100) #colab test를 위해 'n_trials=3'로 설정
#n_trials: parameter optimization trial 수 조정 ## 실제 연구시 'n_trials=100' 으로 설정하였습니다.
trial = study.best_trial
print('MSE: {}'.format(trial.value))
print("Best hyperparameters: {}".format(trial.params))

x=optuna.visualization.matplotlib.plot_param_importances(study)
x.figure.savefig('dual_03param_import_MLM.pdf')

x = optuna.visualization.matplotlib.plot_optimization_history(study)
x.figure.savefig('dual_03opt_his_MLM.pdf')
x = optuna.visualization.matplotlib.plot_slice(study)
x[0].figure.savefig('dual_03slice_MLM.pdf')

##최적화된 파라미터 적용후 학습
import sklearn.model_selection
import lightgbm as ltb

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
# Optuna를 통해 얻은 최적의 하이퍼파라미터
best_params = trial.params
# 모델 생성
model = ltb.LGBMRegressor(
    n_estimators=best_params['n_estimators'],
    max_depth=int(best_params['max_depth']),
    learning_rate=best_params['learning_rate'],
    min_child_samples=best_params['min_child_samples'],
    device_type='cpu',
    n_jobs=3,verbose=-1,
    random_state=42)

# 모델 학습
X = feat_tot #maccs_df #desc_df  #fp_df #new_mordred_df
y = target #AD_GPU Top score / best_score

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

y_train = y_train.values.ravel()
y_test = y_test.values.ravel()

model.fit(X_train, y_train)
# 예측 및 성능 평가
y_pred = model.predict(X_test)

booster = model.booster_
booster.save_model('lgbmOPT_MLM_dual_03.txt')

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
plt.xlabel("HLM", fontsize='xx-large')
plt.ylabel("Predicted MLM", fontsize='xx-large')
plt.grid()
plt.plot(range(0, 100), range(0, 100), "r--", linewidth=0.8, )
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.text(0, 95, f"Pearson's corr:{corr}", fontsize=12)
plt.text(0, 85, f"RMSE:{rmse}", fontsize=12)
plt.axis('square')
plt.savefig("LGBMopt__MLM_dual_03.pdf")
