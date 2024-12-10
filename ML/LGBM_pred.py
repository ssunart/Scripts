###Random Forest와 LightGBM (LGBM) 회귀 모델을 사용하여 
###Molecular Descriptors를 기반으로 AutoDock-GPU 도킹 점수를 예측하는 머신러닝 작업
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
maccs_df = pd.read_csv('100k_maccsKey.csv')
maccs_df


#load_descriptor - RDKit key
#desc_df = pd.read_csv('100k_RDkit195.csv')
#desc_df


#merge_descriptor - ecfp,RDKit,maccs key
###dataframe 합치기
#merged_df = pd.concat([fp_df, maccs_df, desc_df], axis=1)
#merged_df.to_csv('merged100K_desc.csv', index=False)
#merged_df


##model 학습##

X = maccs_df #maccs_df #desc_df  #fp_df #new_mordred_df
y = score_df #AD_GPU Top score / best_score



# In[7]:


import sklearn.model_selection
import sklearn.ensemble
import optuna
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)


# In[ ]:


##Parameter optimization
#def objective(trial):
#    n_estimators = trial.suggest_int('n_estimators', 100, 3000)
#    max_depth = int(trial.suggest_float('max_depth', 1, 100, log=True))
#    learning_rate = float(trial.suggest_float(0.01, 0.1))
#    min_samples_split = trial.suggest_int('min_samples_split', 2, 15)
#    min_samples_leaf = trial.suggest_int('min_samples_leaf', 1, 10)
#    max_features = trial.suggest_categorical('max_features', ['auto', 'sqrt', 'log2'])
#    clf = RandomForestRegressor(
#        n_estimators=n_estimators,
#        max_depth=max_depth,
#        min_samples_split=min_samples_split,
#        min_samples_leaf=min_samples_leaf,
#        max_features=max_features
#    )
#    return cross_val_score(clf, X_train, y_train,
#                           n_jobs=-1, cv=3,
#                           scoring='neg_mean_squared_error').mean()
#study = optuna.create_study(direction='maximize')
#study.optimize(objective, n_trials=100)
#trial = study.best_trial
#print('MSE: {}'.format(trial.value))
#print("Best hyperparameters: {}".format(trial.params))



from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
# Optuna를 통해 얻은 최적의 하이퍼파라미터
##best_params = trial.params
# 모델 생성
#model = RandomForestRegressor(
#    n_estimators=best_params['n_estimators'],
#    max_depth=best_params['max_depth'],
#    min_samples_split=best_params['min_samples_split'],
#    min_samples_leaf=best_params['min_samples_leaf'],
#    max_features=best_params['max_features'],
#    n_jobs=2,
#    random_state=42)


# 모델 학습
#model.fit(X_train, y_train)
# 예측 및 성능 평가
#y_pred = model.predict(X_test)



##no optimization RF

#my_model = RandomForestRegressor(n_jobs=3, random_state=42)
#my_model.fit(X_train, y_train)

#y_pred = my_model.predict(X_test)


##index 중복
#cols=pd.Series(merged_df.columns)
#for dup in cols[cols.duplicated()].unique(): 
#    cols[cols[cols == dup].index.values.tolist()] = [dup + '_' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]
#merged_df.columns=cols

#X = maccs_df #maccs_df #desc_df  #fp_df #new_mordred_df
#y = score_df #AD_GPU Top score / best_score


# In[8]:


##LGBM model##
import sklearn.model_selection
import lightgbm as ltb
X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

import sklearn.ensemble
from sklearn.ensemble import RandomForestRegressor
my_model = ltb.LGBMRegressor(n_estimators=100, n_jobs=2, random_state=42)
my_model.fit(X_train, y_train)

y_pred = my_model.predict(X_test)


# In[ ]:


##결과보기##
from sklearn.metrics import mean_squared_error
mse1 = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error of this model is {mse1:.3f} (docking score unit)")

np.corrcoef(y_test, y_pred)

print(f"Pearson's correlation coefficient of our model is {np.corrcoef(y_test,y_pred)[0,1]:.3f}")

corr = np.corrcoef(y_test,y_pred)[0,1].round(3)
rmse = mean_squared_error(y_test, y_pred) ** 0.5
rmse = rmse.round(3)


# In[ ]:


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
plt.savefig("TotalRplot_maccs_LGBM.pdf")


# In[ ]:
