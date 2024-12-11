###저장된 모델을 불러와 새로운 데이터 예측

import json
import lightgbm as lgb
import numpy as np
import pandas as pd

# 저장된 모델 불러오기
booster = lgb.Booster(model_file='../../lgbmDEF_MLM_dual_03.txt')

mordred_df=pd.read_csv('mordred_desc.csv')
maccs_df=pd.read_csv('maccsKey.csv')

feat_tot=pd.concat([mordred_df, maccs_df],axis=1)
cols=pd.Series(feat_tot.columns)
for dup in cols[cols.duplicated()].unique():
    cols[cols[cols == dup].index.values.tolist()] = [dup + '_' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]

feat_tot.columns=cols
feat_tot
X_test=feat_tot

train_01=pd.read_csv('../../mordred_desc.csv')
train_02=pd.read_csv('../../maccsKey.csv')

train_tot=pd.concat([train_01, train_02],axis=1)
cols=pd.Series(train_tot.columns)
for dup in cols[cols.duplicated()].unique():
    cols[cols[cols == dup].index.values.tolist()] = [dup + '_' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]

train_tot.columns=cols
train_tot
X_train=train_tot

# 훈련 데이터의 열 이름을 가져옵니다.
train_columns = X_train.columns.tolist()

# 테스트 데이터에서 훈련 데이터와 동일한 열만 선택합니다.
X_test_filtered = X_test[train_columns]

# 불러온 모델로 예측을 수행합니다.
y_pred = booster.predict(X_test_filtered)

# 예측값을 확인합니다.
print("예측 결과:", y_pred)

# 예측값을 pandas DataFrame으로 변환
df_pred = pd.DataFrame({'Prediction': y_pred})

# DataFrame을 CSV 파일로 저장
df_pred.to_csv('prediction_MLMresults_dual.csv', index=False)
