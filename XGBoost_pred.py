#XGBoost Regressor)을 사용하여 분자 데이터를 기반으로 예측을 수행하고, 결과를 처리 및 저장하는 전체 워크플로

#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Module load

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


# In[2]:


#model & descriptors load
import xgboost as xgb
merged_df=pd.read_csv('scr_lib_merged_desc.csv')
# Sample 10% of the rows from the dataset
#sample_df = merged_df.sample(frac=0.1)

# Make predictions for test data

loaded_model = xgb.XGBRegressor()
loaded_model.load_model('DSpred_xgbOpt.json')


# In[3]:


merged_df


# In[ ]:


# Sample 10% of the rows from the dataset
#sample_df = merged_df.sample(frac=0.5)
#sample_df


# In[4]:


X_test = merged_df #sample_df
# make predictions for test data
predictions = loaded_model.predict(X_test)


# In[5]:


predictions


# In[6]:


df_pred = pd.DataFrame(predictions, columns=['prediction'])
df_pred.hist()


# In[7]:


df_pred.describe()


# In[8]:


##중복제거 및 별도 csv file 저장
df_molecules=pd.read_csv('scr_lib_pred.csv')
df_id=df_molecules['CatalogID']
df_smiles=df_molecules['SMILES']


# In[9]:


df_molecules


# In[10]:


mergePred_df = pd.concat([df_id, df_smiles, df_pred], axis=1)
mergePred_df


# In[11]:


# Sort DataFrame by the absolute value of the prediction column
mergePred_df['abs_prediction'] = mergePred_df['prediction'].abs()
mergePred_df.sort_values(by='abs_prediction', ascending=False, inplace=True)

# Add rank column based on the absolute prediction value
mergePred_df['rank'] = mergePred_df['abs_prediction'].rank(method='min', ascending=False)

# Remove the temporary absolute value column
mergePred_df.drop(columns='abs_prediction', inplace=True)

# Check the final DataFrame
print(mergePred_df)


# In[12]:


mergePred_df.to_csv('scr_lib_pred.csv', index=False)


# In[ ]:

