###최적의 ML 알고리즘 동시다발 선별하기

raw_data=pd.read_csv('desc_0919/cleaned_results_Total.csv')
feat_01=pd.read_csv('desc_0919/ecfp2_2048.csv')
target=raw_data['HLM']
data=pd.concat([feat_01, target],axis=1)

from pycaret.regression import *
s = setup(data, target = 'HLM', experiment_name = 'MS_18', session_id =123, train_size = 0.8, use_gpu = True)
best = compare_models(n_select =5)
