###Docking 프로그램의 결과 점수(z-score)를 계산하고, 각 프로그램의 결과를 비교 및 분석.
#순위 매기기 및 평균 순위를 기반으로 상위 결과를 도출.
#결과를 CSV 파일로 저장하고 시각적으로 확인.

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from scipy.stats import zscore

# csv 파일 읽기
data = pd.read_csv("docking_results_default.csv")

# z-score 정규화를 수행할 컬럼 선택
columns_to_normalize = ['AutoDock-GPU Top score', 'AutoDock-Vina Top score', 'LeDock Top score']

# 각 컬럼에 대해 z-score 정규화 수행
for column in columns_to_normalize:
    data[column + ' z-score'] = zscore(data[column].dropna())

# 결과를 새로운 csv 파일에 저장
data.to_csv("docking_results_z_scores.csv", index=False)


# In[3]:


import pandas as pd
import numpy as np

# csv 파일 읽기
data = pd.read_csv("docking_results_default.csv")

# z-score 정규화를 수행할 컬럼 선택
columns_to_normalize = ['AutoDock-GPU Top score', 'AutoDock-Vina Top score', 'LeDock Top score']

# 각 컬럼에 대해 z-score 정규화 수행
for column in columns_to_normalize:
    mean = data[column].mean()
    std_dev = data[column].std()
    data[column + ' z-score'] = (data[column] - mean) / std_dev

# 결과를 새로운 csv 파일에 저장
data.to_csv("docking_results_z_scores.csv", index=False)


# In[4]:


import pandas as pd
import numpy as np

# csv 파일 읽기
data = pd.read_csv("docking_results_default.csv")

# z-score 정규화를 수행할 컬럼 선택
columns_to_normalize = ['AutoDock-GPU Top score', 'AutoDock-Vina Top score', 'LeDock Top score']

# 각 컬럼에 대해 z-score 정규화 수행
for column in columns_to_normalize:
    numeric_data = pd.to_numeric(data[column], errors='coerce')
    mean = numeric_data.mean()
    std_dev = numeric_data.std()
    data[column + ' z-score'] = (numeric_data - mean) / std_dev

# 결과를 새로운 csv 파일에 저장
data.to_csv("docking_results_z_scores.csv", index=False)


# In[9]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# csv 파일 읽기
data = pd.read_csv("docking_results_default.csv")

# z-score 정규화를 수행할 컬럼 선택
columns_to_normalize = ['AutoDock-GPU Top score', 'AutoDock-Vina Top score', 'LeDock Top score']

# 각 컬럼에 대해 z-score 정규화 수행
for column in columns_to_normalize:
    numeric_data = pd.to_numeric(data[column], errors='coerce')
    mean = numeric_data.mean()
    std_dev = numeric_data.std()
    data[column + ' z-score'] = (numeric_data - mean) / std_dev

# 각 프로그램별로 z-score를 기준으로 순위 매기기
for column in columns_to_normalize:
    data[column + ' rank'] = data[column + ' z-score'].rank()

# 평균 순위 계산
data['Average Rank'] = data[[col + ' rank' for col in columns_to_normalize]].mean(axis=1)

# 평균 순위를 기준으로 데이터 정렬
sorted_data = data.sort_values(by='Average Rank')

# 결과 출력
print("평균 순위를 기준으로 정렬된 결과:")
print(sorted_data.head(10))  # 상위 10개 결과 출력


# 히스토그램 그리기
for column in columns_to_normalize:
    plt.figure()
    plt.hist(data[column + ' z-score'].dropna(), bins=50, alpha=0.75, label=column + ' z-score')
    plt.xlabel('z-score')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    plt.title(column + ' z-score Distribution')
    plt.grid(True)
    plt.show()

# 결과를 새로운 csv 파일에 저장
data.to_csv("docking_results_z_scores.csv", index=False)


# In[ ]:
