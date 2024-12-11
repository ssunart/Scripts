###두 가지 데이터 분석 및 시각화
##PPAR과 ERRa 리간드의 공통 분자를 찾고, 이들의 결합 친화도를 다양한 방식으로 분석하고 시각화
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# CSV 파일 읽기
ppar_df = pd.read_csv('mw500_0131P.csv')
erra_df = pd.read_csv('mw500_0131E.csv')


# 공통 분자 찾기
common_ligands = ppar_df[ppar_df['CID'].isin(erra_df['CID'])]

# 두 데이터프레임 병합
merged_df = pd.merge(common_ligands, erra_df, on='CID', suffixes=('_PPAR', '_ERRa'))

# 밀도 플롯 생성
plt.figure(figsize=(4, 3))
sns.kdeplot(data=merged_df, x='ensemble_docking_PPAR', y='ensemble_docking_ERRa', cmap="Reds", fill=True, bw_adjust=1.0)

# 하이라이트 영역
#highlight = merged_df[(merged_df['AD_GPU_PPAR'] > -14) & (merged_df['AD_GPU_ERRa'] > -14)]
#plt.scatter(highlight['AD_GPU_PPAR'], highlight['AD_GPU_ERRa'], color='blue')

plt.xlabel('PPAR AK2 Score')
plt.ylabel('ERRa AK2 Score')
plt.title('Density Plot of binding affinities for Common Ligands in PPAR and ERRa')
# x=y 선 추가
plt.plot([-17, -5], [-17, -5], 'k--', color='grey')

# x와 y의 범위 설정
plt.xlim(-17, -5)
plt.ylim(-17, -5)

plt.grid(False)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt

# 시각화
plt.scatter(merged_df['ensemble_docking_PPAR'], merged_df['ensemble_docking_ERRa'], s=5)  # 점 크기 조절
plt.xlabel('PPAR AK2 Score')
plt.ylabel('ERRa AK2 Score')
plt.title('Binding affinities for Common Ligands in PPAR and ERRa')

# x=y 선 추가
plt.plot([-16.0, -6], [-16, -6], 'k--', color='red')  # 점선으로 표시

# x와 y의 범위 설정
plt.xlim(-16, -6)
plt.ylim(-16, -6)

plt.grid(False)
plt.show()

merged_df.to_csv("plot_df.csv")
merged_df['Total_Similarity_PPAR'].hist()
merged_df['Molecular_Weight_PPAR'].hist()
merged_df['ensemble_docking_PPAR'].hist()
merged_df['ensemble_docking_ERRa'].hist()


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



# 히스토그램 생성
plt.figure(figsize=(10, 6))
plt.hist(merged_df['Similarity_SR11023_PPAR'], bins=20, color='skyblue', edgecolor='black', density=False)

# 제목과 레이블 추가
plt.title('Similarity_SR11023')
plt.xlabel('Similarity Score')
plt.ylabel('Density')

# 그리드 추가
#plt.grid(True)

# 표시
plt.show()

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



# 히스토그램 생성
plt.figure(figsize=(10, 6))
plt.hist(merged_df['Total_Similarity_ERRa'], bins=20, color='skyblue', edgecolor='black', density=True)

# 제목과 레이블 추가
plt.title('Total Similarity ERRa Distribution')
plt.xlabel('Similarity Score')
plt.ylabel('Density')

# 그리드 추가
#plt.grid(True)

# 표시
plt.show()

merged_df['Molecular_Weight_PPAR'].hist()

# 히스토그램 생성
plt.figure(figsize=(10, 6))
plt.hist(merged_df['Molecular_Weight_PPAR'], bins=20, color='skyblue', edgecolor='black', density=False)

# 제목과 레이블 추가
plt.title('MW distribution')
plt.xlabel('Molecular weight')
plt.ylabel('Density')

# 그리드 추가
#plt.grid(True)

# 표시
plt.show()
merged_df['ensemble_docking_PPAR'].hist()

# 히스토그램 생성
plt.figure(figsize=(10, 6))
plt.hist(merged_df['ensemble_docking_PPAR'], bins=30, color='skyblue', edgecolor='black', density=True)

# 제목과 레이블 추가
plt.title('AKScore2_PPAR distribution')
plt.xlabel('AKScore 2')
plt.ylabel('Density')

# 그리드 추가
#plt.grid(True)

# 표시
plt.show()

merged_df['ensemble_docking_PPAR'].describe()
merged_df['ensemble_docking_ERRa'].hist()

# 히스토그램 생성
plt.figure(figsize=(10, 6))
plt.hist(merged_df['ensemble_docking_ERRa'], bins=30, color='skyblue', edgecolor='black', density=True)

# 제목과 레이블 추가
plt.title('AKScore2_ERRa distribution')
plt.xlabel('AKScore 2')
plt.ylabel('Density')

# 그리드 추가
#plt.grid(True)

# 표시
plt.show()
