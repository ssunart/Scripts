###PPAR과 ERRa 간 결합 친화도의 분포 및 상관관계 시각화
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 임시 데이터 생성
#np.random.seed(0)
#data = np.random.rand(100, 2)
#df = pd.DataFrame(data, columns=['Feature 1', 'Feature 2'])
merged_df
# 그림(figure) 설정
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# 첫 번째 서브플롯에 히스토그램 표시
ax[0].hist(merged_df['ensemble_docking_PPAR'], bins=20, color='skyblue', edgecolor='black')
ax[0].set_title('AKScore2_PPAR')
ax[0].set_xlabel('Value')
ax[0].set_ylabel('Frequency')

# 두 번째 서브플롯에 산점도 표시
ax[1].scatter(merged_df['ensemble_docking_PPAR'], merged_df['ensemble_docking_ERRa'], alpha=0.7)
ax[1].set_title('AKScore2 PPAR vs AD_GPU_ERRa')
ax[1].set_xlabel('PPAR')
ax[1].set_ylabel('ERRa')

# 서브플롯 간격 조정
plt.tight_layout()

# 표시
plt.show()

import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import scipy.stats

# 데이터프레임을 가정하고, 여기에는 'column_x'와 'column_y'라는 두 열이 있다고 합시다.
# merged_df = pd.DataFrame({
#     'column_x': [데이터를 채워주세요],
#     'column_y': [데이터를 채워주세요]
# })

# 조인트플롯 생성
jointplot = sns.jointplot(data=merged_df, x='ensemble_docking_PPAR', y='ensemble_docking_ERRa', kind="scatter",marginal_kws={'bins': 20}, joint_kws={'s': 10} )

# 상관 계수 및 p-값 계산 후 표시
r, p = scipy.stats.pearsonr(merged_df['ensemble_docking_PPAR'], merged_df['ensemble_docking_ERRa'])
#jointplot.ax_joint.text(0.05, 0.95, f'pearsonr = {r:.2f}; p = {p:.2e}', 
#                        transform=jointplot.ax_joint.transAxes, 
#                        ha='left', 
#                        va='top')
# x와 y축의 범위를 동일하게 설정
ax = jointplot.ax_joint
xlim = ax.get_xlim()
ylim = ax.get_ylim()

# 더 넓은 범위를 찾아서 두 축에 적용
lim = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
ax.set_xlim(lim)
ax.set_ylim(lim)
# x=y 선 추가
plt.plot([-18, -6], [-18, -6], 'k--', color='red')  # 점선으로 표시
# 그래프 표시
plt.show()
