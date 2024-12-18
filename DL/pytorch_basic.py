#pytorch load & version check
import torch
print(torch.__version__)

import nummpy as np #numpy 라이브러리 불러오기
t= np.array([0., 1., 2., 3., 4., 5.])
print(t) # 1차원 벡터가 출력된다

# 1차원 텐서인 벡터의 차원과 크기 출력하기
print('Rank of t:', t.ndim) # 텐서의 차원 표시
print('Shape of t:', t.shape) # 텐서의 크기 표시 (6,0) 은 1x6 의 크기를 가지는 벡터를 의미한다.

#numpy는 python을 이용한 수치계산을 위해 사용되는 array이다.
##numpy에서 각 벡터의 원소에 접근하기 : 인덱싱 활용
print('t[0] t[1] t[-1] =', t[0], t[1], t[-1])

##슬라이싱을 통한 접근
print('t[2:4] t[3:-1] =', t[2:4], t[3:-1]) #마지막 범위는 접근에 포함되지 않는다
## 처음과 끝 번호 생략한 슬라이싱
print("t[:2] t[3:] =", t[:2], t[3:])

#2D array with numpy : Numpy로 2차원 행렬 만들기
t2 = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
print(t2)
t_dim = t2.ndim
t_shape = t2.shape
print(f"Rank of t2= {t_dim}, Shape of t2= {t_shape}")

t3 = np.array([[1],[2],[3]])
print(t3)
t3_dim = t3.ndim
t3_shape= t3.shape
print(f"Rank of t3= {t3_dim}, Shape of t3= {t3_shape}")

import torch
#파이토치로 1차원 텐서 벡터 만들기
t = torch.FloatTensor([0.,1.,2.,3.,4.,5.,6.])
t_dim = t.dim() ## rank, 차원
t_shape = t.shape ## shape
t_size = t.size() ##shape
print(t)
print(f"t는 {t_dim} 차원 텐서로 shape는 {t_shape} 이고 size는 {t_size} 이다.")
