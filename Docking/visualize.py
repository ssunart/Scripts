###여러 도킹 점수(AutoDock-GPU, AutoDock-Vina, Glide SP, Glide XP, Fred, Dock6)에 대한 데이터를 분석
###각 도킹 방법의 피어슨 상관계수(Pearson Correlation) 분포를 히스토그램으로 시각화

#!/usr/bin/env python
# coding: utf-8

# In[52]:


import os, sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


data = pd.read_csv('soowon_docking.csv')
GPU = data['Autodock-GPU']
VINA = data['autodock-vina']
SP = data['GlideSP']
XP = data['GlideXP']
FRED = data['Fred']
DOCK = data['Dock6']

Gpu = GPU.describe()
Vina = VINA.describe()
Sp = SP.describe()
Xp = XP.describe()
Fred = FRED.describe()
Dock = DOCK.describe()


pro_describe = {'AutoDock-GPU' : Gpu, 'AutoDock-Vina' : Vina, 'Glide SP' : Sp, 'Glide XP' : Xp, 'Fred' : Fred, 'Dock6' : Dock}
pro_corr = {'AutoDock-GPU' : GPU, 'AutoDock-Vina' : VINA, 'Glide SP' : SP, 'Glide XP' : XP, 'Fred' : FRED, 'Dock6' : DOCK}
print(pro_describe)
print(pro_corr)
for pro, corr in pro_corr.items():
    avg = pro_describe[pro][1]
    avg = avg.round(2)
    
    std = pro_describe[pro][2]
    std = std.round(2)
    
    fir_q = pro_describe[pro][6]
    fir_q = fir_q.round(2)
    
    sec_q = pro_describe[pro][5]
    sec_q = sec_q.round(2)
    
    thr_q = pro_describe[pro][4]
    thr_q = thr_q.round(2)
    
    # avg = str(pro_describe[pro]).split()[3]
    # std = str(pro_describe[pro]).split()[5]
    # thr_q = str(pro_describe[pro]).split()[9]
    # sec_q = str(pro_describe[pro]).split()[11]
    # fir_q = str(pro_describe[pro]).split()[13]
    
    
    
    
    plt.axvline(0, 0, 100, color='lightgray', linestyle='--', linewidth=0.7,)
#plt.hist(corr, bins=22, histtype='bar')
    sns.histplot(corr, kde=True, color='royalblue', bins=22)
    plt.title(f'{pro}', fontsize=25)
    plt.xlim(-1, 1)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0, 20)
    plt.xlabel('Pearson Correlation', fontsize=20)
    plt.ylabel('Count', fontsize=20)
    plt.text(-0.96, 14.3,f"Avg = {avg}\nStd = {std}\nQ1  = {fir_q}\nQ2  = {sec_q}\nQ3  = {thr_q}", fontsize = 13)
    plt.savefig(f'{pro}_hist_quart_02.png')
    plt.close()


# In[ ]:
