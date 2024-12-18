#Load rdkit & version check
import rdkit
rdkit.__version__
import rdkit.Chem

###Reading the data### 
#reference: https://pubs.acs.org/doi/10.1021/acs.jcim.7b00083
#PubChemQC database contains about 200 million TD-DFT calculation results using the GAMESS QM package.
#The DB contains the following information
#Optimized molecular structures
#HOMO-LUMO gap
#Excitation state information, such as oscillator strength and energies of excited states
#osciilator strength ~ absoption coefficient

import pandas as pd
pubchem = pd.read_csv("https://www.dropbox.com/s/s9uhxw06z42gs8b/PubchemQC_subset_HOMO-LUMO_and_OS.csv?dl=1", on_bad_lines='skip')

#10만 개 data 사용
from tqdm import tqdm # progressive bar 표시를 위해서... 
gap_list = []
os_list = []
smi_list = []
ii = 0
max_mol = 100000 # 일단 10만개의 데이터를 사용하자... 
for gap, os, smi in zip(tqdm(pubchem["HOMO-LUMO_gap(eV)"]), pubchem["Oscillator_Strength"], pubchem["SMILES"]):
  gap = float(gap)
  os = float(os)
  gap_list.append(gap) # HOMO-LUMO gap
  os_list.append(os)  # Oscillator strength 값
  smi_list.append(smi) # 분자의 SMILES 표현식
  ii += 1
  if ii >= max_mol: 
    break
  else:
    continue

print(len(gap_list))
print(gap_list)
print(os_list)
print(smi_list[:10])

###Converting a molecule into an extended circular fingerprint (ECFP)
from rdkit.Chem import AllChem
mol = rdkit.Chem.MolFromSmiles(smi_list[0]) # SMILES -> Mol-type 변수로 읽어들인다. 
vec = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 1024) # ECFP2를 이용. 1024개의 0/1로 표현하겠다.

smi_list[0]
print(vec)
print(vec.ToList())

#RDKit을 이용해서 10만개 분자를 ECFP로 변환
X_ECFP = []
for smi in tqdm(smi_list):
  mol = rdkit.Chem.MolFromSmiles(smi)
  vec = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 1024)
  X_ECFP.append(vec.ToList()) # rdkit 최근 버젼.
  #X_ECFP.append([int(x) for x in vec.ToBitString()]) # rdkit 구 버젼 (2021 이전 버젼...)

print(len(X_ECFP))
print(X_ECFP[:2])
Y = gap_list # 예측하고자하는 목표 물성: HOMO-LUMO gap...

###Spliting training and test set###
#The ratio of training, validation, and test (holdout) data is a little bit arbitrary.
#The most common ratios are 6:2:2, 8:1:1, and 7:1:2 etc.

#Scikit-learn 로드
import sklearn
from sklearn.model_selection import train_test_split

#Train : Test = 8 : 2
train_X, test_X, train_y, test_y = train_test_split(X_ECFP, Y, test_size = 0.2, shuffle=True)
print(len(train_X))
print(len(test_X))

#Training set을 Training / Validation set으로 나눔.
train_X, val_X, train_y, val_y = train_test_split(train_X, train_y, test_size = 0.1, shuffle=True)

#MLP 모델을 구현
import torch
train_X = torch.tensor(train_X, dtype = torch.float32)
val_X = torch.tensor(val_X, dtype = torch.float32)
test_X = torch.tensor(test_X, dtype = torch.float32)

train_y = torch.tensor(train_y, dtype = torch.float32)
val_y = torch.tensor(val_y, dtype = torch.float32)
test_y = torch.tensor(test_y, dtype = torch.float32)

type(train_X)
train_X.shape

from torch.utils.data import TensorDataset, DataLoader

ds_train = TensorDataset(train_X, train_y)
ds_test = TensorDataset(test_X, test_y)

# batch size는 256
loader_train = DataLoader(ds_train, batch_size=256, shuffle=True, drop_last=True)
loader_test = DataLoader(ds_test, batch_size=256, shuffle=False)

# 전체 데이터의 개수. 
print(len(loader_train.dataset))
print(len(loader_test.dataset))

# 몇 개의 batch가 들어있는지? 
print(len(loader_train))
print(len(loader_test))

##Multi-layer perceptron을 이용해서 HOMO-LUMO gap을 예측
import torch.nn as nn

model = nn.Sequential(
    nn.Linear(1024, 128), # input_dim = 1024 // output dim = 128
    nn.LeakyReLU(0.1),    # activation function으로 leakyReLU 사용. https://pytorch.org/docs/stable/nn.html#non-linear-activations-weighted-sum-nonlinearity
    nn.BatchNorm1d(128),  # Batch Normalization
    nn.Linear(128, 64),   # 128 dim -> 64 dim
    nn.LeakyReLU(0.1),    # 
    nn.BatchNorm1d(64),   
    nn.Linear(64, 32),    # 64 -> 32 dim
    nn.LeakyReLU(0.1),
    nn.BatchNorm1d(32),
    nn.Linear(32, 1),     # 32 dim -> output_dim = 1 
    )

print(model)

from torch import optim

# 오차함수 선택 / Selecting a loss function
loss_fn = nn.MSELoss() # 회귀니까, mean squared error 를 선택. / Since this is a regression task, we will use mean squared error.

# 가중치를 학습하기 위한 최적화 기법 선택 / This is an optimizer to train weight parameters. 
optimizer = optim.Adam(model.parameters(), lr=0.01, weight_decay = 1e-5)

def train(epoch, X_val, y_val): 
    loss_list = [] # This list will save loss values 
    model.train()  # Model is set to a training mode

    # 데이터로더에서 미니배치를 하나씩 꺼내 학습을 수행 // At each step, 256 data points are extracted from the loader.
    for data, targets in loader_train: 

        optimizer.zero_grad()  # 경사를 0으로 초기화 / All gradients are set to zero.
        outputs = model(data)  # 데이터를 입력하고 오차 값을 계산 / The predicted values are calculated.
        outputs = outputs[:, 0] # targets와 모양이 동일하도록. / Convert 2D tensor -> 1D tensor so that it has the same shape with the "targets" tensor.

        loss = loss_fn(outputs, targets)  # loss 값 계산. / Calculating loss value.
        loss_list.append(loss.item()) # tensor 자체보다, 그 값만 필요할 경우. / Only values are stored. Otherwise, all necessary informations for gradient calculation will be contained in the tensor.

        loss.backward()   # 역전파 계산을 통한 기울기 값 계산 / Calculating gradient values.
        optimizer.step()  # 역전파 계산한 값으로 가중치를 수정 / Updating gradients with an optimizer.

    val_loss = loss_fn(model(X_val)[:,0], y_val) # calculate validation loss
    print(f"epoch {epoch} // Loss : {loss:.5f} // Val_Loss : {val_loss:.5f}")
    
    return loss.item(), val_loss.item()

###Test set에 대한 성능을 확인하는 함수
def test():
    model.eval()  # 신경망을 추론 모드로 전환
    mse = 0.0
    # 데이터로더에서 미니배치를 하나씩 꺼내 추론을 수행
    with torch.no_grad():  # 추론 과정에는 미분이 필요없음
        for data, targets in loader_test:
            outputs = model(data)  # 데이터를 입력하고 출력을 계산
            #print(outputs.shape)
            #print(targets.shape)

            # 추론 계산
            mse += torch.sum((outputs[:,0] - targets)**2)

    # 정확도 출력
    mse /= len(loader_test.dataset)  # 데이터 총 건수
    #mse = torch.sqrt(mse)
    print(f'\n테스트 데이터에서 예측 정확도: {torch.sqrt(mse).item():.3f} (eV)')

##100 epoch 동안 학습을 진행
n_epoch = 100
train_loss_list = []
val_loss_list = []
for epoch in range(n_epoch): # 100 epoch 동안 학습을 시키자!
    
    train_loss, val_loss = train(epoch, val_X, val_y)
    
    train_loss_list.append(train_loss)
    val_loss_list.append(val_loss)


