###SMILES서열 및 docking score 기반으로 RNN docking socre predictor 생성 및 예측
# General Imports
import torch
import torch.nn as nn
import os
import pandas as pd
import numpy as np
import sklearn
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from matplotlib import pyplot as plt
%matplotlib inline

!cp ~/Downloads/dockingR_default_0529.csv .
docking_result= "dockingR_default_0529.csv" # data load
data = pd.read_csv(docking_result)
new_data = data[data['AD_GPU Top score'].notna()] ## NaN 제거 ## 임시로 AutoDock GPU score를 best_score라고 가정할 것. 데이터 완성시 best score 계산후 그 column을 여기에 치환
#new_data = data.applymap(lambda x: np.random.uniform(-10,0) if pd.isnull(x) else x) ### NaN에 랜덤값(-10~0) 삽입

#smiles_array = data["smiles"]
ptn_seq_array = new_data["Ptn_seq"] ##try
smiles_array = new_data["SMILES"]
target_prop = new_data["AD_GPU Top score"] # 최적 docking 결과 , 추후 bst_score로 치환

charset_ptn = set("".join(list(ptn_seq_array))+"!E")
charPtn_to_int = dict((c,i) for i,c in enumerate(charset_ptn))
int_to_charPtn = dict((i,c) for i,c in enumerate(charset_ptn))
maxPtn_len = max([len(ptn_seq) for ptn_seq in ptn_seq_array])
embed_ptn = maxPtn_len+5
print(charPtn_to_int)
print("Ptn_Character set: ", str(charset_ptn))
print(f"Number of all characters: {len(charset_ptn)}")
print(f"Length of padded Ptn_seq: {embed_ptn}")

charset = set("".join(list(smiles_array))+"!E")
char_to_int = dict((c,i) for i,c in enumerate(charset))
int_to_char = dict((i,c) for i,c in enumerate(charset))
max_len = max([len(smile) for smile in smiles_array])
embed = max_len+5
print(char_to_int)
print("Character set: ", str(charset))
print(f"Number of all characters: {len(charset)}")
print(f"Length of padded SMILES: {embed}")

def add_padding(ptn):
  new_ptn = 'E'*(embed-1 - len(ptn)) + ptn + '!'
  return new_ptn

ptnSeq_list = []
for ptn in ptn_seq_array:
  padded_ptn = add_padding(ptn)
  #print(smi)
  ptnSeq_list.append(padded_ptn)
print(ptnSeq_list[0])
print(ptnSeq_list[1])

def add_padding(smi):
  new_smi = 'E'*(embed-1 - len(smi)) + smi + '!'
  return new_smi

smiles_list = []
for smi in smiles_array:
  padded_smi = add_padding(smi)
  #print(smi)
  smiles_list.append(padded_smi)
print(smiles_list[0])
print(smiles_list[1])

from sklearn.model_selection import train_test_split
smiles_train, smiles_test, y_train, y_test = train_test_split(smiles_list, target_prop, random_state=42, test_size=0.2)
print("No. of training SMILES:", len(smiles_train))
print("No. of test SMILES:    ", len(smiles_test))


def vectorize(smiles):
  int_vec_list = []
  one_hot =  np.zeros((len(smiles), embed , len(charset)), dtype=np.int8)
  for i,smile in enumerate(smiles):
    #encode SMILES
    int_vec = []
    for j,c in enumerate(smile):
      try:
        one_hot[i,j,char_to_int[c]] = 1
        int_vec.append(char_to_int[c])
      except:
        one_hot[i,j,char_to_int['E']]=1
        int_vec.append(char_to_int['E'])
    int_vec_list.append(int_vec)
  return one_hot, int_vec_list

X_train, X_train_int = vectorize(smiles_train)
X_test, X_test_int = vectorize(smiles_test)
print(smiles_train[0])
print(X_train_int[0])
print(X_train[0])
plt.matshow(X_train[0].T)
print(X_train.shape)

from torch.utils.data import Dataset

class My_dataset(Dataset):
  def __init__(self, X, Y):
    self.X = torch.tensor(X, dtype=torch.float)
    self.Y = 1.0*torch.tensor(Y, dtype=torch.float)

  def __len__(self):
    return len(self.X)

  def __getitem__(self, idx):
    return self.X[idx, :,:], self.Y[idx]

dataset=My_dataset(X_train, np.array(y_train))
len(X_train)

train_loader = torch.utils.data.DataLoader(dataset=dataset,
                                          batch_size=8,
                                          shuffle=True, num_workers=1) ##batch size 조절해보기 4, 8, 16, 32

input_shape = X_train.shape[1:]
output_dim = 1
print("input_shape:", input_shape)
print("output_dim:", output_dim)
latent_dim = 64
lstm_dim = 64

import torch
import torch.nn as nn
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class Net(nn.Module):
    def __init__(self, dims, lstm_size, fc_hidden_size, dropout_rate, out_size):
        super(Net, self).__init__()
         
        length = dims[0]
        number_tokens = dims[1]
     
        self.lstm = nn.LSTM(input_size=number_tokens, hidden_size=lstm_size, num_layers=1, batch_first=True, bidirectional=False)
        self.fc1 = nn.Linear(lstm_size, fc_hidden_size) # Output layer
        self.activation = nn.ReLU() # Non-Linear ReLU Layer       
        self.fc_out = nn.Linear(fc_hidden_size, out_size) # Output layer
        self.dropout = nn.Dropout(dropout_rate)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x): # Forward pass: stacking each layer together
        out, (h_n, c_n) = self.lstm(x) #LSTM network reads in one-hot-encoded SMILES, h_n is last output, out is for all timesteps
        out = self.dropout(h_n) #Dropout
        out = self.fc1(out) # Pass into the hidden layer
        out = self.activation(out) # Use ReLU on hidden activation
        out = self.dropout(out) # dropout
        out = self.fc_out(out) # Use a linear layer for the output
        #out = self.sigmoid(out)
        return out
    
epochs = 200
dims = input_shape
lstm_size = 64      # The size of the LSTM layer
fc_hidden_size = 128   # The size of the hidden fully-connected layer
dropout_rate = 0.2  # The dropout rate
#output_size = len(charset)    # This is just a single task, so this will be one
output_size = 1
batch_size = 8       # The mini_batch size during training
learning_rate = 0.001  # The initial learning rate for the optimizer

model = Net(dims, lstm_size, fc_hidden_size, dropout_rate, output_size)
model.to(device)

print(X_train.shape)
print(X_train[0].shape)
print(device)

from torch.optim.lr_scheduler import ReduceLROnPlateau

#criterion = nn.CrossEntropyLoss()
criterion = nn.MSELoss()

optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

#lr_scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=50, 
#                  verbose=True, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=1e-6, eps=1e-08)

import sys
for x, y in train_loader:
  x = x.to(device)
  y = y.to(device)
  print(model(x), y)
  break


from tqdm import tqdm
epochs = 10 #epoch조절 가능
model.train() #Ensure the network is in "train" mode with dropouts active
train_losses = []
validation_losses = []

for e in range(epochs):
  running_loss = 0.0

  for input, labels in tqdm(train_loader):
      input = input.to(device)
      labels = labels.to(device)
      
      # Training pass
      optimizer.zero_grad() # Initialize the gradients, which will be recorded during the forward pass
      #print("input:")
      #print(input.shape)

      output = model(input) #Forward pass of the mini-batch
      output = output.squeeze(0).squeeze(-1)

      #print("output:")
      #print(output.shape)
      #print("labels:")
      #print(labels.shape)
      loss = criterion(output, labels) #Computing the loss
      loss.backward() # calculate the backward pass
      #print(loss.item())

      torch.nn.utils.clip_grad_norm_(model.parameters(), 1)
      optimizer.step() # Optimize the weights
        
      running_loss += loss.item()
  print(output)
  print(labels)
  print(f"Epoch: {e} Running Loss: {running_loss/len(train_loader):.4f}")

test_set = My_dataset(X_test, np.array(y_test))
test_loader = torch.utils.data.DataLoader(dataset=test_set, batch_size=1, shuffle=True)

preds=[]
truth=[]
def test(model, loader):
  model.eval()
  running_loss = 0.0
  for input, labels in tqdm(loader):
    input=input.to(device)
    labels=labels.to(device)
    #print(input.shape)
    output = model(input) # Forward pass of the mini-batch
    output = output.squeeze(0).squeeze(-1)
    loss = criterion(output, labels) #Computing the loss
    #print(output, labels)
    running_loss += loss.item()
    preds.append(output.item())
    truth.append(labels.item())
  return running_loss/len(loader), preds, truth

loss_val, preds, truth = test(model, test_loader)

import matplotlib
import matplotlib.pyplot as plt
print(preds)
print(truth)
plt.scatter(truth, preds, s=10.0, marker='.')
plt.xlim(xmin=-15, xmax=0)
plt.ylim(ymin=-15, ymax=0)
plt.xlabel("True docking score (kcal/mol)", fontsize='xx-large')
plt.ylabel("Predicted docking score (kcal/mol)", fontsize='xx-large')
plt.grid()
#plt.show()

pearson = np.corrcoef(preds, truth)[0,1]
print(f"Pearson Corr. Coef.: {pearson:.2f}")
