import rdkit
from rdkit import Chem
from json import load
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge
import matplotlib.pyplot as plt
from VGA.GroupAdd.Library import GroupLibrary
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

def get_descriptors(smi):
    mol=Chem.MolFromSmiles(smi)
    mol=Chem.AddHs(mol)
    descriptors={}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            around_atoms=[x.GetSymbol() for x in atom.GetNeighbors()]
            Hydrogen_count=0
            heavy_count=0
            for around_atom in around_atoms:
                if around_atom=='H':
                    Hydrogen_count += 1
                else:
                    heavy_count += 1
            descriptor=(Hydrogen_count,heavy_count,atom.GetSymbol())
            if descriptor not in descriptors.keys():
                descriptors[descriptor]=1
            else:
                descriptors[descriptor] += 1
        else:
            continue
    return descriptors

all_smis=load(open('p_smi_Ef_300.json'))[0]
all_Efs=load(open('p_smi_Ef_300.json'))[1]

#remove all the broken DFT calculations
error_num=[6, 8, 10, 16, 17, 26, 27, 28, 29, 31, 38, 39, 42, 45, 49, 54, 61, 72, 73, 81, 84, 87, 88, 89, 91, 95, 101, 103, 109, 114, 116, 119, 121, 132, 139, 141, 144, 145, 146, 147, 151, 156, 157, 159, 161, 164, 165, 168, 171, 173, 175, 177, 182, 183, 184, 185, 186, 191, 194, 195, 197, 199, 201, 202, 205, 206, 210, 212, 214, 215, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 235, 238, 242, 244, 247, 250, 251, 252, 253, 259, 265, 269, 273, 274, 275, 279, 283, 285, 286, 287, 289, 291, 292, 295, 299]
error_n=[]
for e in error_num:
    if e < 51:
        error_n.append(e)
    elif e > 66:
        error_n.append(e-3)
    elif e>=51 and e < 62:
        error_n.append(e-1)
    elif e >=62 and e < 66:
        error_n.append(e-2)

error_num=error_n[::-1]
for num in error_num:
    del all_smis[num]
    del all_Efs[num]
#num_all_descriptors=15
count=[]
train=[]
#print(len(all_smis))
for num,smi in enumerate(all_smis):
    descriptors=get_descriptors(smi)
    if False not in [descriptor in count for descriptor in descriptors.keys()]:
        #print('all_in: ',str(num))
        continue
    else:
        #print('new_descriptor: ',str(num))
        for descriptor in descriptors.keys():
            if descriptor not in count:
                count.append(descriptor)
        train.append(num)
#print(train)
#generate coefficient matrix
coefficient_matrix=np.zeros((len(all_smis),15))
descriptorss=[]
for i in range(len(all_smis)):
    smi=all_smis[i]
    descriptors=get_descriptors(smi)
    for descriptor in descriptors.keys():
        if descriptor not in descriptorss:
            descriptorss.append(descriptor)
        index=count.index(descriptor)
        coefficient_matrix[i][index]=descriptors[descriptor]
#print(count)
#print(all_smis[1])
#print('descriptor: ',descriptorss)
#print(coefficient_matrix[1])

#split train and test set
train_smi=[]
train_Ef=[]
train_matrix=[]

for t in train:
    train_smi.append(all_smis[t])
    train_Ef.append(all_Efs[t])
    train_matrix.append(coefficient_matrix[t])

train_0=train[::-1]
#print(train_smi)
for num in train_0:
    del all_smis[num]
    del all_Efs[num]
    coefficient_matrix=np.delete(coefficient_matrix,num,0)

reactions=load(open('reaction_300.json'))
mean=[]
max=[]
kf=KFold(n_splits=5)
for x in [0,0.1]:
    mean_sum=[]
    max_sum=[]
    mean_high,mean_low=[],[]
    max_high,max_low=[],[]
    #print(type(all_Efs))
    #print(all_Efs)
    #y_train,y_test=[],[]
    for i,j in kf.split(coefficient_matrix,all_Efs):
        #X_train,X_test,y_train,y_test=train_test_split(coefficient_matrix,all_Efs,test_size=0.3)
        #X_train=np.append(X_train,train_matrix,axis=0)
        #print(i)
        X_train,X_test=coefficient_matrix[i],coefficient_matrix[j]
        y_train,y_test=[],[]
        for k in i:
            y_train.append(all_Efs[k])
        for l in j:
            y_test.append(all_Efs[l])

        X_train=np.append(X_train,train_matrix,axis=0)
        
        for Ef in train_Ef:
            y_train.append(Ef)
        rge=Ridge(alpha=x,fit_intercept=False)
        rge_fit=rge.fit(X_train,y_train)

        #all_Efs == all_smis
        y_p = rge_fit.predict(coefficient_matrix)
        react_energies=[]
        t_energies=[]
        for num,reacts in enumerate(reactions[0]):
            pros=reactions[1][num]
            if reacts not in train_smi and True not in [pro in train_smi for pro in pros]:
                pro_idx_0=all_smis.index(pros[0])
                pro_idx_1=all_smis.index(pros[1])
                react_idx=all_smis.index(reacts)
                energy=y_p[pro_idx_0]+y_p[pro_idx_1]-y_p[react_idx]
                react_energies.append(energy)
                t_energies.append(reactions[2][num])
        mean_ae=abs(np.asarray(react_energies)-np.asarray(t_energies))
        max_ae=abs(np.asarray(react_energies)-np.asarray(t_energies))
        mean_sum.append(mean_ae)
        max_sum.append(max_ae)
    #mean_sum = np.mean(mean_sum)
    #print(mean_sum)
    #max_sum = np.mean(max_sum)
    #print(max_sum)
    mean_high,mean_low=np.max(mean_sum),np.min(mean_sum)
    max_high,max_low=np.max(max_sum),np.min(max_sum)
    mean_high_error,mean_low_error=mean_high-mean_sum,mean_sum-mean_low
    max_high_error,max_low_error=max_high-max_sum,max_sum-max_low
    #mean.append(mean_sum)
    #max.append(max_sum)
#print('mean_ae: ',mean_sum)
#print('mean_high_error: ',mean_high_error)
#print('mean_low_error: ',mean_low_error)
#print('max_ae: ',max_sum)
#print('max_high_error: ',max_high_error)
#print('max_low_error: ',max_low_error)

