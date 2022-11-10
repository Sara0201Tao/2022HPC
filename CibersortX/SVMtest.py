import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.svm import SVC

#模型调用
f1 = open('optSVM.model','rb')
s1 = f1.read()
optSVM = pickle.loads(s1)
f1.close()

f2 = open('pca.model','rb')
s2 = f2.read()
pca = pickle.loads(s2)
f2.close()

#读取数据
bulk12 = pd.read_csv("bulk12_sub15_Sbatch_Results.csv")

bulk12_data = bulk12.values[:,2:17]
bulk12_data_pca = pca.transform(bulk12_data)
pred_bulk12 = optSVM.predict(bulk12_data_pca)

pred_bulk12