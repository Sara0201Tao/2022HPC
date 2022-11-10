import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc
from sklearn.manifold import TSNE

train_data_label = pd.read_csv('train_data_label.csv')
Train_data_label = train_data_label.values
np.random.shuffle(Train_data_label)
Train_data = Train_data_label[:,:15]
Train_label = Train_data_label[:,15]

#PCA
pca = PCA()
pca.fit(Train_data)
plt.plot(np.arange(pca.n_components_),pca.explained_variance_ratio_)

Train_data_pca = PCA(n_components=7).fit_transform(Train_data)

#设置SVM的不同kernel以及惩罚系数，通过交叉验证来选择最好的参数
K = ["linear", "poly", "rbf", "sigmoid"]
C = [0.1, 0.5, 1, 1.5, 2, 2.5, 3]

#SVM分类器，输入为样本的训练误差
def my_clf(kernel, C, train_data, train_label):
    clf = SVC(kernel=kernel, C=C)
    clf.fit(train_data,train_label)
    train_err = 1 - np.sum(clf.predict(train_data) == train_label)/len(train_label)
    mark = ''
    if train_err < 0.301 :
        mark = '*************'
    print("kernel:", kernel, "||  C = ", C, "||  train_error:", train_err, mark)
    return train_err

#五则交叉验证选取最优参数，误差小于0.34即正确率大于0.66的一组参数用***mark标记
#输入为原始数据PCA降维后的数据
saved_cv_err = []
for i in range(len(K)):
    for j in range(len(C)):
        clf = SVC(kernel=K[i], C=C[j])
        score = cross_val_score(clf, Train_data_pca, Train_label, cv=5)
        cv_err = 1 - sum(score)/5
        saved_cv_err.append(cv_err)
        mark = ''
        if cv_err < 0.34 :
            mark = '*************'
            
        print("kernel:", K[i], "||  C = ", C[j], "||  5_Fold_CV_error:", cv_err, mark)

#不同超参数的训练误差
saved_train_err = []

for i in range(len(K)):
    for j in range(len(C)):
        train_err = my_clf(kernel=K[i], C=C[j], train_data=Train_data_pca, train_label=Train_label)
        saved_train_err.append(train_err)


svm = SVC(kernel='sigmoid', C=1.5)
svm.fit(Train_data_pca, Train_label)
svm.score(Train_data_pca,Train_label) # 0.83333

saved_cv_err = []

for i in range(10):
    
    np.random.shuffle(Train_data_label)
    Train_data = Train_data_label[:,:15]
    Train_label = Train_data_label[:,15]
    Train_data_pca = PCA(n_components=7).fit_transform(Train_data)
    
    optSVM = SVC(kernel='sigmoid', C=1.5)
    score = cross_val_score(optSVM, Train_data_pca, Train_label, cv=5)
    print("5_Fold_CV_acc:", sum(score)/5)
    cv_err = 1 - sum(score)/5
    saved_cv_err.append(cv_err)
    
    optSVM.fit(Train_data_pca, Train_label)
    print("train_acc:", optSVM.score(Train_data_pca, Train_label), "\n")
    
print('***********************')
print('average 5 fold cv error:', sum(saved_cv_err)/len(saved_cv_err), 'training error: ', 1-optSVM.score(Train_data_pca, Train_label))



pca = PCA(n_components=7)
pca.fit(Train_data)

optSVM = SVC(kernel='sigmoid', C=1.5, probability=True)
pred_label = optSVM.fit(Train_data_pca, Train_label).decision_function(Train_data_pca)

fpr, tpr, thersholds = roc_curve(Train_label, pred_label)
roc_auc = auc(fpr, tpr)
fig_ROC = plt.figure()

plt.figure(figsize=(5,5))

plt.plot(fpr, tpr, color='darkorange',
         lw=2, label='ROC curve (AUC = %0.2f)' % roc_auc) # 假正率为横坐标，真正率为纵坐标做曲线
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([-0.003, 1.003])
plt.ylim([-0.003, 1.003])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve of SVM Classifier on train data')
plt.legend(loc="lower right")
# plt.savefig("ROC.png", dpi=1080)
plt.show()

optSVM = SVC(kernel='sigmoid', C=1.5, probability=True)
optSVM.fit(Train_data_pca, Train_label)

s = pickle.dumps(optSVM)
f = open('optSVM.model', "wb+")
f.write(s)
f.close()
print ("Done\n")

pca = PCA(n_components=7)
pca.fit(Train_data)

ss = pickle.dumps(pca)
ff = open('pca.model', "wb+")
ff.write(ss)
ff.close()
print ("Done\n")




