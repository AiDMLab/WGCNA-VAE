from sklearn import metrics
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df=pd.read_excel("C:\\Users\\dell\\Downloads\\vae\\combine.xlsx",header=0)
X=df.iloc[0:len(df),1:(len(df.columns)-1)]
Y=df.iloc[0:len(df),len(df.columns)-1]
train_x,test_x=X.iloc[0:173,:],X.iloc[173:238,:]
train_y,test_y=Y.iloc[0:173],Y.iloc[173:238]
#train_x,test_x,train_y,test_y=train_test_split(X,Y,test_size=0.4)
model=SVC(probability=True,gamma='auto',kernel="linear")
model.fit(train_x,train_y)
predict=model.predict(test_x)
confusion=metrics.confusion_matrix(test_y,predict)
accuracy=metrics.accuracy_score(test_y,predict)
fpr,tpr,thresholds= metrics.roc_curve(test_y,model.decision_function(test_x))
auc=metrics.auc(fpr,tpr)
print(accuracy)
print(auc)
plt.figure(figsize=(10,10))
plt.rc('font',family="Times New Roman")
plt.plot(fpr,tpr)
plt.plot([0,1],[0,1],color='navy',linestyle='--')
plt.title('ROC curve',fontsize=22,fontweight="bold")
plt.xlabel('False Positive Rate',fontsize=20,fontweight="bold")
plt.ylabel('True Positive Rate',fontsize=20,fontweight="bold")
plt.show()
