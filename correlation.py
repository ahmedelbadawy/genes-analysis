import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

#read csv files and create DataFrames
healthy = pd.read_csv('lusc-rsem-fpkm-tcga_paired.txt', sep='\t')
cancer = pd.read_csv('lusc-rsem-fpkm-tcga-t_paired.txt', sep='\t')
#remove rows with zeros greater than or equal to 25
h1 =healthy[(healthy == 0).sum(1) < 25]
c1 =cancer[(cancer == 0).sum(1) < 25]
li = list(set(h1['Hugo_Symbol']).intersection(set(c1['Hugo_Symbol'])))
h =h1[h1['Hugo_Symbol'].isin(li)]
c =c1[c1['Hugo_Symbol'].isin(li)]
h = h.reset_index(drop=True)
c = c.reset_index(drop=True)
# print(h.shape[0])
# print(c.shape[0])
###########################################################
#get Pearson correlation coefficients and insert them in a list
cc =[]
for i in range(h.shape[0]):
    h_g = h.iloc[i, 2:]
    c_g =c.iloc[i, 2:]
    r , x =pearsonr(h_g ,c_g)
    cc.append(r)
#create adataframe (corr)
corr = h[['Hugo_Symbol' , 'Entrez_Gene_Id']]
corr['cc'] = cc
#Rank genes based on their correlation coefficient (CC)
corr = corr.sort_values(by=['cc'])
print(corr)
#get the highest positive CC
max = corr[corr['cc'] == corr.cc.max()]
#get the lowest negative CC
min = corr[corr['cc'] == corr.cc.min()]
print (' The gene with the highest positive CC is ' ,max.iloc[0,0] , ", cc = ", max.iloc[0,2])
print(' The gene with the lowest negative CC is ' ,min.iloc[0,0] ,", cc = " ,min.iloc[0,2] )
#Plot the expression levels of the highest positive CC
plt.scatter(y = c[c['Hugo_Symbol'] == max.iloc[0,0]].iloc[:, 2:] , x = h[h['Hugo_Symbol'] == max.iloc[0,0]].iloc[:, 2:] , c= 'red' )
plt.ylabel('Expression levels for cancerous')
plt.xlabel('Expression levels for healthy')
plt.xticks(rotation=45)
plt.title('The highest positive CC')
plt.show()
plt.clf()
#Plot the expression levels of the lowest negative CC
plt.scatter(y = c[c['Hugo_Symbol'] == min.iloc[0,0]].iloc[:, 2:] , x = h[h['Hugo_Symbol'] == min.iloc[0,0]].iloc[:, 2:] , c= 'blue' )
plt.ylabel('Expression levels for cancerous')
plt.xlabel('Expression levels for healthy')
plt.xticks(rotation=45)
plt.title('The lowest negative CC')
plt.show()
