import pandas as pd
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
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
rel_val =[]
ind_val =[]
# print(h.shape[0])
# print(c.shape[0])
#get p-values for each paired and ind. and insert them in a list
for i in range(h.shape[0]):
    h_g = h.iloc[i, 2:]
    c_g =c.iloc[i, 2:]
    #to make Samples paired
    p_val_rel = ttest_rel(h_g, c_g).pvalue
    #to make Samples indpendant
    p_val_ind = ttest_ind(h_g, c_g).pvalue
    rel_val.append(p_val_rel)
    ind_val.append(p_val_ind)
# print(len(rel_val))
# print(len(ind_val))
#Apply the FDR multiple tests correction method
p_relval_fdr = multipletests(rel_val, alpha=0.05, method='fdr_bh')[1]
p_indval_fdr = multipletests(ind_val, alpha=0.05, method='fdr_bh')[1]

#get the list of DEGs before and after the FDR correction for Samples paired
sign_relg = pd.DataFrame({'Hugo_Symbol':h['Hugo_Symbol'].tolist(), 'p_values':rel_val, 'p_values_fdr':p_relval_fdr})
sign_relg['significance:p_vlaue'] = sign_relg['p_values'].apply(lambda x: x < 0.05)
sign_relg['significance:p_vlaue_fdr'] = sign_relg['p_values_fdr'].apply(lambda x: x < 0.05)
DIG_rel = sign_relg[sign_relg['significance:p_vlaue']== True]
rel_before = DIG_rel['Hugo_Symbol'].tolist()
DIG_rel_fdr = sign_relg[sign_relg['significance:p_vlaue_fdr']== True]
rel_after = DIG_rel_fdr['Hugo_Symbol'].tolist()
#print the num. of DEGs
print('There were ',len(rel_before),' DEGs in case of paired samples before FDR correction' )
print('There were ',len(rel_after),' DEGs in case of paired samples after FDR correction')
###################################################################################
#get the list of DEGs before and after the FDR correction for Samples indpendant
sign_indg = pd.DataFrame({'Hugo_Symbol':h['Hugo_Symbol'].tolist(), 'p_values':ind_val, 'p_values_fdr':p_indval_fdr})
sign_indg['significance:p_vlaue'] = sign_indg['p_values'].apply(lambda x: x < 0.05)
sign_indg['significance:p_vlaue_fdr'] = sign_indg['p_values_fdr'].apply(lambda x: x < 0.05)
DIG_ind = sign_indg[sign_indg['significance:p_vlaue']== True]
ind_before = DIG_ind['Hugo_Symbol'].tolist()
DIG_ind_fdr = sign_indg[sign_indg['significance:p_vlaue_fdr']== True]
ind_after = DIG_ind_fdr['Hugo_Symbol'].tolist()
#print the num. of DEGs
print('There were ',len(ind_before),' DEGs in case of indpendant samples before FDR correction' )
print('There were ',len(ind_after),' DEGs in case of indpendant samples after FDR correction')


############################################################################################
#Compare the two DEGs sets (paired and independent) after the FDR correction
rel_diff = list(set(rel_after).difference(set(ind_after)))
ind_diff = list(set(ind_after).difference(set(rel_after)))
inter = list(set(rel_after).intersection(set(ind_after)))
#print the num. of DEGs
print('There were ',len(rel_diff),' distinct DEGs in the paired samples after FDR correction.')
print('There were ',len(ind_diff), ' distinct DEGs in the indpendant samples after FDR correction.')
print('There were ',len(inter), ' common DEGs between the paired samples and the independent samples after FDR correction')
