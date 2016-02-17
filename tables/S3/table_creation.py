import numpy as np

##A script for creating tables for each cancer, with the data sorted 

def compare(first,second):
    if float(first[-2])>float(second[-2]):
        return 1
    elif float(first[-2])<float(second[-2]):
        return -1
    else:
        return 0


import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

##create necessary dictionaries

##get annotations for lncrnas
##this file was downloaded from mitranscriptome.org
f=open(os.path.join(BASE_DIR,'lncrna','mitranscriptome.gtf','mitranscriptome.gtf','mitranscriptome.v2.gtf','mitranscriptome.v2.gtf'))
transcript_dict={}
for i in f:
    transcript_dict[i.split('transcript_id')[1].strip().split('";')[0].strip('"')]=i.strip().split()[-1].split('";')[0].strip('"')


for cancer in ['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV',\
               'READ','SKCM','STAD','UCEC']:
    f=open(os.path.join(BASE_DIR,'lncrna','cox',cancer,'coeffs_pvalues_adjusted.txt'))
    cox_results=[i.strip().split() for i in f]
    f=open(os.path.join(BASE_DIR,'lncrna','cox',cancer,'final_lncrnas.txt'))
    expression={}
    data=[eval(i.strip().replace('nan',"'nan'")) for i in f]
    for i in data:
        for j in i:
            if j[1]!='nan':
                expression[j[0]]=expression.get(j[0],[])+[j[1]]
    f=open(os.path.join(BASE_DIR,'tables','S3',cancer+'.txt'),'w')
    for i in sorted(cox_results,cmp=compare):
        f.write(i[0])
        f.write('\t')
        f.write(transcript_dict[i[0]])
        f.write('\t')
        f.write(i[1])
        f.write('\t')
        f.write(i[2])
        f.write('\t')
        f.write(i[3])
        f.write('\t')
        f.write(str(np.median(expression[i[0]])))
        f.write('\t')
        f.write(str(np.mean(expression[i[0]])))
        f.write('\n')
    f.close()
    del expression




