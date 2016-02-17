##A script for creating a table

import numpy as np
## Load necessary modules
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def compare(first,second):
    if float(first[-2])>float(second[-2]):
        return 1
    elif float(first[-2])<float(second[-2]):
        return -1
    else:
        return 0

##load data for each cancer, find total genes in oncolnc, get patient info

f=open(os.path.join(BASE_DIR,'tcga_data','GBM','mrna','unc.edu.0cbec58e-f95e-4c60-a85d-210dc56bdf3c.1545137.rsem.genes.normalized_results'))
f.readline()
TCGA_id_to_gene={}
data=[i.split()[0] for i in f]
for i in data:
    TCGA_id_to_gene[i.split('|')[1]]=i.split('|')[0]



f=open(os.path.join(BASE_DIR,'tables','gene_result.txt'))
f.readline()
current_id_to_gene={}
for i in f:
    x=i.split('\t')
    current_id_to_gene[x[2]]=x[5]


new_ids={}
f=open(os.path.join(BASE_DIR,'tables','new_ids_annotated.txt'))
data=[i.strip().split() for i in f]
for i in data:
    if i[2]!='None':
        new_ids[i[0]]=i[2]

all_ids={}
all_genes={}
allowed_ids={}

##this table only records genes present in oncolnc, and only genes which have a clear current id are allowed in oncolnc
for i in TCGA_id_to_gene:
    if i in current_id_to_gene:
        all_ids[i]=''
        all_genes[current_id_to_gene[i]]=''
        allowed_ids[i]=[i,current_id_to_gene[i]]

for i in TCGA_id_to_gene:      
    if i in new_ids:
        if new_ids[i] not in all_ids:
            all_ids[new_ids[i]]=''
            all_genes[current_id_to_gene[new_ids[i]]]=''
            allowed_ids[i]=[new_ids[i],current_id_to_gene[new_ids[i]]]
        else:
            pass
    else:
        pass


all_cancers=[]

cancers=['BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV',\
       'PAAD','READ','SARC','SKCM','STAD','UCEC']
for cancer in cancers:
    f=open(os.path.join(BASE_DIR,'mrna','cox',cancer,'coeffs_pvalues_adjusted.txt'))
    cox_results=[i.strip().split() for i in f]
    count=0
    for index,i in enumerate(sorted(cox_results,cmp=compare)):
        if i[0] in allowed_ids:
            count+=1
    f=open(os.path.join(BASE_DIR,'mrna','cox',cancer,'coeffs_pvalues.txt'))
    data=[i for i in f]
    genes_in_oncolnc=count
    f=open(os.path.join(BASE_DIR,'mrna','cox',cancer,'patient_info.txt'))
    f.readline()
    data=f.readline().strip().split()
    all_cancers.append([genes_in_oncolnc]+data)




f=open('table_1.txt','w')
for i,j in zip(all_cancers,cancers):
    f.write(j)
    f.write('\t')
    ##write total patients (add males and females)
    f.write(str(int(i[2])+int(i[3])))
    f.write('\t')
    ##write male/female
    f.write(i[2]+'/'+i[3])
    f.write('\t')
    ##write average age at diagnosis
    f.write(i[1])
    f.write('\t')
    ##write events
    f.write(i[4])
    f.write('\t')
    ##write median survival
    f.write(i[5])
    f.write('\t')
    ##write genes in oncolnc
    f.write(str(i[0]))
    f.write('\t')
    f.write('\n')
f.close()







