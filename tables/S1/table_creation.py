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

##need to get the gene ids from a RNA-SEQV2 file, any file will work
f=open(os.path.join(BASE_DIR,'tcga_data','GBM','mrna','unc.edu.0cbec58e-f95e-4c60-a85d-210dc56bdf3c.1545137.rsem.genes.normalized_results'))
f.readline()
TCGA_id_to_gene={}
data=[i.split()[0] for i in f]
for i in data:
    TCGA_id_to_gene[i.split('|')[1]]=i.split('|')[0]


##this gene_result file is from http://www.ncbi.nlm.nih.gov/gene, downloaded Jan. 2016
f=open(os.path.join(BASE_DIR,'tables','gene_result.txt'))
f.readline()
current_id_to_gene={}
for i in f:
    x=i.split('\t')
    current_id_to_gene[x[2]]=x[5]


##I manually curated ids that got changed and created this file

new_ids={}
f=open(os.path.join(BASE_DIR,'tables','new_ids_annotated.txt'))
data=[i.strip().split() for i in f]
for i in data:
    if i[2]!='None':
        new_ids[i[0]]=[i[2],i[3]]
    else:
        new_ids[i[0]]='None'



for cancer in ['BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD',\
               'READ','SARC','SKCM','STAD','UCEC']:
    f=open(os.path.join(BASE_DIR,'mrna','cox',cancer,'coeffs_pvalues_adjusted.txt'))
    data=[i.strip().split() for i in f]
    ids,coeffs,pvalues,adjusted=zip(*data)
    cox_results=zip(ids,[TCGA_id_to_gene[i] for i in ids],coeffs,pvalues,adjusted)
    f=open(os.path.join(BASE_DIR,'mrna','cox',cancer,'final_genes.txt'))
    expression={}
    data=[eval(i.strip()) for i in f]
    for i in data:
        for j in i:
            expression[j[0]]=expression.get(j[0],[])+[j[1]]


    f=open(os.path.join(BASE_DIR,'tables','S1',cancer+'.txt'),'w')
    for i in sorted(cox_results,cmp=compare):
        f.write(i[0])
        f.write('\t')
        f.write(i[1])
        f.write('\t')
        if i[0] in current_id_to_gene:
            f.write(i[0])
            f.write('\t')
            f.write(current_id_to_gene[i[0]])
            f.write('\t')
        elif i[0] in new_ids:
            if new_ids[i[0]]=='None':
                f.write('None')
                f.write('\t')
                f.write('None')
                f.write('\t')
            else:
                f.write(new_ids[i[0]][0])
                f.write('\t')
                f.write(new_ids[i[0]][1])
                f.write('\t')
        else:
            raise
        f.write(i[2])
        f.write('\t')
        f.write(i[3])
        f.write('\t')
        f.write(i[4])
        f.write('\t')
        f.write(str(np.median(expression[i[0]])))
        f.write('\t')
        f.write(str(np.mean(expression[i[0]])))
        f.write('\n')
    f.close()


    del expression




