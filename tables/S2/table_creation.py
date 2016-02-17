import numpy as np

##A script for creating tables for each cancer, with the data sorted 

def compare(first,second):
    if float(first[-2])>float(second[-2]):
        return 1
    elif float(first[-2])<float(second[-2]):
        return -1
    else:
        return 0

count=0
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

##create necessary dictionaries

##get annotations for micrornas
##I created this file from mirbase mature.fa
f=open(os.path.join(BASE_DIR,'mirna','human_mirnas.txt'))
mirna_dict=eval(f.read())
names_to_transcripts={i:j for i,j in zip(mirna_dict.values(),mirna_dict.keys())}


for cancer in ['BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD',\
               'READ','SARC','SKCM','STAD','UCEC']:
    f=open(os.path.join(BASE_DIR,'mirna','cox',cancer,'coeffs_pvalues_adjusted.txt'))
    cox_results=[i.strip().split() for i in f]
    f=open(os.path.join(BASE_DIR,'mirna','cox',cancer,'final_mirnas.txt'))
    expression={}
    data=[eval(i.strip()) for i in f]
    for i in data:
        for j in i:
            expression[j[0]]=expression.get(j[0],[])+[j[1]]
    f=open(os.path.join(BASE_DIR,'tables','S2',cancer+'.txt'),'w')
    if cancer!='GBM':
        for i in sorted(cox_results,cmp=compare):
            f.write(i[0])
            f.write('\t')
            f.write(names_to_transcripts[i[0]])
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
    else:
        ##GBM data was from a microarray without gene identifiers, just names
        ##To update these miRNAs I had to get all possible ids from the mature.fa file and link names to ids with the aliases file
        f2=open(os.path.join(BASE_DIR,'mirna','mature.fa'))
        mature={i.split()[1]:i.split()[0].strip('>') for i in f2 if '>' in i}
        f2=open(os.path.join(BASE_DIR,'mirna','aliases.txt'))
        aliases={}
        for i in f2:
            aliases[i.strip().split()[1]]=i.split()[0]
        names=[i[0] for i in cox_results]
        all_aliases={}
        for i in names:
            for j in aliases:
                    if i in j.split(';'):
                            all_aliases[i]=all_aliases.get(i,[])+[[j,aliases[j]]]
        for i in sorted(cox_results,cmp=compare):
            if len(all_aliases[i[0]])==1:
                if all_aliases[i[0]][0][1] in  mature:
                    count+=1
                    f.write(mature[all_aliases[i[0]][0][1]])
                    f.write('\t')
                    f.write(all_aliases[i[0]][0][1])
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
                else:
                    ##these are either deleted or stem loop sequences, check manually
                    print i[0],all_aliases[i[0]][0][1]
                    f.write(i[0])
                    f.write('\t')
                    f.write('deleted')
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
            else:
                ##these have several possible current names
                f.write(i[0]+'?')
                f.write('\t')
                f.write(''.join([l+m for l,m in zip([k[1] for k in all_aliases[i[0]]],['\\']*len(all_aliases[i[0]]))]))
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

print count
