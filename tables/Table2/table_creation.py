##A script for creating a table

## Load necessary modules
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


##load data for each cancer, find total genes in oncolnc, get patient info
f=open(os.path.join(BASE_DIR,'mirna','cox','BLCA','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','BLCA','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
BLCA=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'mirna','cox','BRCA','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','BRCA','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
BRCA=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','CESC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','CESC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
CESC=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','COAD','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','COAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
COAD=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','ESCA','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','ESCA','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
ESCA=[genes_in_oncolnc]+data



##not all GBM data is included in oncolnc, need to count genes in oncolnc
f=open(os.path.join(BASE_DIR,'mirna','mature.fa'))
transcript_to_names={i.split()[1]:i.split()[0].strip('>') for i in f if '>' in i}
names_to_transcripts={i:j for i,j in zip(transcript_to_names.values(),transcript_to_names.keys())}


f=open(os.path.join(BASE_DIR,'mirna','cox','GBM','coeffs_pvalues.txt'))
data=[i.split() for i in f]
f2=open(os.path.join(BASE_DIR,'mirna','aliases.txt'))
aliases={}
for i in f2:
    aliases[i.strip().split()[1]]=i.split()[0]
names=[i[0] for i in data]
all_aliases={}
for i in names:
    for j in aliases:
            if i in j.split(';'):
                    all_aliases[i]=all_aliases.get(i,[])+[[j,aliases[j]]]
count=0
for i in data:
    if len(all_aliases[i[0]])==1:
        if all_aliases[i[0]][0][1] in transcript_to_names:
            count+=1
genes_in_oncolnc=count
f=open(os.path.join(BASE_DIR,'mirna','cox','GBM','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
GBM=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'mirna','cox','HNSC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','HNSC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
HNSC=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'mirna','cox','KIRC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','KIRC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
KIRC=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','KIRP','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','KIRP','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
KIRP=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','LAML','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','LAML','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LAML=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'mirna','cox','LGG','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','LGG','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LGG=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'mirna','cox','LIHC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','LIHC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LIHC=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'mirna','cox','LUAD','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','LUAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LUAD=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','LUSC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','LUSC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LUSC=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','SKCM','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','SKCM','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
SKCM=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','OV','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','OV','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
OV=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','PAAD','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','PAAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
PAAD=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','READ','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','READ','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
READ=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','SARC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','SARC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
SARC=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'mirna','cox','STAD','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','STAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
STAD=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'mirna','cox','UCEC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'mirna','cox','UCEC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
UCEC=[genes_in_oncolnc]+data


all_cancers=[BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LAML,LGG,LIHC,LUAD,LUSC,OV,PAAD,READ,SARC,SKCM,STAD,UCEC]

names=['BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD',\
       'READ','SARC','SKCM','STAD','UCEC']




f=open('table_2.txt','w')
for i,j in zip(all_cancers,names):
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







