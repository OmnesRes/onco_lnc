##A script for creating a table


## Load necessary modules
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


##load data for each cancer, find total genes in oncolnc, get patient info
f=open(os.path.join(BASE_DIR,'lncrna','cox','BLCA','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','BLCA','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
BLCA=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','BRCA','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','BRCA','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
BRCA=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','CESC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','CESC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
CESC=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','COAD','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','COAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
COAD=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','GBM','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','GBM','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
GBM=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','HNSC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','HNSC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
HNSC=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','KIRC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','KIRC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
KIRC=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','KIRP','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','KIRP','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
KIRP=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','LAML','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','LAML','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LAML=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','LGG','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','LGG','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LGG=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','LIHC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','LIHC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LIHC=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','LUAD','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','LUAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LUAD=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','LUSC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','LUSC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
LUSC=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','SKCM','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','SKCM','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
SKCM=[genes_in_oncolnc]+data

f=open(os.path.join(BASE_DIR,'lncrna','cox','OV','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','OV','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
OV=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','READ','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','READ','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
READ=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','STAD','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','STAD','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
STAD=[genes_in_oncolnc]+data


f=open(os.path.join(BASE_DIR,'lncrna','cox','UCEC','coeffs_pvalues.txt'))
data=[i for i in f]
genes_in_oncolnc=len(data)
f=open(os.path.join(BASE_DIR,'lncrna','cox','UCEC','patient_info.txt'))
f.readline()
data=f.readline().strip().split()
UCEC=[genes_in_oncolnc]+data


all_cancers=[BLCA,BRCA,CESC,COAD,GBM,HNSC,KIRC,KIRP,LAML,LGG,LIHC,LUAD,LUSC,OV,READ,SKCM,STAD,UCEC]

names=['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV',\
       'READ','SKCM','STAD','UCEC']

f=open('table_3.txt','w')
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







