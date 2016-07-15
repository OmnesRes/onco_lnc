## A script for finding every cox coefficient and pvalue for every COAD lncRNA in the beta MiTranscriptome data set (normalized counts)

from rpy2 import robjects as ro
import numpy as np
import os
ro.r('library(survival)')
import re

##This call will only work if you are running python from the command line.
##If you are not running from the command line manually type in your paths.
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


f=open(os.path.join(BASE_DIR,'tcga_data','COAD','clinical','nationwidechildrens.org_clinical_follow_up_v1.0_coad.txt'))
##get the column indexes needed
columns=f.readline().split('\t')
patient_column=columns.index('bcr_patient_barcode')
alive_column=columns.index('last_contact_days_to')
death_column=columns.index('death_days_to')
f.readline()
f.readline()
data=[i.split('\t') for i in f]
## A patient can be listed multiple times in the file. The most recent listing (furthest down in the file), contains the most recent
## follow up data.  This code checks if the patient has already been loaded into the list, and if so, takes the more recent data.
## This required an empty value in the list initialization.
## Data is: [[Patient ID, time(days), Vital status],[Patient ID, time(days), Vital status],...]
clinical1=[['','','']]
for i in data:
    if clinical1[-1][0]==i[patient_column]:
        if re.search('^[0-9]+$',i[death_column]):
            clinical1[-1]=[i[patient_column],int(i[death_column]),'Dead']
        elif re.search('^[0-9]+$',i[alive_column]):
            clinical1[-1]=[i[patient_column],int(i[alive_column]),'Alive']
        else:
            pass
    else:
        if re.search('^[0-9]+$',i[death_column]):
            clinical1.append([i[patient_column],int(i[death_column]),'Dead'])
        elif re.search('^[0-9]+$',i[alive_column]):
            clinical1.append([i[patient_column],int(i[alive_column]),'Alive'])
        else:
            pass


## Removing the empty value.
clinical=clinical1[1:]


## Sex and age information were taken from the "clinical_patient" file.  A dictionary was created for sex.
more_clinical={}



sex_dict={}
sex_dict['MALE']=0
sex_dict['FEMALE']=1

## The "clinical_patient" file can also contain patients not listed in the follow_up files.
## In these cases the clinical data for these patients gets appended to a new clinical list.

f=open(os.path.join(BASE_DIR,'tcga_data','COAD','clinical','nationwidechildrens.org_clinical_patient_coad.txt'))
##get the column indexes needed
columns=f.readline().split('\t')
sex_column=columns.index('gender')
age_column=columns.index('age_at_initial_pathologic_diagnosis')
patient_column=columns.index('bcr_patient_barcode')
alive_column=columns.index('last_contact_days_to')
death_column=columns.index('death_days_to')
f.readline()
f.readline()
clinical4=[]
data=[i.split('\t') for i in f]
for i in data:
    try:
        more_clinical[i[patient_column]]=[0,sex_dict[i[sex_column]],int(i[age_column])]
        if re.search('^[0-9]+$',i[death_column]):
            clinical4.append([i[patient_column],int(i[death_column]),'Dead'])
        elif re.search('^[0-9]+$',i[alive_column]):
            clinical4.append([i[patient_column],int(i[alive_column]),'Alive'])
        else:
            pass
    except:
        pass


new_clinical=[]

##It is possible that the clinical data in the clinical_patient file is more up to date than the follow_up files
##All the clinical data is merged checking which data is the most up to date


for i in clinical4:
    if i[0] not in [j[0] for j in clinical]:
        new_clinical.append(i)
    else:
        if i[1]<=clinical[[j[0] for j in clinical].index(i[0])][1]:
            new_clinical.append(clinical[[j[0] for j in clinical].index(i[0])])
        else:
            new_clinical.append(i)

##also do the reverse since clinical can contain patients not included in clinical4
for i in clinical:
    if i[0] not in [j[0] for j in new_clinical]:
        new_clinical.append(i)


## only patients who had a follow up time greater than 0 days are included in the analysis
clinical=[i for i in new_clinical if i[1]>0]


final_clinical=[]

## A new list containing both follow up times and sex and age is constructed.
## Only patients with sex and age information are included.
## Data is [[Patient ID, time (days), vital status, 0, sex, age at diagnosis],...]


for i in clinical:
    if i[0] in more_clinical:
        final_clinical.append(i+more_clinical[i[0]])



##In a separate script I parsed the mitranscriptome.expr.counts.tsv file and extracted the COAD patient and expression values.
##From this file I will load the expression data.
##There are duplicated transcripts and the possibility of a patient having multiple sequencing files.

##create a dictionary to check for duplicated data
lncrna_dict={}
##I have the list of transcripts saved in a file
f=open(os.path.join(BASE_DIR,'lncrna','transcripts.txt'))
transcripts=eval(f.read())

f=open(os.path.join(BASE_DIR,'tcga_data','COAD','lncrna','COAD.txt'))
##patient list is at the top of the file

patients=f.readline().strip().split()
lncrnas=[[]]*len(patients)
for i,j in zip(transcripts,f):
    if i not in lncrna_dict:
        data=eval(j.strip())
        for index, k in enumerate(data):
            lncrnas[index]=lncrnas[index]+[[i,float(k)]]
        lncrna_dict[i]=''


##create a dictionary mapping patient to all of their lncrna expression data
patient_dict={}
for index, i in enumerate(patients):
    patient_dict[i[:12]]=patient_dict.get(i[:12],[])+[lncrnas[index]]


##find which patients have complete clinical data, order the data, and average data if necessary
##it's possible there are expression data for patients without clinical data, and clinical data without expression data

##create a new clinical list called clinical_and_files for consistency with previous scripts
clinical_and_files=[]
for i in final_clinical:
    if i[0] in patient_dict:
        clinical_and_files.append(i)


ordered_lncrnas=[]
for i in clinical_and_files:
    temp=[]
    for j in patient_dict[i[0]]:
        temp.append(j)
    if len(temp)==1:
        ordered_lncrnas.append(temp[0])
    else:
        values=[]
        for k in temp:
            values.append([kk[1] for kk in k])
        ordered_lncrnas.append(zip([z[0] for z in temp[0]],list(sum([np.array(kkk) for kkk in values])/float(len(temp)))))


## Only want lncras that meet an expression cutoff
## It is not known what expression level of lncrnas is needed for function, so a soft value for median was chosen.
## I don't want to perform an analysis with all 0 expression however, so zeros are still counted.
## A cutoff of .1 and no more than a fourth of the patients containing no expression was chosen
final_lncrnas=[[]]*len(ordered_lncrnas)
for i in range(len(ordered_lncrnas[0])):
    temp=[]
    for j in ordered_lncrnas:
        temp.append(j[i])
    count=0
    for k in temp:
        if k[1]==0:
            count+=1
    median=np.median([ii[1] for ii in temp])
    if count<len(ordered_lncrnas)/4.0 and median>.1:
        for index, kk in enumerate(temp):
            final_lncrnas[index]=final_lncrnas[index]+[kk]


## This will write the final lncrnas to a medium sized file ~10-50MB which could be useful for further analyses, this step can be skipped.
f=open(os.path.join(BASE_DIR,'lncrna','cox','COAD','final_lncrnas.txt'),'w')
for i in final_lncrnas:
    f.write(str(i))
    f.write('\n')
f.close()




##Performing Cox regression on all of the lncrnas in final_lncrnas

death_dic={}
death_dic['Alive']=0
death_dic['Dead']=1

coeffs=[]
pvalues=[]
lncrnas=[] ##This list tracks the lncrna names
for i in range(len(final_lncrnas[0])):
    kaplan=[]
    lncrnas.append(final_lncrnas[0][i][0])
    for k,j in zip(clinical_and_files,final_lncrnas): ## These lists contain the clinical information and lncrna data in the same order.
        kaplan.append([k[1],k[2],k[3],k[4],k[5],j[i][1]])
    data=[ii[-1] for ii in kaplan] ## Grabbing all the lncrna values for the current lncrna being analyzed
    ro.globalenv['expression']=ro.FloatVector(data)
    res=ro.r('round(qnorm((rank(expression, na.last="keep")-0.5)/sum(!is.na(expression))), digit=5)') ## Perform inverse normal transformation
    inverse_norm=list(res) ## Convert robject to python list
    ## Prepare the variables for rpy2
    ro.globalenv['lncrna']=ro.FloatVector(inverse_norm)
    ro.globalenv['times']=ro.IntVector([ii[0] for ii in kaplan])
    ro.globalenv['died']=ro.IntVector([death_dic[ii[1]] for ii in kaplan])
    ro.globalenv['sex']=ro.IntVector([ii[3] for ii in kaplan])
    ro.globalenv['age']=ro.IntVector([ii[4] for ii in kaplan])
    res=ro.r('coxph(Surv(times,died) ~ lncrna + sex + age)')  ## Perform Cox regression
    ## Parse the string of the result with python for the lncrna coefficient and pvalue
    for entry in str(res).split('\n'):
        try:
            if entry.split()[0]=='lncrna':
                coeff=entry.split()[1]
                pvalue=entry.split()[-1]
                break
        except:
            pass
    coeffs.append(coeff)
    pvalues.append(pvalue)


## This will write the results to a tab delimited file with lncrna name, cox coefficient, and pvalue.
f=open(os.path.join(BASE_DIR,'lncrna','cox','COAD','coeffs_pvalues.txt'),'w')
for i,j,k in zip(lncrnas,coeffs,pvalues):
    f.write(i)
    f.write('\t')
    f.write(j)
    f.write('\t')
    f.write(k)
    f.write('\n')
f.close()




