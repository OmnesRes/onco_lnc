## A script for finding every cox coefficient and pvalue for every miRNA in BRCA Tier 3 data downloaded Jan. 6th 2016


## Load necessary modules
from rpy2 import robjects as ro
import numpy as np
import os
ro.r('library(survival)')
import re


##This call will only work if you are running python from the command line.
##If you are not running from the command line manually type in your paths.
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))



## There were three clinical files with nonredundant data.  V4.0 is in general the most uptodate, but it is possible
## for data in the other files to be more uptodate.  As a result, clinical data will be merged.


f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt'))
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
        if re.search('^[0-9]+$',i[alive_column]):
            clinical1[-1]=[i[patient_column],int(i[alive_column]),'Alive']
        elif re.search('^[0-9]+$',i[death_column]):
            clinical1[-1]=[i[patient_column],int(i[death_column]),'Dead']
        else:
            pass
    else:
        if re.search('^[0-9]+$',i[alive_column]):
            clinical1.append([i[patient_column],int(i[alive_column]),'Alive'])
        elif re.search('^[0-9]+$',i[death_column]):
            clinical1.append([i[patient_column],int(i[death_column]),'Dead'])
        else:
            pass


## Removing the empty value.
clinical1=clinical1[1:]


f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt'))
##get the column indexes needed
columns=f.readline().split('\t')
patient_column=columns.index('bcr_patient_barcode')
alive_column=columns.index('last_contact_days_to')
death_column=columns.index('death_days_to')
f.readline()
f.readline()
data=[i.split('\t') for i in f]
clinical2=[['','','']]
for i in data:
    if clinical2[-1][0]==i[patient_column]:
        if re.search('^[0-9]+$',i[alive_column]):
            clinical2[-1]=[i[patient_column],int(i[alive_column]),'Alive']
        elif re.search('^[0-9]+$',i[death_column]):
            clinical2[-1]=[i[patient_column],int(i[death_column]),'Dead']
        else:
            pass
    else:
        if re.search('^[0-9]+$',i[alive_column]):
            clinical2.append([i[patient_column],int(i[alive_column]),'Alive'])
        elif re.search('^[0-9]+$',i[death_column]):
            clinical2.append([i[patient_column],int(i[death_column]),'Dead'])
        else:
            pass


##removing the empty value
clinical2=clinical2[1:]


##merging the data
new_clinical=[]

for i in clinical2:
    if i[0] not in [j[0] for j in clinical1]:
        new_clinical.append(i)
    else:
        if i[1]<=clinical1[[j[0] for j in clinical1].index(i[0])][1]:
            new_clinical.append(clinical1[[j[0] for j in clinical1].index(i[0])])
        else:
            new_clinical.append(i)



for i in clinical1:
    if i[0] not in [j[0] for j in new_clinical]:
        new_clinical.append(i)




f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_follow_up_v1.5_brca.txt'))
##get the column indexes needed
columns=f.readline().split('\t')
patient_column=columns.index('bcr_patient_barcode')
alive_column=columns.index('last_contact_days_to')
death_column=columns.index('death_days_to')
f.readline()
f.readline()
data=[i.split('\t') for i in f]
clinical3=[['','','']]
for i in data:
    if clinical3[-1][0]==i[patient_column]:
        if re.search('^[0-9]+$',i[alive_column]):
            clinical3[-1]=[i[patient_column],int(i[alive_column]),'Alive']
        elif re.search('^[0-9]+$',i[death_column]):
            clinical3[-1]=[i[patient_column],int(i[death_column]),'Dead']
        else:
            pass
    else:
        if re.search('^[0-9]+$',i[alive_column]):
            clinical3.append([i[patient_column],int(i[alive_column]),'Alive'])
        elif re.search('^[0-9]+$',i[death_column]):
            clinical3.append([i[patient_column],int(i[death_column]),'Dead'])
        else:
            pass


##removing the empty value
clinical3=clinical3[1:]


##merging the data
newer_clinical=[]

for i in clinical3:
    if i[0] not in [j[0] for j in new_clinical]:
        newer_clinical.append(i)
    else:
        if i[1]<=new_clinical[[j[0] for j in new_clinical].index(i[0])][1]:
            newer_clinical.append(new_clinical[[j[0] for j in new_clinical].index(i[0])])
        else:
            newer_clinical.append(i)

for i in new_clinical:
    if i[0] not in [j[0] for j in newer_clinical]:
        newer_clinical.append(i)


## Grade, sex, and age information were taken from the "clinical_patient" file.  A dictionary was created for sex and grade.
more_clinical={}
grade_dict={}
grade_dict['Infiltrating Ductal Carcinoma']=1
grade_dict['Metaplastic Carcinoma']=3
grade_dict['Mucinous Carcinoma']=4
grade_dict['Medullary Carcinoma']=5
grade_dict['Infiltrating Lobular Carcinoma']=6




sex_dict={}
sex_dict['MALE']=0
sex_dict['FEMALE']=1

## The "clinical_patient" file can also contain patients not listed in the follow_up files.
## In these cases the clinical data for these patients gets appended to a new clinical list.

f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','clinical','nationwidechildrens.org_clinical_patient_brca.txt'))
columns=f.readline().split('\t')
grade_column=columns.index('histological_type')
sex_column=columns.index('gender')
age_column=columns.index('age_at_diagnosis')
patient_column=columns.index('bcr_patient_barcode')
alive_column=columns.index('last_contact_days_to')
death_column=columns.index('death_days_to')
f.readline()
f.readline()
data=[i.split('\t') for i in f]
clinical4=[]
for i in data:
    try:
        more_clinical[i[patient_column]]=[grade_dict[i[grade_column]],sex_dict[i[sex_column]],int(i[age_column])]
        if re.search('^[0-9]+$',i[alive_column]):
            clinical4.append([i[patient_column],int(i[alive_column]),'Alive'])
        elif re.search('^[0-9]+$',i[death_column]):
            clinical4.append([i[patient_column],int(i[death_column]),'Dead'])
        else:
            pass

    except:
        pass


newest_clinical=[]

##It is possible that the clinical data in the clinical_patient file is more up to date than the follow_up files
##All the clinical data is merged checking which data is the most up to date

for i in clinical4:
    if i[0] not in [j[0] for j in newer_clinical]:
        newest_clinical.append(i)
    else:
        if i[1]<=newer_clinical[[j[0] for j in newer_clinical].index(i[0])][1]:
            newest_clinical.append(newer_clinical[[j[0] for j in newer_clinical].index(i[0])])
        else:
            newest_clinical.append(i)


##also do the reverse since clinical can contain patients not included in clinical4
for i in newer_clinical:
    if i[0] not in [j[0] for j in newest_clinical]:
        newest_clinical.append(i)


## only patients who had a follow up time greater than 0 days are included in the analysis
clinical=[i for i in newest_clinical if i[1]>0]


final_clinical=[]

## A new list containing both follow up times and grade, sex, and age is constructed.
## Only patients with grade, sex, and age information are included.
## Data is [[Patient ID, time (days), vital status, grade, sex, age at diagnosis],...]

for i in clinical:
    if i[0] in more_clinical:
        final_clinical.append(i+more_clinical[i[0]])



## Need to map the miRNA files to the correct patients
## The necessary information is included in the FILE_SAMPLE_MAP.txt file
f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','FILE_SAMPLE_MAP_mirna.txt'))
f.readline()
data=[i.strip().split() for i in f if i!='\n']


#### 01 indicates a primary tumor, and only primary tumors are included in this analysis
TCGA_to_mirna={}
for i in data:
    ##normalized files were used
    if 'isoform.quantification' in i[0]:
        if i[1].split('-')[3][:-1]=='01':
            x=''.join([k+j for k,j in zip(['','-','-'],i[1].split('-')[:3])])
            TCGA_to_mirna[x]=TCGA_to_mirna.get(x,[])+[i[0]]



clinical_and_files=[]
## I only care about patients that contained complete clinical information
for i in final_clinical:
    if TCGA_to_mirna.has_key(i[0]):
        ## The miRNA files are added to the clinical list
        ## Data structure: [[Patient ID, time (days), vital status, grade, sex, age at diagnosis,[miRNA files]],...]
        clinical_and_files.append(i+[TCGA_to_mirna[i[0]]])
    else:
        pass

## A list of lists of miRNAs is constructed, the order of miRNA lists is same as the clinical_and_files data
## The order of mirnas within the lists is defined by me (they are sorted).
## I use my reannonated read counts derived from the isoform files.
## Data structure: [[mirnas for patient 1], [mirnas for patient 2], ....]

f=open(os.path.join(BASE_DIR,'mirna','mirna_list.txt'))
mirna_list=[i.strip() for i in f]
mirnas=[]

for i in clinical_and_files:
    temp=[]
    for j in i[-1]:
        f=open(os.path.join(BASE_DIR,'tcga_data','BRCA','mirna',j.split('.txt')[0]+'new.txt'))
        mirna_dict={mirna:counts for mirna,counts in [[i.split()[0],float(i.strip().split()[-1])] for i in f]}
        temp.append([[mirna,mirna_dict.get(mirna,0)] for mirna in mirna_list])
    ## In the case that the patient only contained 1 primary tumor miRNA file.
    if len(temp)==1:
        mirnas.append(temp[0])
    ## If the patient contained more than 1 primary tumor miRNA file
    ## this list comprehension will average the files for any number of files.
    else:
        values=[]
        for k in temp:
            values.append([kk[1] for kk in k])
        mirnas.append(zip([z[0] for z in temp[0]],list(sum([np.array(kkk) for kkk in values])/float(len(temp)))))


## Only want mirnas that meet an expression cutoff
## A cutoff of .5 reads per million mirna mapped and no more than a fourth of the patients containing no expression was chosen
final_mirnas=[[]]*len(mirnas)
for i in range(len(mirnas[0])):
    temp=[]
    for j in mirnas:
        temp.append(j[i])
    count=0
    for k in temp:
        if k[1]==0:
            count+=1
    median=np.median([ii[1] for ii in temp])
    if count<len(mirnas)/4.0 and median>.5:
        for index, kk in enumerate(temp):
            final_mirnas[index]=final_mirnas[index]+[kk]


## This will write the final mirnas to a file (1-20 MB) which could be useful for further analyses, this step can be skipped.
f=open(os.path.join(BASE_DIR,'mirna','cox','BRCA','final_mirnas.txt'),'w')
for i in final_mirnas:
    f.write(str(i))
    f.write('\n')
f.close()

##Performing Cox regression on all of the mirnas in final_mirnas

death_dic={}
death_dic['Alive']=0
death_dic['Dead']=1



coeffs=[]
pvalues=[]
mirnas=[]  ##This list tracks the mirna names
for i in range(len(final_mirnas[0])):
    kaplan=[]
    mirnas.append(final_mirnas[0][i][0])
    for k,j in zip(clinical_and_files,final_mirnas):  ## These lists contain the clinical information and miRNA data in the same order.
        kaplan.append([k[1],k[2],k[3],k[4],k[5],j[i][1]])
    data=[ii[-1] for ii in kaplan]  ## Grabbing all the mirna values for the current mirna being analyzed
    ro.globalenv['expression']=ro.FloatVector(data)
    res=ro.r('round(qnorm((rank(expression, na.last="keep")-0.5)/sum(!is.na(expression))), digit=5)') ## Perform inverse normal transformation
    inverse_norm=list(res)  ## Convert robject to python list
    ## Prepare the variables for rpy2
    ro.globalenv['mirna']=ro.FloatVector(inverse_norm)
    ro.globalenv['times']=ro.IntVector([ii[0] for ii in kaplan])
    ro.globalenv['died']=ro.IntVector([death_dic[ii[1]] for ii in kaplan])
            
    ##ductal
    ductal=[]
    for ii in kaplan:
        if ii[2]==1:
            ductal.append(1)
        else:
            ductal.append(0)
            
    ##metaplastic
    metaplastic=[]
    for ii in kaplan:
        if ii[2]==3:
            metaplastic.append(1)
        else:
            metaplastic.append(0)


    ##mucinous
    mucinous=[]
    for ii in kaplan:
        if ii[2]==4:
            mucinous.append(1)
        else:
            mucinous.append(0)


    ##medullary
    medullary=[]
    for ii in kaplan:
        if ii[2]==5:
            medullary.append(1)
        else:
            medullary.append(0)

    ##lobular
    lobular=[]
    for ii in kaplan:
        if ii[2]==6:
            lobular.append(1)
        else:
            lobular.append(0)


    
    ro.globalenv['ductal']=ro.IntVector(ductal)
    ro.globalenv['metaplastic']=ro.IntVector(metaplastic)
    ro.globalenv['mucinous']=ro.IntVector(mucinous)
    ro.globalenv['medullary']=ro.IntVector(medullary)
    ro.globalenv['lobular']=ro.IntVector(lobular)
    ro.globalenv['sex']=ro.IntVector([ii[3] for ii in kaplan])
    ro.globalenv['age']=ro.IntVector([ii[4] for ii in kaplan])
    res=ro.r('coxph(Surv(times,died) ~ mirna + ductal + metaplastic + mucinous + medullary + lobular + sex + age)') ## Perform Cox regression
    ## Parse the string of the result with python for the mirna coefficient and pvalue
    for entry in str(res).split('\n'):
        try:
            if entry.split()[0]=='mirna':
                coeff=entry.split()[1]
                pvalue=entry.split()[-1]
                break
        except:
            pass
    coeffs.append(coeff)
    pvalues.append(pvalue)


## This will write the results to a tab delimited file with mirna name, cox coefficient, and pvalue.
f=open(os.path.join(BASE_DIR,'mirna','cox','BRCA','coeffs_pvalues.txt'),'w')
for i,j,k in zip(mirnas,coeffs,pvalues):
    f.write(i)
    f.write('\t')
    f.write(j)
    f.write('\t')
    f.write(k)
    f.write('\n')
f.close()




