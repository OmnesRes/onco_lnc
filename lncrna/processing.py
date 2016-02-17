##this script will parse the mitranscriptome counts file for primary tumors and save results in separate files
##WARNING, I believe this script requires 6GB of RAM
##you will have to move the resulting files to their correct locations in tcga_data

f=open('mitranscriptome.expr.counts.tsv')
f.readline()
transcripts=[]
for i in f:
    transcripts.append(i.split('\t')[0].strip())

f=open('transcripts.txt','w')
f.write(str(transcripts))
f.close()

f=open('library_info.txt')
f.readline()
cancers=['READ', 'GBM','BLCA', 'UCEC', 'CESC', 'LIHC', 'HNSC', 'STAD', 'SKCM', 'COAD', 'LUAD', 'LUSC', 'OV', 'KIRP', 'LGG', 'LAML', 'BRCA', 'KIRC']
cancer_dict={}
file_to_patient={}

for i in f:
    x=i.split('\t')
    if x[10].strip() in cancers:
        if x[10].strip()=='LAML':
            if x[11].strip().split('-')[3][:2]=='03':
                cancer_dict[x[10].strip()]=cancer_dict.get(x[10].strip(),[])+[[x[11].strip(),x[0].strip()]]
                file_to_patient[x[0].strip()]=x[11].strip()
        elif x[10].strip()=='SKCM':
            if x[11].strip().split('-')[3][:2]=='01' or x[11].strip().split('-')[3][:2]=='06':
                cancer_dict[x[10].strip()]=cancer_dict.get(x[10].strip(),[])+[[x[11].strip(),x[0].strip()]]
                file_to_patient[x[0].strip()]=x[11].strip()
        else:
            if x[11].strip().split('-')[3][:2]=='01':
                cancer_dict[x[10].strip()]=cancer_dict.get(x[10].strip(),[])+[[x[11].strip(),x[0].strip()]]
                file_to_patient[x[0].strip()]=x[11].strip()



file_dict={}
for i in cancer_dict:
    for j in cancer_dict[i]:
        file_dict[j[1]]=i


cancer_to_indexes={}


                
f=open('mitranscriptome.expr.counts.tsv')
x=f.readline().split('\t')
for index, i in enumerate(x):
    if i.strip() in file_dict:
        cancer_to_indexes[file_dict[i.strip()]]=cancer_to_indexes.get(file_dict[i.strip()],[])+[[index,i.strip()]]

        
expression_dict={}

for j in f:
    x=j.split('\t')
    for k in cancer_to_indexes:
        temp=[]
        for index in cancer_to_indexes[k]:
            temp.append(x[index[0]])
        expression_dict[k]=expression_dict.get(k,[])+[temp]



for i in cancer_to_indexes:
    w=open(i+'.txt','w')
    for j in cancer_to_indexes[i]:
        w.write(file_to_patient[j[1]])
        w.write('\t')
    w.write('\n')
    for k in expression_dict[i]:
        w.write(str(k))
        w.write('\n')
    w.close()




            
