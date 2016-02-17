## Add a column of BH adjusted pvalues for coeffs_pvalues.txt


## Load necessary modules
from rpy2 import robjects as ro

##open file and read data
f=open('coeffs_pvalues.txt')
data=[i.strip().split() for i in f]
genes,coeffs,pvalues=zip(*data)

#prepare data for R
ro.globalenv['pvalues']=ro.FloatVector([float(i) for i in pvalues])

#perform BH adjustment
res=ro.r('p.adjust(pvalues,"BH")')

#extract data
adjusted=list(res)

f=open('coeffs_pvalues_adjusted.txt','w')

for i,j,k,l in zip(genes,coeffs,pvalues,adjusted):
    f.write(i)
    f.write('\t')
    f.write(j)
    f.write('\t')
    f.write(k)
    f.write('\t')
    f.write(str(l))
    f.write('\n')
f.close()




