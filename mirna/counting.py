##a script for getting updated readcounts for human mirnas
##the id is used when available
##when unannotated the coordinates are used

f=open('human_mirnas.txt')
mirna_dict=eval(f.read())
f.close()
f=open('coordinates_dict.txt')
coordinates_dict=eval(f.read())
f.close()

##files.txt is a list of the files you want to recount
f=open('files.txt')
data=[i.split('.txt')[0] for i in f]

for k in data:
    f=open(k+'.txt')
    f.readline()
    count_dict={}
    for i in f:
        if 'MIMA' in i:
            x=i.strip().split('\t')
            try:
                count_dict[mirna_dict[x[-1].split(',')[1]]]=count_dict.get(mirna_dict[x[-1].split(',')[1]],0)+int(x[2])
            except:
                chromosome='chr'+x[1].split(':')[1]
                start=int(x[1].split(':')[2].split('-')[0])
                end=int(x[1].split(':')[2].split('-')[1])
                sign=x[1].split(':')[-1]
                overlaps=[]
                for ii in coordinates_dict:
                    for j in coordinates_dict[ii]:
                        if chromosome==j[0]:
                            if abs(start-j[1])<=3 and abs(end-j[2])<=3 and sign==j[3]:
                                overlaps.append(ii)
                                break
                if len(overlaps)==1:
                    count_dict[overlaps[0]]=count_dict.get(overlaps[0],0)+int(x[2])
                if len(overlaps)>1:
                    print overlaps, x
                    raise
        elif 'unannotated' in i:
            x=i.strip().split('\t')
            chromosome='chr'+x[1].split(':')[1]
            start=int(x[1].split(':')[2].split('-')[0])
            end=int(x[1].split(':')[2].split('-')[1])
            sign=x[1].split(':')[-1]
            overlaps=[]
            for ii in coordinates_dict:
                for j in coordinates_dict[ii]:
                    if chromosome==j[0]:
                        if abs(start-j[1])<=3 and abs(end-j[2])<=3 and sign==j[3]:
                            overlaps.append(ii)
                            break
            if len(overlaps)==1:
                count_dict[overlaps[0]]=count_dict.get(overlaps[0],0)+int(x[2])
            if len(overlaps)>1:
                print overlaps, x
                raise
        else:
            pass
    w=open(k+'new.txt','w')
    total=sum(count_dict.values())/1000000.0
    for i in count_dict:
        w.write(i)
        w.write('\t')
        w.write(str(count_dict[i]/total))
        w.write('\n')
    w.close()





    
