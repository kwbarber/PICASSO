"""
File to analyze raw data (exported .txt file from microarray image by GenePix)
for double base pair substitutions in dsDNA sequence targeted by sgRNA
as in Fig. 2C and Fig. S2D
@author: kwbarber
"""

import csv

barcode='CTCACGGAGAAATCTATGGA' ##sgRNA spacer barcode

afile= ##enter file to be analyzed

imat=[]
file_reader=open(afile)
csv_reader=csv.reader(file_reader,delimiter=',')
for i in csv_reader:        
    imat.append(i)


Name=imat[0].index('Name')
Sequence=imat[0].index('Sequence')
Intensity=imat[0].index('F635 Median - B635')


lmat=['A','C','G','T']

nmat=[] ##imports probe sequences for normalization, here using on-target sequence sans PAM
for i in imat[1:]:
    if 'F-NP' in i[Name] and barcode in i[Sequence]:
        nmat.append([i[Name],i[Sequence][-34:-14],i[Intensity]])
norm_factor=0
for i in nmat:
    norm_factor+=int(i[2])
norm_factor=norm_factor/len(nmat)
if norm_factor<0:
    norm_factor=0

wmat=[]
for i in imat[1:]:
    if 'F' in i[Name] and 'vmut' in i[Name]:
        wmat.append([i[Name],i[Sequence][-34:-14],i[Intensity]])


WT_list=[] ##this is WT sequence intensity
for i in imat[1:]:
    if 'F-sub-1' in i[Name]:
        if i[Sequence][-34:-14]==barcode:
            WT_list.append(int(i[Intensity])-norm_factor) ##normalized using NP sequence

pmat=[]
for i in wmat:
    counter=0
    pcount=[]
    for j in barcode:
        if i[1][counter]!=j:
            pcount.append(counter+1)
        counter+=1
    if int(i[2]) - norm_factor >= 0:  ##normalized intensity using NP probe control
        pmat.append([pcount,int(i[2]) - norm_factor])  ##normalized intensity using NP probe control
    else:
        pmat.append([pcount,0])  ##negative values converted to zero

    
amat=[]
for i in pmat:
    breaker=0
    for j in pmat:
        if i[0]==j[0] and i[1]!=j[1]:
            aval=(int(i[1])+int(j[1]))/2
            amat.append([i[0],aval])
            breaker=1
    if breaker==0:
        amat.append([i[0],int(i[1])])
amat.sort()

final_list=[]
for i in amat:
    if i not in final_list:
        final_list.append(i)
        
final_mat=[]
counter=1
while counter<21:
    temp_list=[]
    for i in final_list:
        if i[0][1]==counter:
            temp_list.append(i[1])
    final_mat.append(temp_list)
    counter+=1

max_val=0 ##use this and final_norm_mat for values normalized to maximum value
for i in final_mat:
    for j in i:
        if j>max_val:
            max_val=j

final_norm_mat=[]
for i in final_mat:
    temp_mat=[]
    for j in i:
        temp_mat.append(j/max_val)
    final_norm_mat.append(temp_mat)
    
##replace 'xxx' below with file destination
outWrite=open('xxx.csv','w')
counter=0
for i in final_norm_mat[1:]:
    outWrite.write(str(i).replace(']','').replace('[','').replace("'","").replace(' ','')+'\n')
outWrite.close()



