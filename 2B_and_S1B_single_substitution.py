"""
File to analyze raw data (exported .txt file from microarray image by GenePix)
for single base pair substitutions in dsDNA sequence targeted by sgRNA
as in Fig. 2B and Fig. S1B

@author: kwbarber
"""

import csv

class obv_bio(object): ##includes normal molecular biology functions
    
    def __init__(self):
        self.temp_str=''
        
    def revC(self, seq):
        self.temp_str=''
        self.rev_str=''
        for i in seq:
            
            if i =='A':
                self.temp_str+='T'
            if i =='T':
                self.temp_str+='A'
            if i =='G':
                self.temp_str+='C'
            if i =='C':
                self.temp_str+='G'
                
            if i =='a':
                self.temp_str+='t'
            if i =='t':
                self.temp_str+='a'
            if i =='g':
                self.temp_str+='c'
            if i =='c':
                self.temp_str+='g'
                
        for i in reversed(self.temp_str):
            self.rev_str+=i
        
afile= ##enter file to be analyzed here

imat=[]
file_reader=open(afile)
csv_reader=csv.reader(file_reader,delimiter=',')
for i in csv_reader:        
    imat.append(i)


barcode='AGATGTGTCGGGCAGTCTCT' #sgRNA spacer of sequence with no substitutions

lmat=['A','C','G','T']

Name=imat[0].index('Name')
Sequence=imat[0].index('Sequence')
Intensity=imat[0].index('F635 Median - B635')

wmat=[]
for i in imat[1:]:
    if 'F' in i[Name] and 'smut' in i[Name]:
        wmat.append([i[Name],i[Sequence][-34:-14],i[Intensity]])


WT_list=[]
for i in imat[1:]:
    if 'F-sub-1' in i[Name]:
        if i[Sequence][-34:-14]==barcode:
            WT_list.append(int(i[Intensity]))

pmat=[]
pcount=0
temp_list=[]
while pcount<=19:
    for i in lmat:
        temp_list=[]
        for j in wmat:
            if j[1][pcount]==i:
                temp_list.append(int(j[2]))
        if len(temp_list)>3:
            pmat.append([pcount+1,i,WT_list])
        else:
            pmat.append([pcount+1,i,temp_list])
    pcount+=1
        
final_mat=[]
for i in pmat:
    temp_sum=0
    for j in i[2]:
        temp_sum+=j
    try:
        final_mat.append([i[0],i[1],temp_sum/len(i[2])]) 
    except ZeroDivisionError:
        final_mat.append([i[0],i[1],'']) 


prism_mat=[[],[],[],[]]
for i in final_mat:   
    if i[1]=='A':
        prism_mat[0].append(i[-1])
    if i[1]=='C':
        prism_mat[1].append(i[-1])
    if i[1]=='G':
        prism_mat[2].append(i[-1])
    if i[1]=='T':
        prism_mat[3].append(i[-1])
    
##replace xxx below with output file location
outWrite=open('xxx.csv'),'w')
counter=0
for i in prism_mat:
    outWrite.write(str(i)+'\n')
outWrite.close()



