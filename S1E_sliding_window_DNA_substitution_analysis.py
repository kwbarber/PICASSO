"""
File to analyze raw data (exported .txt file from microarray image by GenePix)
performing sliding window dsDNA target substitution analysis, as in Fig. S1E
"""

import csv
from statistics import mean
from statistics import stdev

barcode='CTCACGGAGAAATCTATGGA' ##sgRNA spacer for analysis
afile = ##enter file to be analyzed

class obv_bio(object): ##normal molecular biology software functions

    def __init__(self):
        self.temp_str = ''

    def revC(self, seq):
        self.temp_str = ''
        self.rev_str = ''
        for i in seq:

            if i == 'A':
                self.temp_str += 'T'
            if i == 'T':
                self.temp_str += 'A'
            if i == 'G':
                self.temp_str += 'C'
            if i == 'C':
                self.temp_str += 'G'

            if i == 'a':
                self.temp_str += 't'
            if i == 't':
                self.temp_str += 'a'
            if i == 'g':
                self.temp_str += 'c'
            if i == 'c':
                self.temp_str += 'g'

        for i in reversed(self.temp_str):
            self.rev_str += i

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

max_list=[]
for i in imat[1:]:
    if 'F' in i[Name] and 'sub' in i[Name]:
        max_list.append(int(i[Intensity])-norm_factor) ##normalized using NP sequence
max_val=max(max_list)

wmat=[]
for i in imat[1:]:
    if 'F' in i[Name] and 'sub' in i[Name]:
        wmat.append([i[Name],i[Sequence][-34:-14],str((float(i[Intensity])-norm_factor)/max_val)]) ##normalization factor subtracted

seq_dict={}
for i in wmat:
    seq_dict[i[0]]=[]
for i in wmat:
    if float(i[-1])<0:
        seq_dict[i[0]].append(0)
    else:
        seq_dict[i[0]].append(float(i[-1]))

final_mat=[]
for i in seq_dict:
    final_mat.append([i,str(seq_dict[i]).replace(',','|')])

##replace xxx below with output file path
outWrite=open('xxx.txt','w')
counter=0
for i in reversed(final_mat):
    outWrite.write(str(i[1:]).replace(',','').replace('[','').replace(']','').replace('|','\t').replace("'","")+'\n')
outWrite.close()
