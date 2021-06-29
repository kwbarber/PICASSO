"""
Combines microarray spot technical replicates from raw data (.txt) exported from GenePix 

@author: kwbarber
"""

six_file= ##enter .csv file with 635nm laser values here
five_file= ##enter .csv file with 532nm laser values here
four_file= ##enter .csv file with 488nm laser values here

import csv
from statistics import mean
from statistics import stdev

imat_six=[]
file_reader=open(sixfile)
csv_reader=csv.reader(file_reader,delimiter=',')
for i in csv_reader:
    imat_six.append(i)
file_reader.close()

imat_five=[]
file_reader=open(fivefile)
csv_reader=csv.reader(file_reader,delimiter=',')
for i in csv_reader:
    imat_five.append(i)
file_reader.close()

imat_four=[]
file_reader=open(fourfile)
csv_reader=csv.reader(file_reader,delimiter=',')
for i in csv_reader:
    imat_four.append(i)
file_reader.close()

Name6=imat_six[0].index('Name')
Sequence6=imat_six[0].index('Sequence')
Intensity6=imat_six[0].index('F635 Median - B635')
Column6=imat_six[0].index('Column')
Row6=imat_six[0].index('Row')

Name5=imat_five[0].index('Name')
Sequence5=imat_five[0].index('Sequence')
Intensity5=imat_five[0].index('F532 Median - B532')
Column5=imat_five[0].index('Column')
Row5=imat_five[0].index('Row')

Name4=imat_four[0].index('Name')
Sequence4=imat_four[0].index('Sequence')
Intensity4=imat_four[0].index('F488 Median - B488')
Column4=imat_four[0].index('Column')
Row4=imat_four[0].index('Row')

class obv_bio(object): ##includes normal molecular biology functions

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

    def rev_transl(self, seq):
        ##for single sequence input, this attribute relies on "restr_avoid" attribute as well
        ##for lists, performs codon optimization without regard to restr site presence
        ##individual sequences will be codon optimized and checked for restriction site presence iteratively

        ##non-rare codon dictionary
        self.codon_dict = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'C': ['TGT', 'TGC'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGC', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC'],
            'K': ['AAA', 'AAG'],
            'L': ['TTG', 'CTT', 'CTC', 'CTG', 'TTA'],
            'M': ['ATG'],
            'N': ['AAT', 'AAC'],
            'P': ['CCT', 'CCA', 'CCG'],
            'Q': ['CAA', 'CAG'],
            'R': ['CGC', 'CGT'],
            'S': ['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'V': ['GTT', 'GTC', 'GTG', 'GTA'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            '*': ['TGA', 'TAA', 'TAG'],
            's': ['TAG'],  ##these are for amber suppression for phospho incorporation
            't': ['TAG'],
            'y': ['TAG'],
            'X': ['XXX']
        }

        self.temp_counter = 0
        self.temp_str = ''
        self.peptide_seq = ''
        for i in seq:
            self.temp_str += i
            if len(self.temp_str) == 3:
                for j in self.codon_dict:
                    if self.temp_str in self.codon_dict[j]:
                        self.peptide_seq += j
                self.temp_str = ''

    def translate(self, seq):
        self.codon_dict = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'C': ['TGT', 'TGC'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGC', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC'],
            'K': ['AAA', 'AAG'],
            'L': ['TTG', 'CTT', 'CTC', 'CTG', 'TTA'],
            'M': ['ATG'],
            'N': ['AAT', 'AAC'],
            'P': ['CCT', 'CCA', 'CCG'],
            'Q': ['CAA', 'CAG'],
            'R': ['CGC', 'CGT'],
            'S': ['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'V': ['GTT', 'GTC', 'GTG', 'GTA'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            '*': ['TGA', 'TAA', 'TAG'],
            's': ['TAG'],  ##these are for amber suppression for phospho incorporation
            't': ['TAG'],
            'y': ['TAG'],
            'X': ['XXX']
        }

        self.temp_seq=str(seq)
        self.seq_len=int(len(seq)/3)
        self.pep_seq=''
        for i in range(self.seq_len):
            for j in self.codon_dict:
                if self.temp_seq[i*3:(i*3)+3] in self.codon_dict[j]:
                    self.pep_seq+=j


lmat=['A','C','G','T']

wmat_six=[]
for i in imat_six[1:]:
    wmat_six.append([i[Name6],i[Sequence6][-34:-14],i[Intensity6]])

wmat_five=[]
for i in imat_five[1:]:
    wmat_five.append([i[Name5],i[Sequence5][-34:-14],i[Intensity5]])

wmat_four=[]
for i in imat_four[1:]:
    wmat_four.append([i[Name4],i[Sequence4][-34:-14],i[Intensity4]])

int_dict={}
for i in wmat_four:
    if i[0] not in int_dict:
        int_dict[i[0]] = [[float(i[-1])],[],[]]
    else:
        int_dict[i[0]][0].append(float(i[-1]))

for i in wmat_five:
    int_dict[i[0]][1].append(float(i[-1]))

for i in wmat_six:
    int_dict[i[0]][2].append(float(i[-1]))

final_mat=[['Seq','I488','SD488','I532','SD532','I635','SD635','I488rawval','I532rawval','I635rawval']]
for i in int_dict:
    final_mat.append([i,mean(int_dict[i][0]),stdev(int_dict[i][0]),mean(int_dict[i][1]),stdev(int_dict[i][1]),mean(int_dict[i][2]),stdev(int_dict[i][2]),str(int_dict[i][0]).replace(',','|'),str(int_dict[i][1]).replace(',','|'),str(int_dict[i][2]).replace(',','|')])


##replace 'xxx' with output file:
outWrite=open('xxx.csv','w')

counter=0
for i in final_mat:
    outWrite.write(str(i).replace(']','').replace('[','').replace(' ','').replace("'","")+'\n')
outWrite.close()
