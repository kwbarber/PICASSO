"""
Takes microarray raw data (exported as .txt file from GenePix software) and splits it into 4 subarrays, for CustomArray microarrays using 4 subarrays

@author: kwbarber
"""

import glob
import os
import xlsxwriter
import pandas as pd

file_base= ##direct to folder containing all files to be split into subarrays

if not os.path.exists(file_base+'/_SUBARRAYSPLIT'):
    os.mkdir(os.path.join(file_base,'_SUBARRAYSPLIT'))

for i in glob.glob(os.path.join(file_base, '*.txt')): ##reads in all .txt files from current diretory
    file_reader=open(i)
    subarray_index=''
    headers=''
    row1=1
    row2=1
    row3=1
    row4=1
    print('Reading file '+str(i))

    workbook1 = xlsxwriter.Workbook(str(i).replace('.txt','_S1.xlsx').replace(current_directory,current_directory+'/_SUBARRAYSPLIT'),{'strings_to_numbers': True})
    worksheet1 = workbook1.add_worksheet()
    workbook2 = xlsxwriter.Workbook(str(i).replace('.txt','_S2.xlsx').replace(current_directory,current_directory+'/_SUBARRAYSPLIT'),{'strings_to_numbers': True})
    worksheet2 = workbook2.add_worksheet()
    workbook3 = xlsxwriter.Workbook(str(i).replace('.txt','_S3.xlsx').replace(current_directory,current_directory+'/_SUBARRAYSPLIT'),{'strings_to_numbers': True})
    worksheet3 = workbook3.add_worksheet()
    workbook4 = xlsxwriter.Workbook(str(i).replace('.txt','_S4.xlsx').replace(current_directory,current_directory+'/_SUBARRAYSPLIT'),{'strings_to_numbers': True})
    worksheet4 = workbook4.add_worksheet()

    for j in file_reader:
        if 'Array #' in j:
            subarray_index=j.replace('"','').replace('\n','').split('\t').index('Array #') ##identifies column with subarray identifier
            headers=j.replace('"','').split('\t')

            row=0

            col=0
            for k in headers:
                worksheet1.write(row,col,k)
                col+=1

            col=0
            for k in headers:
                worksheet2.write(row,col,k)
                col+=1

            col=0
            for k in headers:
                worksheet3.write(row,col,k)
                col+=1

            col=0
            for k in headers:
                worksheet4.write(row,col,k)
                col+=1


        if subarray_index!='' and j.replace('\n','').replace('"','').split('\t')[subarray_index]!='-' and 'Array #' not in j:
            subarray_ID=int(j.replace('"','').split('\t')[subarray_index])

            if subarray_ID==1:
                col=0
                for k in j.replace('"','').split('\t'):
                    worksheet1.write(row1,col,k)
                    col+=1
                row1+=1
            if subarray_ID==2:
                col=0
                for k in j.replace('"','').split('\t'):
                    worksheet2.write(row2,col,k)
                    col+=1
                row2+=1
            if subarray_ID==3:
                col=0
                for k in j.replace('"','').split('\t'):
                    worksheet3.write(row3,col,k)
                    col+=1
                row3+=1
            if subarray_ID==4:
                col=0
                for k in j.replace('"','').split('\t'):
                    worksheet4.write(row4,col,k)
                    col+=1
                row4+=1

    workbook1.close()
    workbook2.close()
    workbook3.close()
    workbook4.close()

for i in glob.glob(os.path.join(file_base+'/_SUBARRAYSPLIT/', '*.xlsx')): ##reads in all .xlsx files to replace with .csv
    if '~' not in i:
        data_xls = pd.read_excel(i, 'Sheet1', index_col=None)
        data_xls.to_csv(i.replace('.xlsx','.csv'), encoding='utf-8', index=False)
        os.remove(i)
