# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 17:10:42 2021

@author: Bastien

Main script for ethovision file analysis

Use the module analysis_ethovision_shank3opto.py

"""

import functions_analysis_ethovision as etho

import glob 
import pandas as pd
from multiprocessing import Pool
import os
#-----------------------------------------------------------------------------

#                            CONVERTING TO CSV

#-----------------------------------------------------------------------------

#Set the working directory
os.chdir('C:\\Users\\redon\\Documents\\Data\\Operant\\OC_DatShank3_M_Autoopto_June21\\Raw_data_202109_DatShank_Operant')
'''
#Convert all xlsx files into csv
list_file = glob.glob("*.xlsx")
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        p.map(etho.xlsx_to_csv, list_file)
'''
#-----------------------------------------------------------------------------

#                            EXTRACTING DATA

#-----------------------------------------------------------------------------

#Set the working directory
#os.chdir('C:\\Users\\redon\\Documents\\Data\\operant\\OC_DatShank3_M_Autoopto_June21\\Ethovision_raw_45min\\csv')

list_csv = glob.glob('*.csv')

#Classify csv files based on their protocol
list_file_SH, list_file_FR1, list_file_FR3, list_file_PR, list_file_reversal = etho.csv_classify_protocol(list_csv)

#Generate SH dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_SH = p.map(etho.processing, list_file_SH)
        
    data_SH = pd.concat(list_data_SH, axis = 0, ignore_index = True)
    data_SH = data_SH.astype({'Protocol': "category",
                                'Duration': 'float64',
                                'Genotype': "category",
                                'mouse_id': 'category',
                                'NPA': 'float64', 'NPI': 'float64', 
                                'Door Open': 'float64', 'Light': 'float64', 
                                'Mean ITI': 'float64', 'Total time light' : 'float64', 
                                'Mean time light' : 'float64', "DO IZ": "float64", 
                                "DC IZ" : 'float64', 'Time Zone NPA': 'float64', 
                                'First NPA': 'float64', 'First Light': 'float64'}, copy =False)


    #Save the dataframe as csv
    data_SH.to_csv('df\\data\\df_data_SH.csv')
    
#Generate FR1 dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_FR1 = p.map(etho.processing, list_file_FR1)
        
    data_FR1 = pd.concat(list_data_FR1, axis = 0, ignore_index = True)
    data_FR1 = data_FR1.astype({'Protocol': "category",
                                'Duration': 'float64',
                                'Genotype': "category",
                                'mouse_id': 'category',
                                'NPA': 'float64', 'NPI': 'float64', 
                                'Door Open': 'float64', 'Light': 'float64', 
                                'Mean ITI': 'float64', 'Total time light' : 'float64', 
                                'Mean time light' : 'float64', "DO IZ": "float64", 
                                "DC IZ" : 'float64', 'Time Zone NPA': 'float64', 
                                'First NPA': 'float64', 'First Light': 'float64'}, copy =False)
    
    data_FR1['Protocol'] = data_FR1['Protocol'].cat.set_categories(['FR1.1','FR1.2','FR1.3',
                                                            'FR1.4','FR1.5','FR1.6',
                                                            'FR1.7','FR1.8','FR1.9',
                                                            'FR1.10','FR1.11',
                                                            'FR3.1','FR3.2','FR3.3',
                                                            'FR3.4','FR3.5','FR3.6',
                                                            'FR3.7','FR3.8','FR3.9',
                                                            'FR3.10','FR3.11', 'FR3.12',
                                                            'FR3.13','FR3.14',
                                                            'FR3.15','FR3.16'], ordered = True)

    #Save the dataframe as csv
    data_FR1.to_csv('df\\data\\df_data_FR1.csv')

#Generate FR3 dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_FR3 = p.map(etho.processing, list_file_FR3)
        
    data_FR3 = pd.concat(list_data_FR3, axis = 0, ignore_index = True)
    data_FR3 = data_FR3.astype({'Protocol': "category",
                                'Duration': 'float64',
                                'Genotype': "category",
                                'mouse_id': 'category',
                                'NPA': 'float64', 'NPI': 'float64', 
                                'Door Open': 'float64', 'Light': 'float64', 
                                'Mean ITI': 'float64', 'Total time light' : 'float64', 
                                'Mean time light' : 'float64', "DO IZ": "float64", 
                                "DC IZ" : 'float64', 'Time Zone NPA': 'float64', 
                                'First NPA': 'float64', 'First Light': 'float64'}, copy =False)
    data_FR3['Protocol'] = data_FR3['Protocol'].cat.set_categories(['FR1.1','FR1.2','FR1.3',
                                                            'FR1.4','FR1.5','FR1.6',
                                                            'FR1.7','FR1.8','FR1.9',
                                                            'FR1.10','FR1.11',
                                                            'FR3.1','FR3.2','FR3.3',
                                                            'FR3.4','FR3.5','FR3.6',
                                                            'FR3.7','FR3.8','FR3.9',
                                                            'FR3.10','FR3.11', 'FR3.12',
                                                            'FR3.13','FR3.14',
                                                            'FR3.15','FR3.16'], ordered = True)
    
    #Save the dataframe as csv
    data_FR3.to_csv('df\\data\\df_data_FR3.csv')


#Generate PR dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_PR = p.map(etho.processing, list_file_PR)
        
    data_PR = pd.concat(list_data_PR, axis = 0, ignore_index = True)
    data_PR = data_PR.astype({'Protocol': "category",
                              'Duration': 'float64',
                              'Genotype': "category",
                              'mouse_id': 'category',
                              'NPA': 'float64', 'NPI': 'float64', 
                              'Door Open': 'float64', 'Light': 'float64', 
                              'Mean ITI': 'float64', 'Total time light' : 'float64', 
                              'Mean time light' : 'float64', "DO IZ": "float64", 
                              "DC IZ" : 'float64', 'Time Zone NPA': 'float64', 
                              'First NPA': 'float64', 'First Light': 'float64'}, copy =False)
   

    #Save the dataframe as csv
    data_PR.to_csv('df\\data\\df_data_PR.csv')


#Generate reversal dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_reversal = p.map(etho.processing, list_file_reversal)
        
    data_reversal = pd.concat(list_data_reversal, axis = 0, ignore_index = True)
    data_reversal = data_reversal.astype({'Protocol': "category",
                              'Duration': 'float64',
                              'Genotype': "category",
                              'mouse_id': 'category',
                              'NPA': 'float64', 'NPI': 'float64', 
                              'Door Open': 'float64', 'Light': 'float64', 
                              'Mean ITI': 'float64', 'Total time light' : 'float64', 
                              'Mean time light' : 'float64', "DO IZ": "float64", 
                              "DC IZ" : 'float64', 'Time Zone NPA': 'float64', 
                              'First NPA': 'float64', 'First Light': 'float64'}, copy =False)
   

    #Save the dataframe as csv
    data_reversal.to_csv('df\\data\\df_data_reversal.csv')

list_files = [list_file_SH, list_file_FR1, list_file_FR3, list_file_PR, list_file_reversal]

list_name = ['list_file_SH', 'list_file_FR1', 'list_file_FR3', 'list_file_PR', 'list_file_reversal']

i = 0
for file in list_files:
    name = list_name[i]
    df_file = pd.DataFrame(file)
    df_file.to_csv('df\\files\\'+ name + '.csv')
    i+=1