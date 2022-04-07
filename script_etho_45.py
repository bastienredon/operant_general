# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:38:47 2021

Main script for Ethovision analysis for 45 min sessions

@author: redon
"""

import analysis_ethovision_shank3opto as etho
import glob 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
import re
from multiprocessing import Pool

#Set the working directory

os.chdir('C:\\Users\\redon\\Desktop\\Ethovision')

#Convert all xlsx files into csv
list_file = glob.glob("*.xlsx")
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        p.map(etho.xlsx_to_csv, list_file)
     
list_csv = glob.glob('*.csv')

#Classify csv files based on their protocol
list_file_SH, list_file_FR1, list_file_FR3, list_file_PR = etho.csv_classify_protocol(list_csv)


#-----------------------------------------------------------------------------

#                            ANALYSE - 20 min

#-----------------------------------------------------------------------------

#Generate FR1 dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_FR1 = p.map(etho.processing, list_file_FR1)
        
    data_FR1 = pd.concat(list_data_FR1, axis = 0, ignore_index = True)
    data_FR1 = data_FR1.astype({'Protocol': "category",
                        'Genotype': "category",
                        'mouse_id': 'category',
                        'NPA': 'float64', 'NPI': 'float64', 
                        'Door Open': 'float64', 'Light': 'float64', 
                        'Mean ITI': 'float64', 'Total time light' : 'float64', 
                        'Mean time light' : 'float64', "DO IZ": "float64", 
                        "DC IZ" : 'float64', 'Time Zone NPA': 'float64'}, copy =False)
    data_FR1['Protocol'] = data_FR1['Protocol'].cat.set_categories(['FR1.1','FR1.2','FR1.3',
                                                            'FR1.4','FR1.5','FR1.6',
                                                            'FR1.7','FR1.8','FR1.9',
                                                            'FR1.10','FR1.11',
                                                            'FR3.1','FR3.2','FR3.3',
                                                            'FR3.4','FR3.5','FR3.6',
                                                            'FR3.7','FR3.8','FR3.9',
                                                            'FR3.10','FR3.11'], ordered = True)

#Generate FR3 dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_FR3 = p.map(etho.processing, list_file_FR3)
        
    data_FR3 = pd.concat(list_data_FR3, axis = 0, ignore_index = True)
    data_FR3 = data_FR3.astype({'Protocol': "category",
                        'Genotype': "category",
                        'mouse_id': 'category',
                        'NPA': 'float64', 'NPI': 'float64', 
                        'Door Open': 'float64', 'Light': 'float64', 
                        'Mean ITI': 'float64', 'Total time light' : 'float64', 
                        'Mean time light' : 'float64', "DO IZ": "float64", 
                        "DC IZ" : 'float64', 'Time Zone NPA': 'float64'}, copy =False)
    data_FR3['Protocol'] = data_FR3['Protocol'].cat.set_categories(['FR1.1','FR1.2','FR1.3',
                                                            'FR1.4','FR1.5','FR1.6',
                                                            'FR1.7','FR1.8','FR1.9',
                                                            'FR1.10','FR1.11',
                                                            'FR3.1','FR3.2','FR3.3',
                                                            'FR3.4','FR3.5','FR3.6',
                                                            'FR3.7','FR3.8','FR3.9',
                                                            'FR3.10','FR3.11'], ordered = True)
#Generate PR dataframe with discrete values
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_PR = p.map(etho.processing, list_file_PR)
        
    data_PR = pd.concat(list_data_PR, axis = 0, ignore_index = True)
    data_PR = data_PR.astype({'Protocol': "category",
                        'Genotype': "category",
                        'mouse_id': 'category',
                        'NPA': 'float64', 'NPI': 'float64', 
                        'Door Open': 'float64', 'Light': 'float64', 
                        'Mean ITI': 'float64', 'Total time light' : 'float64', 
                        'Mean time light' : 'float64', "DO IZ": "float64", 
                        "DC IZ" : 'float64', 'Time Zone NPA': 'float64'}, copy =False)
    
    
    
