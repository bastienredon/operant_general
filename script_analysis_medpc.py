# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 16:24:15 2021

@author: redon
"""

import pandas as pd
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import analysis_medassociate as med #module analysis medassociate 

os.chdir('C:\\Users\\redon\\Documents\\Data\\Operant\\extinction')

#-----------------------------------------------------------------------------

#                   PART 1: EXTRACT AND ORGANISE DATA

#-----------------------------------------------------------------------------

#Extract all extinction files:
list_experiment = glob.glob('*.txt')

#Create one dictionary per session to host experimental data (1 dict/mouse/session)
dict_EXT1 = {}
dict_EXT2 = {}
dict_EXT3 = {}
dict_EXT4 = {}
dict_EXT5 = {}
dict_EXT6 = {}
dict_EXT7 = {}
dict_EXT8 = {}

#Initiate the DataFrame hosting discrete data (number of NP, genotype, mouse id and protocol)
df_final = pd.DataFrame(columns = ['Protocol', 'Stimulus', 'Genotype', 'ID','Door_Open','NPA_total', 'NPI','NPI_total'])

#Iterate through each file to extract data
for file in list_experiment:
    
    #Extract for a given file all experimental information
    df_extinction, dict_result, genotype, protocol, mouseid, extinction_type = med.generate_general_df(file)
   
    #Append the suitable information to the general dataframe
    df_final = df_final.append(df_extinction, ignore_index = True)
    dict_result['Protocol'] = protocol
    dict_result['Genotype'] = genotype
    dict_result['Stimulus'] = extinction_type
    key = str(mouseid+'_'+protocol)
    
    #Then add the data dictionary to the suitable session dictionary as value (key = mouse_id)
    if protocol == 'EXT1':
       
        dict_EXT1[key] = dict_result
        
    elif protocol == 'EXT2':
        
        dict_EXT2[key] = dict_result
        
    elif protocol == 'EXT3':
        
        dict_EXT3[key] = dict_result
        
    elif protocol == 'EXT4':
        
        dict_EXT4[key] = dict_result
        
    elif protocol == 'EXT5':
        
        dict_EXT5[key] = dict_result
    
    elif protocol == 'EXT6':
        
        dict_EXT6[key] = dict_result
        
    elif protocol == 'EXT7':
        
        dict_EXT7[key] = dict_result
        
    elif protocol == 'EXT8':
        
        dict_EXT8[key] = dict_result
        
list_dict_ext = [dict_EXT1, dict_EXT2, dict_EXT3, dict_EXT4, dict_EXT5, dict_EXT6, dict_EXT7, dict_EXT8]

dict_complete = {**dict_EXT1, **dict_EXT2, **dict_EXT3, **dict_EXT4, **dict_EXT5, **dict_EXT6, **dict_EXT7, **dict_EXT8}
        
        
#-----------------------------------------------------------------------------

#                   PART 2: PLOT AND ANALYSE MEAN DATA

#-----------------------------------------------------------------------------       
'''

#Subset suitable line of df_final according to respective genotype       
(df_wt_social, df_wt_nonsocial, 
 df_he_social, df_he_nonsocial, 
 df_ho_social, df_ho_nonsocial) = med.segregate_genotypes(df_final)

#Plot NPA of both wt and ho mice in extinction with a social stimuli
med.plot_npa_extinction(df_wt_social, df_ho_social, 'Extinction_social')

#Plot NPA of both wt and ho mice in extinction without a social stimuli
med.plot_npa_extinction(df_wt_nonsocial, df_ho_nonsocial, 'Extinctionnosocial')


#-----------------------------------------------------------------------------

#                   PART 3: PLOT AND ANALYSE TIMECOURSE DATA

#-----------------------------------------------------------------------------

#Plot timecourse in rasterplots (1/session with all subjects)
for dictionary in list_dict_ext:
    med.plot_timecourse_npa_persession(dictionary)

fig_binned = 0
for dictionary in list_dict_ext:    
    med.plot_npa_timecourse_binned(dictionary, fig_binned)
    fig_binned +=1
   
fig = 0   
for dictionary in list_dict_ext:
    med.plot_npa_timecourse_cumul(dictionary, fig)
    fig +=1
    
med.plot_timecourse_npa_pergenotype(dict_complete)
'''