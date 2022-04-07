# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:23:38 2021
@author: redon

Main script for behavioral analysis
"""
import os 
import glob
import pandas as pd
import extracting_data_ethovision 
import functions_analysis_ethovision as etho

new_analysis = input('Start new analysis: press Y, if not press N \n')

#-----------------------------------------------------------------------------

#                   Extract the data or open data stored

#-----------------------------------------------------------------------------
if new_analysis == 'Y':
    print('Extraction script is being executed...\n')
    os.chdir('C:\\Users\\redon\\Documents\\code')
    %run extracting_data_ethovision.py
    print('Extraction completed')
    
elif new_analysis == 'N':
 
    print('Opening 45min-session data in progress...\n')
    os.chdir('C:\\Users\\redon\\Documents\\Data\\Operant\\OC_DatShank3_M_Autoopto_June21\\Raw_data_202109_DatShank_Operant')
    df_data_SH = pd.read_csv('df\\data\\df_data_SH.csv', index_col = [0]) 
    df_data_FR1 = pd.read_csv('df\\data\\df_data_FR1.csv', index_col = [0]) 
    df_data_FR3 = pd.read_csv('df\\data\\df_data_FR3.csv', index_col = [0]) 
    df_data_PR = pd.read_csv('df\\data\\df_data_PR.csv', index_col = [0]) 
    df_data_reversal = pd.read_csv('df\\data\\df_data_reversal.csv', index_col = [0]) 
    print('Done!\n')
    
    print('Opening 45min-session lists in progress...\n') 
    list_file_SH = pd.read_csv('df\\files\\list_file_SH.csv', index_col = [0], header  = 0, names = ['Files'])
    list_file_SH = list(list_file_SH.iloc[:,0])
    list_file_FR1 = pd.read_csv('df\\files\\list_file_FR1.csv', index_col = [0], header  = 0, names = ['Files']) 
    list_file_FR1 = list(list_file_FR1.iloc[:,0])
    list_file_FR3 = pd.read_csv('df\\files\\list_file_FR3.csv', index_col = [0], header  = 0, names = ['Files']) 
    list_file_FR3 = list(list_file_FR3.iloc[:,0])
    list_file_PR = pd.read_csv('df\\files\\list_file_PR.csv', index_col = [0], header  = 0, names = ['Files']) 
    list_file_PR = list(list_file_PR.iloc[:,0])
    list_file_reversal = pd.read_csv('df\\files\\list_file_reversal.csv', index_col = [0], header  = 0, names = ['Files']) 
    list_file_reversal = list(list_file_reversal.iloc[:,0])
    print('Done!\n')
    
else:
    print('What?')
  


#-----------------------------------------------------------------------------

#                             Plot timecourse NPA 

#-----------------------------------------------------------------------------

#Plot 45min session data
os.chdir('C:\\Users\\redon\\Documents\\Data\\Operant\\OC_DatShank3_M_Autoopto_June21\\Raw_data_202109_DatShank_Operant')

bins_20min = range(2,22,2)
bins_45min = range(5,55,5)

#FR1
df_npa_wt_FR1, df_npa_he_FR1, df_npa_ho_FR1 = etho.timecourse_npa(list_file_FR1)
#Calculate binned stats for each genotype
df_stats_wt_FR1 = etho.timecourse_npa_binned(df_npa_wt_FR1, bins_20min)
df_stats_he_FR1 = etho.timecourse_npa_binned(df_npa_he_FR1, bins_20min)
df_stats_ho_FR1 = etho.timecourse_npa_binned(df_npa_ho_FR1, bins_20min)
#Plot each binned stats into a single graph
etho.plot_timecourse_npa(df_stats_wt_FR1, df_stats_he_FR1, df_stats_ho_FR1, 'FR1', bins_20min, y_top = 50)

#FR3
df_npa_wt_FR3, df_npa_he_FR3, df_npa_ho_FR3 = etho.timecourse_npa(list_file_FR3)
#Calculate binned stats for each genotype
df_stats_wt_FR3 = etho.timecourse_npa_binned(df_npa_wt_FR3, bins_20min)
df_stats_he_FR3 = etho.timecourse_npa_binned(df_npa_he_FR3, bins_20min)
df_stats_ho_FR3 = etho.timecourse_npa_binned(df_npa_ho_FR3, bins_20min)
#Plot each binned stats into a single graph
etho.plot_timecourse_npa(df_stats_wt_FR3, df_stats_he_FR3, df_stats_ho_FR3, 'FR3', bins_20min, y_top = 50)


#FR3 post PR
df_npa_wt_rever, df_npa_he_rever, df_npa_ho_rever = etho.timecourse_npa(list_file_reversal)
#Calculate binned stats for each genotype
df_stats_wt_rever = etho.timecourse_npa_binned(df_npa_wt_rever, bins_45min)
df_stats_he_rever = etho.timecourse_npa_binned(df_npa_he_rever, bins_45min)
df_stats_ho_rever = etho.timecourse_npa_binned(df_npa_ho_rever, bins_45min)
#Plot each binned stats into a single graph
etho.plot_timecourse_npa(df_stats_wt_rever, df_stats_he_rever, df_stats_ho_rever, 'Reversal Learning', bins_45min)

