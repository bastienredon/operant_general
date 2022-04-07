# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 14:26:17 2021

@author: redon
"""

import os
import glob
import pandas as pd
import numpy as np
from multiprocessing import Pool
import matplotlib.pyplot as plt

import functions_analysis_ethovision as etho

os.chdir('C:\\Users\\redon\\Documents\\Data\\Operant\\reversal')
'''
list_xlsx = glob.glob('*xlsx')

for file in list_xlsx:
    etho.xlsx_to_csv(file)
'''
list_csv = glob.glob('*.csv')


if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        list_data_reversal = p.map(etho.processing, list_csv)
        
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
    data_reversal['Protocol'] = data_reversal['Protocol'].cat.set_categories(['R_FR1.1', 'R_FR1.2', 
                                                                              'R_FR1.3', 'R_FR1.4',
                                                                              'R_FR1.5', 'R_FR1.6'], ordered = True)
    

def firstlight_session(df):
    
    first_light = '-'
    for i in range(len(df.index)):
        if df['Light'].iloc[i] == 1:
             first_light = df['Trial time'].iloc[i]
             break
         
    return first_light

def firstlight_mouse(df_data_mouse):
    
    latency_mouse = 0
    i = 0
    
    
    for value in df_data_mouse['First Light']:
        print(value)
        if not np.isnan(value):
            latency_mouse += value
            break
        else: 
            latency_mouse += df_data_mouse.at[i, 'Duration']
            i += 1
            
            
            
    return latency_mouse

def firstlight_batch(df_data):
    
    df_reversal = pd.DataFrame(columns = [ 'Genotype','Latency_firstlight'])
    df_reversal = df_reversal.astype({'Latency_firstlight' : 'float64'}, copy =False)
    dict_genotype, list_mice = genotype_dict(df_data)
    
    for mouse in list_mice:
        mask_mouse = (df_data['mouse_id'] == mouse)
        latency_mouse = firstlight_mouse(df_data[mask_mouse].reset_index())
        df_reversal.loc[mouse, 'Latency_firstlight'] = latency_mouse
        df_reversal.loc[mouse, 'Genotype'] = dict_genotype[mouse]
        
    return df_reversal
    

def genotype_dict(df_data):
    
    dict_genotype = {}
    
    list_mice = df_data['mouse_id'].unique().tolist()
    index = df_data.index
    
    for mouse in list_mice:
        mouse_mask = (df_data['mouse_id'] == mouse)
        mouse_index = index[mouse_mask].tolist()
        dict_genotype[mouse] = df_data.at[mouse_index[0], 'Genotype']
        
    return dict_genotype, list_mice
        
def plot_reversal(df_reversal, title = 'Reversal learning'):
    
    x_axis_indiv = df_reversal['Genotype']
    y_axis_indiv = df_reversal['Latency_firstlight']
    
    df_means = df_reversal.groupby(['Genotype']).mean()
    
    x_axis_means = df_means.index.tolist()
    y_axis_means = df_means.Latency_firstlight.tolist()
    
    plt.bar(x_axis_means, y_axis_means, alpha = 0.1)
    plt.plot(x_axis_indiv, y_axis_indiv, '.')
    plt.ylabel('Latency (s)')
    plt.xlabel('Genotype')
    plt.title(title)
    