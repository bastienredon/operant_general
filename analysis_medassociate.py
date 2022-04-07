# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 10:43:56 2021

@author: redon

Analysis of MedAssociate files (text files non organized)
 Variables:
\ A = fixed ratio
\ B = session duration max
\ C = nb np max in session
\ D = timer for session
\ E = timer for table
\ F = count noldus
\ G = array noldus
\ H = idx array noldus
\ I = count inactive
\ J = array inactive
\ K = idx array inactive
\ L = count inactive all
\ M = array inactive all
\ N = idx array inactive all
\ O = temporary count active all
\ P = count active
\ Q = array active
\ R = idx array active
\ S = count active all
\ T = array active all
\ U = idx array active all
"""

import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
from scipy.stats import sem as scp_sem

os.chdir('C:\\Users\\redon\\Documents\\Data\\Operant\\extinction')
file = '2021-08-03_18h38m_Subject 4235-88.txt'

#-----------------------------------------------------------------------------


#FOLLOWING FUNCTIONS ARE MEANT TO EXTRACT, COMPUTE AND ORGANIZE DATA


#-----------------------------------------------------------------------------
 
def detect_floats(sequence):
    '''Function computing medassociate file data by detecting list of numbers in non-organized txt file

    Arg:
        
        sequence = chunk of non organized txt file containing several float to extract
            type Str

    Returns

        list_floats = List containing each single floats detected in the sequence
            type List of floats

    '''
    #Pattern to match any decimal number whatever the size
    pattern = re.compile(r'[0-9]+\.?[0-9]+') #any several number ([0-9]+) followed by a period or not(\.?), and followed by several number
    
    #Find each single decimal number in the sequence
    list_floats_str = re.findall(pattern, sequence)
    
    #Convert the list of string into list of floats to allow future computation
    list_floats = [float(x) for x in list_floats_str]
    
    return list_floats

def import_txt(file):
    """Function allowing to (i) import medpc file into a dataframe 
    and to (ii) convert columns of data into arrays of floats for computation
    WARNING: use detect_floats()
    
    Arg:
        
        file = path or name of the file to open, must be txt or csv-like
            type Str
            
    Return:
        
        array_results = array version of df_medpc_all
            type Array 
        
        df_medpc_all = 2 columns dataframe containing experimental data
            type DataFrame
        
    """

    df_medpc_all = pd.read_csv(file, 
                               sep = ':', #Allow to have 2 column but will raise error on dates
                               skiprows=[0], 
                               on_bad_lines = 'warn', #Allow to overcome error with dates by simply warning when one is encountered
                               header = None).fillna('empty')

    #Converting the dataframe into array
    array_results = df_medpc_all.values
    
    for i in range(len(array_results)):
        array_results[i,1] = detect_floats(array_results[i,1])
    
    return array_results, df_medpc_all

def extract_arrays(array_result):
    """Function iterating through arrays to extract and organize the all data 
    into a dictionary with the variable's name (key) and its content (value)
    
    Arg:
        
        array_results = array containing data from the dataframe 
        
    Return:
        
        dict_result = dictionnary gathering variable names (keys) and values in an array (value) 
            type Dict
        
        list_variables = list of variable names a string 
            type List of Str
        
        list_values = list of variable values in list
            type List of lists
    """
    
    #Initiate empty lists to host data
    list_variables = []
    list_values = []
    
    #Define a regex that would match any uppercase letter
    letter = re.compile(r'[A-Z]')
    
    #Iterate through the row of the array:
    for j in range(len(array_result)-1):
        
        #If column 0 element at index j and j+1 are capital letters
        if re.match(letter, array_result[j,0]) != None  and re.match(letter, array_result[j+1,0]) != None:
            #Add the variable name (element of index j in column 0) to list_variables 
            list_variables.append(re.match(letter, array_result[j,0]).group())
            #Add its value (element of index j in column 1) to list_values
            list_values.append(array_result[j,1])
        
        #If column 0 element is a capital letter and the following is not 
        elif re.match(letter, array_result[j,0]) != None  and re.match(letter, array_result[j+1,0]) == None:
            #Add the variable name (element of index j in column 0) to list_variables
            list_variables.append(re.match(letter, array_result[j,0]).group())
            #Add its value (element of index j in column 1) to list_values to initiate the array
            list_values.append(array_result[j,1])
         
        #If column 0 element is not a capital letter    
        elif re.match(letter, array_result[j,0]) == None:
            #Append its value (element of index j in column 1) to the last element of list_values
            list_values[-1] = list_values[-1] + array_result[j,1]
                
    #Create the dictionnary with all extracted data: name as keys and values as values
    zip_iterator = zip(list_variables, list_values)
    dict_result = dict(zip_iterator)
    
    #Iterate through element within the dictionary:
    for key, value in dict_result.items():
        #If no value is found, set to nan 
        if len(value) == 0:
            dict_result[key] = 0
            
        #if only one element is present in the value, extract it from the value list 
        elif len(value) == 1:
            dict_result[key] = value[0]
            
    return dict_result, list_variables, list_values
   

def generate_general_df(file):
    '''Function reading and extracting data from medassociate file and organize
    these data into a single DataFrame
    
    Arg:
        file = name or path for a .txt file to analyse
            type Str or Path
            
    Return: 
        df_extinction = General df with discrete value of protocol, genotype, ID, NPA and NPI
            type DataFrame
            
        dict_result = Dictionary gathering for each mouse (key) all its data (value)
            type Dict
            
        genotype = Genotype of the subject
            type Str
            
        protocol = Name of the session analysed
            type Str
            
        mouse_id = Identification of the experimental subject
            type Str
    '''
    #Open and read the raw file
    df_info = pd.read_csv(file, 
                          sep = ':', 
                          skiprows = [0],  
                          on_bad_lines = 'warn', 
                          nrows = 25,
                          header = None, 
                          index_col = [0]).fillna('0')
    
    #Create a dataframe meant to host the file's data
    df_extinction = pd.DataFrame(columns = ['Protocol',
                                            'Stimulus',
                                            'Genotype', 
                                            'ID',
                                            'Door_Open',
                                            'NPA_total', 
                                            'NPI',
                                            'NPI_total'])
    
    #Create a regex to capture genotype in the file 
    regex_extinction_type = re.compile(r'\w{3,9}') 
    extinction_type = re.search(regex_extinction_type, df_info.loc['Group',1]).group()
    
    #Create a regex to capture genotype in the file 
    regex_genotype = re.compile(r'\w{2}')
    genotype = re.search(regex_genotype, df_info.loc['Group',1]).group()
    
    #Create a regex to capture genotype in the file 
    regex_mouseid = re.compile(r'\d{4}.\d{,2}')
    mouseid = re.search(regex_mouseid, df_info.loc['Subject',1]).group()
    
    #Identify the session analysed
    regex_session = re.compile(r'\w{3}\d{1}')
    protocol = re.search(regex_session, df_info.loc['Experiment',1]).group()
    
    #Add the variable values to the df
    df_extinction.loc[0, 'Protocol'] = protocol
    df_extinction.loc[0, 'Stimulus'] = extinction_type
    df_extinction.loc[0, 'Genotype'] =  genotype
    df_extinction.loc[0, 'ID'] = mouseid
    df_extinction.NPA_total =  df_info.loc['S'].item()
    df_extinction.NPI = df_info.loc['I'].item()
    df_extinction.NPI_total = df_info.loc['L'].item()
    df_extinction.Door_Open =  df_info.loc['P'].item()
    
    #create a list of array
    array_results, df_medpc_all = import_txt(file)
    dict_result, list_variables, list_values = extract_arrays(array_results)
    
    #Define the data type within the df to facilitate future computations
    df_extinction = df_extinction.astype({'Protocol' : 'category', 
                                          'Stimulus': 'string',
                                          'Genotype' : 'string',
                                          'ID' : 'string',
                                          'Door_Open' : 'float64',
                                          'NPA_total' : 'float64',
                                          'NPI' : 'float64',
                                          'NPI_total' : 'float64'})
    
    
    return df_extinction, dict_result, genotype, protocol, mouseid, extinction_type
 


#-----------------------------------------------------------------------------


#FOLLOWING FUNCTIONS ARE MEANT TO WORK ON DISCRETE MEAN VALUES


#-----------------------------------------------------------------------------

def segregate_genotypes(df):
    
    #Create mask to subset genotypes
    wt = (df.Genotype == 'wt')
    he = (df.Genotype == 'he')
    ho = (df.Genotype == 'ho')
    
    
    social = (df.Stimulus == 'social')
    nonsocial = (df.Stimulus == 'nonsocial')
    
    
    
    
    df_wt_social = df[social & wt]
    df_wt_nonsocial = df[nonsocial & wt]
    df_he_social = df[social & he]
    df_he_nonsocial = df[nonsocial & he]
    df_ho_social = df[social & ho]
    df_ho_nonsocial = df[nonsocial & ho]
    
    return df_wt_social, df_wt_nonsocial, df_he_social, df_he_nonsocial, df_ho_social, df_ho_nonsocial
    

def plot_npa_extinction(df_wt, df_ho, title):
    
    plt.figure(dpi = 800)
    
    #Plot individual values
    plt.plot(df_wt.Protocol,
             df_wt.NPA, 'b.',
             label = 'Dat-cre$^{het}$,Shank3$^{wt}$')
    plt.plot(df_ho.Protocol,
             df_ho.NPA, 'r.',
             label = 'Dat-cre$^{het}$,Shank3$^{ho}$')
    
    #Compute mean wt
    mean_values_wt = df_wt.groupby(['Protocol']).mean()
    sem_values_wt = df_wt.groupby(['Protocol']).sem()
    
    #Compute mean ho
    mean_values_ho = df_ho.groupby(['Protocol']).mean()
    sem_values_ho = df_ho.groupby(['Protocol']).sem()
    
    #Plot mean
    plt.plot(mean_values_wt.index, 
             mean_values_wt.NPA, 'b-',
             mean_values_ho.NPA, 'r-',
             linewidth = 0.7)
    '''
    #Plot SEM
    plt.plot(mean_values_wt.index, 
             mean_values_wt.NPA - sem_values_wt.NPA, 'b--',
             mean_values_ho.NPA - sem_values_ho.NPA, 'r--',
             linewidth = 0.4)
    plt.plot(mean_values_wt.index, 
             mean_values_wt.NPA + sem_values_wt.NPA,'b--',
             mean_values_ho.NPA + sem_values_ho.NPA,'r--',
             linewidth = 0.4)
    
    plt.fill_between(mean_values_wt.index,
                     mean_values_wt.NPA - sem_values_wt.NPA,
                     mean_values_wt.NPA + sem_values_wt.NPA,
                     color = 'b',
                     alpha = 0.3)
    
    plt.fill_between(mean_values_ho.index,
                     mean_values_ho.NPA - sem_values_ho.NPA,
                     mean_values_ho.NPA + sem_values_ho.NPA,
                     color = 'r',
                     alpha = 0.3)
    '''
    plt.ylim(0, 500)
    plt.ylabel('Number of nose pokes')
    plt.title(title)
    plt.legend()
    plt.savefig(title)
 
#-----------------------------------------------------------------------------


#FOLLOWING FUNCTIONS ARE MEANT TO WORK ON TIMECOURSE


#-----------------------------------------------------------------------------   
 
def plot_timecourse_npa_persession(dict_EXTx):
    
    
    i = 0
    
    #Initiate the figure
    plt.figure(dpi = 800)
    
    #Create a list of labels that will be appended at each iteration
    list_labels = []

    #iterate through the results of each subject
    for mouse, result in dict_EXTx.items():
        
        #add the mouse_id to the list of labels
        list_labels.append(mouse)
        
        #Create y axis :
        if type(result['Q']) == float or type(result['Q']) == int:
            y_axis_npa = i * np.ones(1)
        else:
            y_axis_npa = i * np.ones(len(result['Q']))
        
        #Convert timecourse from seconds to hours
        array_npa_hour = np.array(result['Q']) /60/60
        #Convert end time from minutes to hours
        end_session = result['D']/60
        
        #Plot npa timecourse:
        if result['Genotype'] =='wt':
            plt.plot(array_npa_hour, y_axis_npa, 'b|', markersize = 7, label = mouse)
        if result['Genotype'] =='he':
            plt.plot(array_npa_hour, y_axis_npa, 'g|', markersize = 7, label = mouse)
        if result['Genotype'] =='ho':
            plt.plot(array_npa_hour, y_axis_npa, 'r|', markersize = 7, label = mouse)
            
        #Plot the session end as black point for each animal:
        plt.plot(end_session, i, 'k.', markersize = 4)
        
        i += 1
    #Create the list of ticks matching the number of labels
    list_ticks = list(range(len(list_labels)))
    
    
    #Plot the ticks on the y axis
    plt.yticks(list_ticks, list_labels)
    
    plt.title(result['Protocol'])
    plt.xlabel('Time (h)')
    plt.ylabel('Mice ID')
    
    
def plot_timecourse_npa_pergenotype(dict_EXTx):
  
    plt.figure(num = 100, dpi = 800)
    plt.figure(num = 200, dpi = 800)
    plt.figure(num = 300, dpi = 800)
    
    y_wt = 0
    y_he = 0
    y_ho = 0
    
    list_labels_wt = []
    list_labels_he = []
    list_labels_ho = []
    
    
     #iterate through the results of each subject
    for mouse, result in dict_EXTx.items():
        #Convert timecourse from seconds to hours
        array_npa_hour = np.array(result['Q']) /60/60
        #Convert end time from minutes to hours
        end_session = result['D']/60
        
        #Plot npa timecourse:
        if result['Genotype'] =='wt':
            
            #add the mouse_id to the list of labels
            list_labels_wt.append(mouse)
            #Create y axis :
            if type(result['Q']) == float or type(result['Q']) == int:
                y_axis_npa = y_wt * np.ones(1)
            else:
                y_axis_npa = y_wt * np.ones(len(result['Q']))
                
            plt.figure(100)    
            plt.plot(array_npa_hour, y_axis_npa, 'b|', markersize = 7, label = mouse)
            plt.plot(end_session, y_wt, 'k.', markersize = 4)
            plt.title(result['Genotype'])
            
            y_wt += 1
            
        if result['Genotype'] =='he':
            #add the mouse_id to the list of labels
            list_labels_he.append(mouse)
            
            #Create y axis :
            if type(result['Q']) == float or type(result['Q']) == int:
                y_axis_npa = y_he * np.ones(1)
            else:
                y_axis_npa = y_he * np.ones(len(result['Q']))
                
            plt.figure(200)    
            plt.plot(array_npa_hour, y_axis_npa, 'g|', markersize = 7, label = mouse)
            plt.plot(end_session, y_he, 'k.', markersize = 4)
            plt.title(result['Genotype'])
            y_he  += 1
            
        if result['Genotype'] =='ho':
            #add the mouse_id to the list of labels
            list_labels_ho.append(mouse)
            
            #Create y axis :
            if type(result['Q']) == float or type(result['Q']) == int:
                y_axis_npa = y_ho * np.ones(1)
            else:
                y_axis_npa = y_ho * np.ones(len(result['Q']))
                
            plt.figure(300)    
            plt.plot(array_npa_hour, y_axis_npa, 'r|', markersize = 7, label = mouse)
            plt.plot(end_session, y_ho, 'k.', markersize = 4)
            plt.title(result['Genotype'])
            
            y_ho += 1
            
    list_ticks_wt = list(range(y_wt))
    list_ticks_he = list(range(y_he))
    list_ticks_ho = list(range(y_ho))
    
    plt.figure(100) 
    plt.yticks(list_ticks_wt, list_labels_wt)
    plt.xlabel('Time (h)')
    plt.ylabel('Mice ID')
    plt.tight_layout()
    
    plt.figure(200) 
    plt.yticks(list_ticks_he, list_labels_he)
    plt.xlabel('Time (h)')
    plt.ylabel('Mice ID')
    plt.tight_layout()
    
    plt.figure(300) 
    plt.yticks(list_ticks_ho, list_labels_ho)
    plt.xlabel('Time (h)')
    plt.ylabel('Mice ID')
    plt.tight_layout()
    
def plot_npa_timecourse_cumul(dictionary, fig_number):
    
    fig_social = fig_number + 1000
    fig_nonsocial = fig_number + 2000
    
    plt.figure(num = fig_social, dpi = 800)
    plt.figure(num = fig_nonsocial, dpi = 800)
    
    for subject, results in dictionary.items():
        
        genotype = results['Genotype']
        condition = results['Stimulus']
        
        x_timecourse_npa = np.array(results['Q'])/60/60
        
        y_timecourse_npa = np.ones(x_timecourse_npa.size).cumsum(axis = 0)
        
        if genotype == 'wt' and condition == 'social':
            plt.figure(fig_social)
            plt.plot(x_timecourse_npa, y_timecourse_npa, 'b-', alpha = 0.3)
        elif genotype == 'he' and condition == 'social':
            plt.figure(fig_social)
            plt.plot(x_timecourse_npa, y_timecourse_npa, 'g-', alpha = 0.3)
        elif genotype == 'ho' and condition == 'social':
            plt.figure(fig_social)
            plt.plot(x_timecourse_npa, y_timecourse_npa, 'r-', alpha = 0.3)
            
        if genotype == 'wt' and condition == 'nonsocial':
            plt.figure(fig_nonsocial)
            plt.plot(x_timecourse_npa, y_timecourse_npa, 'b-', alpha = 0.3)
        elif genotype == 'he' and condition == 'nonsocial':
            plt.figure(fig_nonsocial)
            plt.plot(x_timecourse_npa, y_timecourse_npa, 'g-', alpha = 0.3)
        elif genotype == 'ho' and condition == 'nonsocial':
            plt.figure(fig_nonsocial)
            plt.plot(x_timecourse_npa, y_timecourse_npa, 'r-', alpha = 0.3)
            
    plt.figure(fig_social)
    plt.title('Cumulative NPA - Extinction with social stimulus') 
    plt.ylabel('Number of Nose Pokes (cumulative)')
    plt.xlabel('Session time (h)')       
    plt.legend()
    
    plt.figure(fig_nonsocial)
    plt.title('Cumulative NPA - Extinction without social stimulus') 
    plt.ylabel('Number of Nose Pokes (cumulative)')
    plt.xlabel('Session time (h)')       
    plt.legend()
    
    
def plot_npa_timecourse_binned(dictionary, fig_number):
    
    fig_social = fig_number + 10000
    fig_nonsocial = fig_number + 20000
    
    plt.figure(num = fig_social, dpi = 800)
    plt.figure(num = fig_nonsocial, dpi = 800)
    
    list_timecourse_wt_social = []
    list_timecourse_he_social = []
    list_timecourse_ho_social = []
    
    list_timecourse_wt_nonsocial = []
    list_timecourse_he_nonsocial = []
    list_timecourse_ho_nonsocial = []
    
    x_axis_binned = np.arange(1,19)
  
    
    for subject, results in dictionary.items():
        
        genotype = results['Genotype']
        condition = results['Stimulus']
        protocol = results['Protocol']
        
        x_timecourse_npa = np.array(results['Q'])/60/60
        y_timecourse_npa = np.ones(x_timecourse_npa.size).cumsum(axis = 0)
        
        timecourse_npa_binned = binned_statistic(x_timecourse_npa, 
                                                 y_timecourse_npa, 
                                                 'count', 
                                                 bins = 18, 
                                                 range = (0,18)) #Highly important to homogenize since session duration is variable
        
        #Plot individual values:   
        if genotype == 'wt' and condition == 'social':
            plt.figure(fig_social)
            plt.plot(x_axis_binned, timecourse_npa_binned[0], 'b.', alpha = 0.5)
            list_timecourse_wt_social.append(timecourse_npa_binned[0])
        elif genotype == 'he' and condition == 'social':
            plt.figure(fig_social)
            plt.plot(x_axis_binned, timecourse_npa_binned[0], 'g.', alpha = 0.5)
            list_timecourse_he_social.append(timecourse_npa_binned[0])
        elif genotype == 'ho' and condition == 'social':
            plt.figure(fig_social)
            plt.plot(x_axis_binned, timecourse_npa_binned[0], 'r.', alpha = 0.5)
            list_timecourse_ho_social.append(timecourse_npa_binned[0])
            
            
        if genotype == 'wt' and condition == 'nonsocial':
            plt.figure(fig_nonsocial)
            plt.plot(x_axis_binned, timecourse_npa_binned[0], 'b.', alpha = 0.5)
            list_timecourse_wt_nonsocial.append(timecourse_npa_binned[0])
        elif genotype == 'he' and condition == 'nonsocial':
            plt.figure(fig_nonsocial)
            plt.plot(x_axis_binned, timecourse_npa_binned[0], 'g.', alpha = 0.5)
            list_timecourse_he_nonsocial.append(timecourse_npa_binned[0])
        elif genotype == 'ho' and condition == 'nonsocial':
            plt.figure(fig_nonsocial)
            plt.plot(x_axis_binned, timecourse_npa_binned[0], 'r.', alpha = 0.5)
            list_timecourse_ho_nonsocial.append(timecourse_npa_binned[0])
   
    #Convert list to np array and compute the mean 
    mean_wt_social = np.array(list_timecourse_wt_social).mean(axis = 0)
    sem_wt_social = scp_sem(np.array(list_timecourse_wt_social), axis = 0, ddof= 0)
    sem_sup_wtsocial = mean_wt_social + sem_wt_social
    sem_inf_wtsocial = mean_wt_social - sem_wt_social
    
    mean_he_social = np.array(list_timecourse_he_social).mean(axis = 0)
    sem_he_social = scp_sem(np.array(list_timecourse_he_social), axis = 0, ddof= 0)
    sem_sup_hesocial = mean_he_social + sem_he_social
    sem_inf_hesocial = mean_he_social - sem_he_social
    
    mean_ho_social = np.array(list_timecourse_ho_social).mean(axis = 0)
    sem_ho_social = np.array(scp_sem(np.array(list_timecourse_ho_social), axis = 0, ddof= 0))
    sem_sup_hosocial = mean_ho_social + sem_ho_social
    sem_inf_hosocial = mean_ho_social - sem_ho_social
    
    mean_wt_nonsocial = np.array(list_timecourse_wt_nonsocial).mean(axis = 0)
    sem_wt_nonsocial = np.array(scp_sem(np.array(list_timecourse_wt_nonsocial), axis = 0, ddof= 0))
    sem_sup_wtnonsocial = mean_wt_nonsocial + sem_wt_nonsocial
    sem_inf_wtnonsocial = mean_wt_nonsocial - sem_wt_nonsocial
    
    mean_he_nonsocial = np.array(list_timecourse_he_nonsocial).mean(axis = 0)
    sem_he_nonsocial = np.array(scp_sem(np.array(list_timecourse_he_nonsocial), axis = 0, ddof= 0))
    sem_sup_henonsocial = mean_he_nonsocial + sem_he_nonsocial
    sem_inf_henonsocial = mean_he_nonsocial - sem_he_nonsocial
    
    mean_ho_nonsocial = np.array(list_timecourse_ho_nonsocial).mean(axis = 0)
    sem_ho_nonsocial = np.array(scp_sem(np.array(list_timecourse_ho_nonsocial), axis = 0, ddof= 0))
    sem_sup_hononsocial = mean_ho_nonsocial + sem_ho_nonsocial
    sem_inf_hononsocial = mean_ho_nonsocial - sem_ho_nonsocial
    
    list_result = [mean_wt_social, sem_wt_social, sem_sup_wtsocial, sem_inf_wtsocial, 
                   mean_he_social, sem_he_social, sem_sup_hesocial, sem_inf_hesocial, 
                   mean_ho_social, sem_ho_social, sem_sup_hosocial, sem_inf_hosocial, 
                   mean_wt_nonsocial, sem_wt_nonsocial, sem_sup_wtnonsocial, sem_inf_wtnonsocial,  
                   mean_he_nonsocial, sem_he_nonsocial, sem_sup_henonsocial, sem_inf_henonsocial, 
                   mean_ho_nonsocial, sem_ho_nonsocial, sem_sup_hononsocial, sem_inf_hononsocial]

    
    plt.figure(fig_social)
    plt.plot(x_axis_binned, mean_wt_social, 
             'b-', label = 'wt')
    plt.plot(x_axis_binned, mean_he_social, 
             'g-', label = 'he')
    plt.plot(x_axis_binned, mean_ho_social, 
             'r-', label = 'ho')
    
    plt.plot(x_axis_binned, sem_inf_wtsocial, 
             x_axis_binned, sem_sup_wtsocial, 
             'b--',
             alpha = 0.1)
    plt.plot(x_axis_binned, sem_inf_hesocial, 
             x_axis_binned, sem_sup_hesocial, 
             'g--',
             alpha = 0.1)
    plt.plot(x_axis_binned, sem_inf_hosocial, 
             x_axis_binned, sem_sup_hosocial,  
             'r--',
             alpha = 0.1)
    
    plt.fill_between(x_axis_binned, 
                     list_result[2], list_result[3], 
                     color = 'b',
                     alpha = 0.1)
    plt.fill_between(x_axis_binned, 
                     list_result[6], list_result[7], 
                     color = 'g',
                     alpha = 0.1)
    plt.fill_between(x_axis_binned, 
                     list_result[10], list_result[11],  
                     color = 'r',
                     alpha = 0.1)
    
    plt.title(str('Extinction social '+protocol+' binned')) 
    plt.ylabel('Number of Nose Pokes (binned)')
    plt.xlabel('Session time (h)')       
    plt.legend()
    
    plt.figure(fig_nonsocial)
    plt.plot(x_axis_binned, mean_wt_nonsocial, 
             'b-', label = 'Dat-cre$^{het}$,Shank3$^{wt}$')
    plt.plot(x_axis_binned, mean_he_nonsocial, 
             'g-', label = 'Dat-cre$^{het}$,Shank3$^{het}$')
    plt.plot(x_axis_binned, mean_ho_nonsocial, 
             'r-', label = 'Dat-cre$^{het}$,Shank3$^{ho}$')
    
    plt.plot(x_axis_binned, sem_inf_wtnonsocial, 
             x_axis_binned, sem_sup_wtnonsocial,
             'b--',
             alpha = 0.1)
    plt.plot(x_axis_binned, sem_inf_henonsocial, 
             x_axis_binned, sem_sup_henonsocial,
             'g--',
             alpha = 0.1)
    plt.plot(x_axis_binned, sem_inf_hononsocial, 
             x_axis_binned, sem_sup_hononsocial,
             'r--',
             alpha = 0.1)
    
    plt.fill_between(x_axis_binned, 
                     list_result[14], list_result[15],
                     color = 'b',
                     alpha = 0.1)
    plt.fill_between(x_axis_binned, 
                     list_result[18], list_result[19],
                     color = 'g',
                     alpha = 0.1)
    plt.fill_between(x_axis_binned,
                     list_result[22], list_result[23],
                     color = 'r',
                     alpha = 0.1)
    
    plt.title(str('Extinction nonsocial '+protocol+' binned')) 
    plt.ylabel('Number of Nose Pokes (binned)')
    plt.xlabel('Session time (h)')       
    plt.legend()
    
    return list_result