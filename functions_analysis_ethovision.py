# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 16:54:37 2021 @author: redon

Contains the following functions:
    -xlsx_to_csv()
    -read_exp_info()
    -open_rawdata()
    -exp_info()
    -obtain_info()
    -csv_classify_protocol()
    -def session_dict()
    -row_difference()
    -basic_variables()
    -intertrial_interval()
    -basic_analysis_1file()
    -plot_results()
    -processing()
    -timecourse_npa()
    -timecourse_npa_binned()
    -plot_timecourse_npa()
    -plot_results_protocol()
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from multiprocessing import Pool
import re
import datetime
from scipy.stats import binned_statistic

os.chdir('C:\\Users\\redon\\Documents\\Data\\Operant\\OC_DatShank3_M_Autoopto_June21\\Raw_data_202109_DatShank_Operant')

def xlsx_to_csv(file):
    """Function automatically converting all .xlsx files of a folder into .csv file
    
    Arg:
        file = path to the file to be converted
        
    Return:
        None
    """
    
    #os.chdir('C:\\Users\\redon\\Documents\\Data\\202104_BG_IL_NP_FR_Door_OptoDOIZ_2Arenas\\Export Files')
    
    #Define the destination folder for the csv to be created
    #dest = 'C:\\Users\\redon\\Documents\\Data\\202104_BG_IL_NP_FR_Door_OptoDOIZ_2Arenas\\test'
    
    #Convert the excel file into a dataframe
    df = pd.read_excel(file)
       
    #Save the dataframe into a csv file
    df.to_csv(file+'.csv', index=False)    
       
def read_exp_info(file):
    """Function reading all information about the experiment. Use to orientate future computations.
    
    Arg:
        file = name of the file path to read
            type Str
            
    Return:
        df_exp_info = dataframe containing all data about experimental conditions
            type DataFrame
            
        protocol = name of the protocol used (SH =shaping, FR = fixed ratio, PR = progressive ratio, Ext = extinction)
            type Str
            
        genotype = name of the animal's Shank3 genotype (wt = wild type, he = heterozygote, ho = homozygote)
            type Str
        
        mouse_id = identification number of the animal
            type Str
        
        trial_duration = Duration of the session analysed
            type Datetime
    """
    
    #Parse the data into a DataFrame and define header and index column
    df_exp_info = pd.read_csv(file,
                              header = None,
                              names = ['Index','Info'], #the column with the exp info is called 'Info'
                              index_col = [0], #the first column is used as indexer
                              sep = ",",
                              decimal = ".",
                              usecols = [0,1], 
                              nrows=(42))
    
    #extract the relevant informations:
    duration_delta = pd.to_timedelta(df_exp_info.at['Trial duration', 'Info'])
    duration = duration_delta.total_seconds()
    protocol = df_exp_info.at['Protocol','Info']
    genotype = df_exp_info.at['Genotype_Shank3','Info']
    mouse_id = df_exp_info.at['ID', 'Info']
    trial_duration = df_exp_info.at['Trial duration', 'Info']
    
    return df_exp_info, duration, protocol, genotype, mouse_id, trial_duration

def open_rawdata(file):
    """Function allowing to import raw data from an excel file to a pandas DataFrame
    WARNING: Check the working directory!!
    
    Arg:
        file = name of the excel file 
            type str
        
    Return:
        df_exp_info = DataFrame containing the experiment informations (first rows from Ethovision's raw data)
            type DataFrame
            
        df_session  = DataFrame containing the raw data with Session Time in indices (step of 0.04 s)
            type DataFrame
    """
    #Parse experiment description into a DataFrame:
    df_exp_info = pd.read_csv(file,
                              header = None,
                              names = ['Index','Info'], #the column with the exp info is called 'Info'
                              index_col = [0], #the first column is used as indexer
                              sep = ",",
                              decimal = ".",
                              usecols = [0,1], 
                              nrows=(42))
    

    #Parse the data into a DataFrame and define header and index column
    df_session = pd.read_csv(file,
                             header = 42,
                             sep = ",",
                             decimal = ".",
                             skiprows=([0, 44]),
                             na_values=("-"))
    
    return df_exp_info, df_session

def exp_info(df_exp_info):
    """Function accessing the experiment informations (see return parameters) for a given animal's session
    
    Arg:
        df_info = DataFrame containing the experiment informations (first rows from Ethovision's raw data)
                type DataFrame
                
    Return:
        protocol = describe the protocol used during the session (SH, FR, PR, Ext)
                type str
            
        genotype = indicate the Shank3 genotype of the experimental animal
                type str
        
        mouse_id = indicate the exact identification number of the experimental animal  
                type str
    """
    
    #Define the DataFrame to analyse
    df = df_exp_info
    
    #Extract Protocol, Genotype and mouse id 
    duration_delta = pd.to_timedelta(df.at['Trial duration', 'Info'])
    duration = duration_delta.total_seconds()
    protocol = df.at['Protocol','Info']
    genotype = df.at['Genotype_Shank3','Info']
    mouse_id = df.at['ID', 'Info']

    return duration, protocol, genotype, mouse_id

def obtain_info(file):
    """Function automatizing the opening of experiment file and extracting experiment informations
    WARNING: uses open_rawdata() and exp_info()
    
    Arg:
        file = name of the file to analyse
            type Pathname
        
    Return:
        df_info = DataFrame containing the experiment informations (first rows from Ethovision's raw data)
                type DataFrame
            
        df_data  = DataFrame containing the raw data with Session Time in indices (step of 0.04 s)
                type DataFrame
        
        protocol = describe the protocol used during the session (SH, FR, PR, Ext)
                type str
            
        genotype = indicate the Shank3 genotype of the experimental animal
                type str
        
        mouse_id = indicate the exact identification number of the experimental animal  
                type str
    """
    #First, extract info and data into separate dataframes
    df_info, df_data = open_rawdata(file)
    
    #Extract experiment informations
    duration, protocol, genotype, mouse_id = exp_info(df_info)
    
    #Print experiment informations
    print('protocol:', protocol,'\n',
          'genotype:', genotype,'\n',
          'mouse_id:', mouse_id)
    
    return df_info, df_data, duration, protocol, genotype, mouse_id

def csv_classify_protocol(list_file):
    '''Function classifying a list of raw csv files based on their protocol
    
    Arg:
        list_file = List containing all raw csv files from Ethovision
            type List
            
    Return: 
        list_file_SH, list_file_FR1, list_file_FR3, list_file_PR = Lists of filne names 
            type Lists
    '''
    #Classify files on the protocol (SH, FR or PR)      
    list_csv = glob.glob('*.csv')
    
    #initiate the lists for each protocol
    list_file_SH=[]
    list_file_FR1=[]
    list_file_FR3=[]
    list_file_PR=[]
    list_file_reversal = []
    
    #iterate through the list of file names
    for file in list_csv:
        
        #extract their relevant informations using read_exp_info function
        df_exp_info, duration, protocol, genotype, mouse_id, trial_duration = read_exp_info(file)
        
        #Create patterns matching relevant protocol names
        pattern_SH = re.compile(r'SH.+')
        pattern_FR1 = re.compile(r'FR1\.[0-9]+')
        pattern_FR3 = re.compile(r'FR3\.[0-9]+')
        pattern_PR = re.compile(r'PR')
        pattern_reversal = re.compile(r'R_FR1\.[0-9]')
        
        #Classify files based on their protocol
        if re.match(pattern_SH, protocol) != None :
            list_file_SH.append(file)
        elif re.match(pattern_FR1, protocol) != None:
            list_file_FR1.append(file)
        elif re.match(pattern_FR3, protocol) != None:
            list_file_FR3.append(file)
        elif re.match(pattern_PR, protocol) != None:
            list_file_PR.append(file)
        elif re.match(pattern_reversal, protocol) != None:
            list_file_reversal.append(file)
       
            
    return list_file_SH, list_file_FR1, list_file_FR3, list_file_PR, list_file_reversal

def session_dict(list_file):
    '''Create a dictionary with keys = sessions (e.g. FR1.1) and values = list of filename 
    Arg:
        list_file = List of filename to be classified
            type list
        
    Return:
        dict_session = dictionary described above
            type Dict
    '''
    
    dict_session = {}
    
    for file in list_file:
        
        df_exp_info, duration, protocol, genotype, mouse_id, trial_duration = read_exp_info(file)
        
        if protocol in dict_session.keys():
            dict_session[protocol].append(file)
        else:
            dict_session[protocol] = [file]
        
    return dict_session


def row_difference(df, column):
    """Function computing the difference between value of index n and n-1 
    in a given column and returning the number of value = -1 (nb of events)
    
    Arg:
        df = Dataframe containing data of interest
                    type DataFrame
        
        column = name of the column to compute
                    type str
        
        parameter = name to give to the computation output
                    type str
                    
    Return:
        count = computation output
                    type int              
    """
   
    #Attribute the column to a variable col_name
    col_name = column
    
    
    #Initialize the list containing the result
    list_diff = []
    
    
    #Iterate through rows for the column of interest
    i=1
    
    for i in range(len(df.index)):
        
        #Compute the difference between row n and n-1 elementwise
        list_diff.append(col_name.iloc[i] - col_name.iloc[i-1])
     
        #increment i
        i +=1
        
    #Count the number of occurences of -1 (end of event) 
    count = int(list_diff.count(-1))
    
    return count

def basic_variables(df):
    """Function computing the basic operant conditioning variables (NPactive, NPinactive, Door opening, Light)
    WARNING: uses row_difference()
    
    Arg: 
        df = DataFrame containing the raw data from Ethovision excel sheet
            type DataFrame
        
    Return:
        npa = Number of active nosepokes 
            type int
        
        
        npi = Number of inactive nosepokes 
            type int
        
        door_open = Number of door opening 
            type int
        
        light = Number of optostimulation received
            type int        
    """
    df_data = df
    
    #Calculate number of NP active (based on row_difference function)
    npa = row_difference(df_data, df_data['NP Active'])
    
    #Calculate number of NP inactive (based on row difference function)
    npi  = row_difference(df_data, df_data['NP Inactive'])
    
    #Calculate number of door openin (based on row_difference function)
    door_open = row_difference(df_data, df_data['Door Open'])
                                                
    #Calculate number of light stimulation effectively received (based on row_difference function)
    light = row_difference(df_data, df_data['Light'])
    
    #Calculate the total time of stimulation for a given session
    total_time_stimulation = df_data['Light'].sum() * 0.04 #each row is 40ms
    
    #Calculate the mean time of stimulation per self admin (x s over 7s possible)
    if light == 0:
       mean_time_stimulation = np.nan 
    else:
        mean_time_stimulation = total_time_stimulation / light
    
    #Time spent in interaction zone when the door is either open (do) or closed (dc)
    do_iz = df_data['DO IZ'].sum() * 0.04
    dc_iz = df_data['DC IZ'].sum() * 0.04
    
    #Time spent in npa zone
    time_zone_npa = df_data['Zone NPA'].sum() * 0.04
    
    #Calculate the latency to the first NP in the session
    first_npa = np.nan
    for i in range(len(df_data.index)):
        if df_data['NP Active'].iloc[i] == 1:
            first_npa = df_data['Trial time'].iloc[i]
            break
        
    #Calculate the latency to the first self stimulation in the session
    first_light = np.nan
    for i in range(len(df_data.index)):
        if df_data['Light'].iloc[i] == 1:
            first_light = df_data['Trial time'].iloc[i]
            break
        
    return npa, npi, door_open, light, total_time_stimulation, mean_time_stimulation, do_iz, dc_iz, time_zone_npa, first_npa, first_light
  
def intertrial_interval(df):
    """Function computing all intertrial interval defined as the time ellapse between the door closing and the first nose poke active
    WARNING: uses row_difference()
        
    Arg:
        df : DataFrame containing the session's raw data (raw excel file from Ethovision)
            type DataFrame
            
    Return:
        array_events = 2d numpy array containing pairs (on axis 0) of time for (i) door closing and (ii) the next nose poke active
            type 2d numpy array
            
        array_ITI = 1d numpy array containing each ITI of the session
            type 1d numpy array                
    """
    
    #Initialise columns of interest
    col_npa = df['NP Active']
    col_door = df['Door Open']
    
    #Calculate number of door openin (based on row_difference function)
    #door_open = row_difference(df, col_door)
    
    ##Initiate a 2d array to store the pairs of values
    array_events = np.array([[],[]])
    
    #Initiate variables needed to iterate through rows
    i = 1
    j = 0
    
    #First step: fill the 2d array array_events with
    #in axis 0: time of door closing
    #in axis 1: time of the first NPA after door closing
    
    #iterate through rows of the dataframe
    for i in range(len(df.index)):
        
        #Compute the difference between row n and n-1 elementwise in the door column
        diff_rows_door = col_door.iloc[i] - col_door.iloc[i-1]
        
        #Detect the time at which the door closes
        if diff_rows_door == -1: 
             
             k = 1
             
             #Then iterate through the following rows to detect the first NP post trial         
             for k in range(len(df.index)-i):
                 
                #Compute the difference between row n and n-1 in the npa column
                diff_rows_column = col_npa.iloc[i+k] - col_npa.iloc[i + k-1]
                
                #Detect the time at which the first NPA is made:
                if diff_rows_column == 1:
                    
                    #Add this time in axis zero or array_events
                    time_door = df['Trial time'].iloc[i]
                    #Add this time to array_events axis 1 to complete the pair endtrial-npa
                    time_npa = df['Trial time'].iloc[i + k]
                    #Compile data into a 2d array
                    array_interval = np.array([[time_door],[time_npa]])
                    
                    #Update array_events with the new apir of values
                    array_events =  np.append(array_events, array_interval, axis = 1)
                    
                    j += 1
                    i += 1
                    
                    break
            
                else:
                    k += 1

        else:
            #increment i
             i +=1
          
    #Initiate the array receiving the ITI
    array_ITI = np.zeros(array_events[0].size)
    
    n = 0
    
    #Calculate each ITI corresponding to the difference between each element of array_ITI on axis 0
    for n in range(array_events[0].size):
        
        #Time in s of the n time door closing 
        time_door_closed = array_events[0, n]
        
        #Time in s of the first npa after the n time door closing
        time_npa = array_events[1, n]
        
        #Calculate the difference and add it to array_ITI
        iti = time_npa - time_door_closed
        
        #Update array_ITI with the new value
        array_ITI[n] = iti
        
        
        n += 1
        

    return array_events, array_ITI

def basic_analysis_1file(file):
    """Extract every discrete informations from an Ethovision csv
    
    Arg: 
        file = Path to the file to analyse
            type Str
    return
        df_info, df_data = dataframes containing (i)exp informations and (ii) variables status every 40ms
            type Dataframes
            
        duration, protocol, genotype, mouse_id = Main experiment informations for the analysed session
            types Time, Str, Str, Str
            
        npa, npi, door_open, light, total_time_stimulation, mean_time_stimulation, do_iz, dc_iz, time_zone_npa, first_npa, first_light, array_ITI
    """
    
    df_info, df_data, duration, protocol, genotype, mouse_id = obtain_info(file)
    
    npa, npi, door_open, light, total_time_stimulation, mean_time_stimulation, do_iz, dc_iz, time_zone_npa, first_npa, first_light = basic_variables(df_data)
    
    array_events, array_ITI = intertrial_interval(df_data)
    
    return df_info, df_data, duration, protocol, genotype, mouse_id, npa, npi, door_open, light, total_time_stimulation, mean_time_stimulation, do_iz, dc_iz, time_zone_npa, first_npa, first_light, array_ITI

def processing(file_name):
    
    data = pd.DataFrame(columns =['Protocol', 'Duration', 'Genotype','mouse_id', 'NPA', 'NPI', 'Door Open', 'Light', "Mean ITI", 'Total time light', 'Mean time light', 'Time Zone NPA', "First NPA", 'First Light'])
    
    df_info, df_data, duration, protocol, genotype, mouse_id, npa, npi, door_open, light, total_time_stimulation, mean_time_stimulation, do_iz, dc_iz, time_zone_npa, first_npa, first_light, array_ITI = basic_analysis_1file(file_name)
    
    
    data.loc[0, 'Protocol'] = protocol
    data.loc[0, 'Duration'] = duration
    data.loc[0, 'Genotype'] = genotype
    data.loc[0, 'mouse_id'] = mouse_id
    data.loc[0, 'NPA'] = npa
    data.loc[0, 'NPI'] = npi
    data.loc[0, 'Door Open'] = door_open
    data.loc[0, 'Light'] = light
    
    if array_ITI.size == 0:
        data.loc[0, 'Mean ITI'] = np.nan
    elif array_ITI.size == 1:
        data.loc[0, 'Mean ITI'] = array_ITI[0]
    elif array_ITI.size > 1:
        data.loc[0, 'Mean ITI'] = array_ITI.mean()
        
    data.loc[0, 'Total time light'] = total_time_stimulation
    data.loc[0, 'Mean time light'] = mean_time_stimulation 
    data.loc[0, 'DO IZ'] = do_iz
    data.loc[0, 'DC IZ'] = dc_iz
    data.loc[0, 'Time Zone NPA'] = time_zone_npa
    data.loc[0, 'First NPA'] = first_npa
    data.loc[0, 'First Light'] = first_light
    
    return data

def timecourse_npa(list_csv_genotype):
    '''Generate npa timecourse dataframes (1/genotype), using row difference, 
    a value of 1 correspond to the initiation of a NP and -1 to its termination
    
    Arg:
        list_csv_genotype = List containing the name of the csv files for each session and individuals
            type List of str
    Return:
        df_npa_wt, df_npa_he, df_npa_ho  = dataframes containing the npa timecourses per genotype
            type Dataframes
    '''
    #Read one of the csv to extract the 'Trial time' column (common to all genotypes):
    df_info_initiation, df_data_initiation = open_rawdata(list_csv_genotype[0])
    col_time = df_data_initiation.loc[:,'Trial time']
    
    #Initiate dataframes for each genotype:
    df_npa_wt = pd.DataFrame(col_time)
    df_npa_he = pd.DataFrame(col_time)
    df_npa_ho = pd.DataFrame(col_time)
    
    for session in list_csv_genotype:
        
        #Call the exp dataframe:
        df_info, df_data = open_rawdata(session)
        id_protocol = str(df_info.at['ID', 'Info'] + '_' + df_info.at['Protocol', 'Info'])
        diff_npa = df_data['NP Active'].diff()
        
        if df_info.at['Genotype_Shank3', 'Info'] == 'wt':        
            df_npa_wt[id_protocol] = diff_npa
            
        elif df_info.at['Genotype_Shank3', 'Info'] == 'he':        
            df_npa_he[id_protocol] = diff_npa
            
        elif df_info.at['Genotype_Shank3', 'Info'] == 'ho':        
            df_npa_ho[id_protocol] = diff_npa
            
    df_npa_wt = df_npa_wt.set_index('Trial time')
    df_npa_he = df_npa_he.set_index('Trial time')
    df_npa_ho = df_npa_ho.set_index('Trial time')
        
    return df_npa_wt, df_npa_he, df_npa_ho    

def timecourse_npa_binned(df_npa, bins_df):
    '''Compute the sum binned statistics across all columns of a df
    
    Arg:
        df_npa = Dataframe of npa timecourse based on diff() method. Index is the session time
            type DataFrame
            
        bins_df = list of bins label. e.g. list(range(0,22,2)) for 10 bins between 0 and 20
            type List 
            
    Return:
        df_stats = Dataframe which index correspond to bins and columns to each mouse binned stats
            type DataFrame
    '''
    
    #Create a mask where 1 == True:
    npa_time = (df_npa == 1)
    df_npa_masked = df_npa[npa_time]
    
    #Initiate the df hosting the stats
    df_stats = pd.DataFrame(index = bins_df)
    
    i=0
    #Iterate through columns
    for col in df_npa_masked.columns:
        
        #Compute the binned statistics across indexes the the given column
        stats, edge, number = binned_statistic(df_npa_masked.index,
                                               df_npa_masked.iloc[:,i].fillna(0),
                                               'sum',
                                               bins = 10)
        #Add the calculated stats to the dataframe
        df_stats[col] = stats
        i+=1
    
    return df_stats
        
def plot_timecourse_npa(df_stats_wt, df_stats_he, df_stats_ho, protocol, bins, y_top = 0):
    '''Plot binned npa timecourse for each genotype
    
    Arg:
        df_stats_wt, df_stats_he, df_stats_ho = timecourse dfs (based on timecourse_npa_binned())
            type Dataframes
        protocol = Name of the protocol to apply an informative title to the plot
            type Str
            
    '''
    #Set the x, y and error for wt
    x_wt = df_stats_wt.index
    y_wt = df_stats_wt.mean(axis = 1)
    error_wt = df_stats_wt.sem(axis = 1)
    
    #Set the x, y and error for he
    x_he = df_stats_he.index
    y_he = df_stats_he.mean(axis = 1)
    error_he = df_stats_he.sem(axis = 1)
    
    #Set the x, y and error for ho
    x_ho = df_stats_ho.index
    y_ho = df_stats_ho.mean(axis = 1)
    error_ho = df_stats_ho.sem(axis = 1)
    
    #Plot each genotype mean and positive/negative error (sem)
    plt.figure(dpi = 800)
    plt.errorbar(x_wt, y_wt, error_wt, 
                 fmt = 'b-o', ms = 2, lw = 1,
                 ecolor = 'b', elinewidth = 0.3,
                 label = 'Dat-cre$^{het}$,Shank3$^{wt}$')
    plt.errorbar(x_he, y_he, error_he, 
                 fmt = 'g-o', ms = 2, lw = 1,
                 ecolor = 'g', elinewidth = 0.3,
                 label = 'Dat-cre$^{het}$,Shank3$^{het}$')
    plt.errorbar(x_ho, y_ho, error_ho, 
                 fmt = 'r-o', ms = 2, lw = 1,
                 ecolor = 'r', elinewidth = 0.3,
                 label = 'Dat-cre$^{het}$,Shank3$^{ho}$')
    
    #Format the x and y axis
    font_labels = {'family': 'arial',
                   'color':  'black',
                   'weight': 'normal',
                   'size': 12}
    plt.xticks(bins)
    plt.xlabel('Time within session (min)', fontdict = font_labels)
    
    if y_top != 0:
        plt.ylim(top = y_top)
    plt.ylim(0)
    plt.ylabel('Number of Nose Pokes', fontdict = font_labels)
    
    #Format title
    font_title = {'family': 'arial',
                  'color':  'black',
                  'weight': 'bold',
                  'size': 16}
    plt.title('Within session NPA timecourse' + ' ' + protocol,
              fontdict = font_title)
    
    #Add the legend
    plt.legend(fontsize = 8, frameon =False) 
    
    #Remove top and right spines
    plt.gca().spines[['top', 'right']].set_visible(False)
    
def plot_results_protocol(df, protocol, first_session = 0, last_session = 0, comment = ''):
    '''Function plotting individual and mean values (+-SEM) for each genotype in each session
   WARNING: WORK ON DISCRETE VALUES, ONLY SUITABLE FOR FR1 AND FR3
    Arg:
        df = Dataframe containing all the data for each animal in each session
           type DataFrame
       
        protocol = Name of the protocol to analyse (SH, FR1, FR3, PR)
            type Str
            
        first_session, last_session = Session number to analyse (e.g. for FR1 to FR5, first = 1, last = 5)
            type Int
            
        comment = Allows to add precision in the figure title
            type Str
      
    Return:
        None
    '''
    #Create a mask to select specific protocols
    list_FR1 = ['FR1.1', 'FR1.2', 'FR1.3', 'FR1.4', 'FR1.5', 'FR1.6', 'FR1.7', 'FR1.8', 'FR1.9', 'FR1.10']
    list_FR3 = ['FR3.1', 'FR3.2', 'FR3.3', 'FR3.4', 'FR3.5', 'FR3.6', 'FR3.7', 'FR3.8', 'FR3.9', 'FR3.10']
    
    if protocol == 'FR1':
        list_protocol = list_FR1[first_session - 1 : last_session]
        mask_protocol = (df['Protocol'].isin(list_protocol))
        df_toplot = df[mask_protocol]
        
    elif protocol == 'FR3':
        list_protocol = list_FR3[first_session - 1 : last_session]
        mask_protocol = (df['Protocol'].isin(list_protocol))
        df_toplot = df[mask_protocol]
    else:
        df_toplot = df
    #Create a mask to select specific genotype  
    wt =(df_toplot['Genotype'] == 'wt')
    he =(df_toplot['Genotype'] == 'he')
    ho =(df_toplot['Genotype'] == 'ho')
    
    #calculate the mean value for each genotype during each session
    wt_mean_values = df_toplot[wt].groupby(['Protocol']).mean().dropna()
    he_mean_values = df_toplot[he].groupby(['Protocol']).mean().dropna()
    ho_mean_values = df_toplot[ho].groupby(['Protocol']).mean().dropna()
    
    #calculate the SEM for each genotype during each session
    wt_sem_values = df_toplot[wt].groupby(['Protocol']).sem().dropna()
    he_sem_values = df_toplot[he].groupby(['Protocol']).sem().dropna()
    ho_sem_values = df_toplot[ho].groupby(['Protocol']).sem().dropna()
    
    dict_parameters = {'NPA' : 'Number of Nose Pokes',
                       'NPI': 'Number of Nose Pokes',
                       'Door Open' : 'Number of door openings',
                       'Light' : 'Number of optostimulations',
                       'Mean ITI' : 'Intertrial interval (s)',
                       'Total time light' : 'Time (s)',
                       'Mean time light': 'Time (s)',
                       'DO IZ' : 'Time (s)',
                       'DC IZ' : 'Time (s)',
                       'Time Zone NPA' : 'Time (s)',
                       'First NPA': 'Latency (s)',
                       'First Light': 'Latency (s)'}
    
    #For each parameter, plot the data
    for parameter, description  in dict_parameters.items():
        
        plt.figure(dpi = 800)
   
        #Plotting individual values npa
        plt.plot(df_toplot[wt].loc[:, 'Protocol'], df_toplot[wt].loc[:, parameter], "b.", 
                 markersize = 3,
                 label = 'Dat-cre$^{het}$,Shank3$^{wt}$')
        plt.plot(df_toplot[he].loc[:, 'Protocol'], df_toplot[he].loc[:, parameter], "g.",
                 markersize = 3,
                 label = 'Dat-cre$^{het}$,Shank3$^{het}$')
        plt.plot(df_toplot[ho].loc[:, 'Protocol'], df_toplot[ho].loc[:, parameter], "r.",
                 markersize = 3,
                 label = 'Dat-cre$^{het}$,Shank3$^{ho}$')
        
        #Plotting mean values npa
        plt.plot(wt_mean_values.index, wt_mean_values.loc[:, parameter], 'b-',
                 he_mean_values.index, he_mean_values.loc[:, parameter], 'g-',
                 ho_mean_values.index, ho_mean_values.loc[:, parameter], 'r-',
                 alpha = 0.7)
    
        #Plotting SEM and filling area between SEM
        
        plt.plot(wt_mean_values.index, wt_mean_values.loc[:,parameter] + wt_sem_values.loc[:,parameter], 'b--',
                 he_mean_values.index, he_mean_values.loc[:,parameter] + he_sem_values.loc[:,parameter], 'g--',
                 ho_mean_values.index, ho_mean_values.loc[:,parameter] + ho_sem_values.loc[:,parameter], 'r--',
                 linewidth = 0.1,
                 alpha = 0.7)
        plt.plot(wt_mean_values.index, wt_mean_values.loc[:,parameter] - wt_sem_values.loc[:,parameter], 'b--',
                 he_mean_values.index, he_mean_values.loc[:,parameter] - he_sem_values.loc[:,parameter], 'g--',
                 ho_mean_values.index, ho_mean_values.loc[:,parameter] - ho_sem_values.loc[:,parameter], 'r--',
                 linewidth = 0.1,
                 alpha = 0.7)
        
    
        plt.fill_between(wt_mean_values.index, 
                         wt_mean_values.loc[:,parameter] - wt_sem_values.loc[:,parameter], 
                         wt_mean_values.loc[:,parameter] + wt_sem_values.loc[:,parameter], 
                         color = 'b', 
                         alpha = 0.1)
        plt.fill_between(he_mean_values.index, 
                         he_mean_values.loc[:,parameter] - he_sem_values.loc[:,parameter], 
                         he_mean_values.loc[:,parameter] + he_sem_values.loc[:,parameter], 
                         color = 'g', 
                         alpha = 0.1)
        plt.fill_between(ho_mean_values.index, 
                         ho_mean_values.loc[:,parameter] - ho_sem_values.loc[:,parameter], 
                         ho_mean_values.loc[:,parameter] + ho_sem_values.loc[:,parameter], 
                         color = 'r', 
                         alpha = 0.1)
        font = {'family':'Arial',
                'color' : 'black',
                'weight' : 'bold',
                'size' : 12}
        
        font_title = {'family': 'arial',
                  'color':  'black',
                  'weight': 'bold',
                  'size': 16}
        plt.title(str(parameter) + ' ' + protocol + ' ' + comment,
                  fontdict = font_title)
        plt.xlabel('Session', fontdict = font)
        plt.xticks(rotation = 20)
        plt.ylabel(str(description), 
                   fontdict = font)
        plt.legend(frameon = False, fontsize = 8)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
    
        plt.ylim(bottom = 0)
        
        #The maximum of light stimulation being 7sec per sequence
        if parameter == 'Mean time light':
            plt.ylim(top = 7)
        #plt.savefig(parameter)
        
def plot_timecourse_session(list_file, bins, protocol, y_top = 0):
    
    df_npa_wt, df_npa_he, df_npa_ho = timecourse_npa(list_file)
    
    df_stats_wt = timecourse_npa_binned(df_npa_wt, bins)
    df_stats_he = timecourse_npa_binned(df_npa_he, bins)
    df_stats_ho = timecourse_npa_binned(df_npa_ho, bins)
    
    plot_timecourse_npa(df_stats_wt, df_stats_he, df_stats_ho, protocol, bins, y_top)