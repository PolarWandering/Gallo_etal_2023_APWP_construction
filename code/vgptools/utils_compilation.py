import os
import numpy as np
from pmagpy import pmag, ipmag
import matplotlib.pyplot as plt
import pandas as pd


def split_datasheet (df_files, file_idx):
    """
    Reads in datasheets and splits them into pole and vgp collections to be filtered, collated and compiled.
    Input: standard vgp datasheet (as described in datasheet template)
    Output: separate dataframes comprised of the study-level poles and site-level vgps extracted from the datasheet
    """
    df = pd.read_excel(df_files['path'][file_idx]) #, skip_blank_lines=True

    df_poles = pd.read_excel(df_files['path'][file_idx], 
                             skiprows = df[df.iloc[:,0]=='Study level data'].index[0]+2,
                             nrows  = df[df.isnull().all(1)].index[1] -3)

    df_vgps = pd.read_excel(df_files['path'][file_idx], 
                            skiprows = df[df.iloc[:,0]=='Site level data'].index[0]+2)
    
    df_vgps = df_vgps.astype({'in_study_pole':int})
    
    return (df_poles, df_vgps)

def recalc_vgps(df_vgps):
    df_vgps['VGP_lon_recalc'] = df_vgps.apply(lambda row: pmag.dia_vgp(row.dec, row.inc, 1, row.slat, row.slon)[0], axis =1)
    df_vgps['VGP_lat_recalc'] = df_vgps.apply(lambda row: pmag.dia_vgp(row.dec, row.inc, 1, row.slat, row.slon)[1], axis =1)
    return df_vgps

def go_reverse(df_vgps):
    """
    Determines a new series where the vgps and directions are all of reversed polarity.
    Input: dataframe with directions and site coordinates. For the last tens of Ma it is reasonable to leave this as is, but this is
    not a safe assumption in deeper time (and thus this mean reverse pole may need to be set on a case-by-case basis).
    Output: the original dataframe and a new series with directions and vgps reported in reverse polarity
    """
    
    df_vgps['inc_reverse'] = np.where(df_vgps['inc'] > 0, -df_vgps['inc'], df_vgps['inc'])
    df_vgps['dec_reverse'] = np.where(df_vgps['inc'] > 0,(df_vgps['dec'] - 180.) % 360., df_vgps['dec'])
        
    df_vgps['vgp_lat_SH'] = np.where(df_vgps['VGP_lat_recalc'] > 0, -df_vgps['VGP_lat_recalc'], df_vgps['VGP_lat_recalc'])
    df_vgps['vgp_lon_SH'] = np.where(df_vgps['VGP_lat_recalc'] > 0,(df_vgps['VGP_lon_recalc'] - 180.) % 360., df_vgps['VGP_lon_recalc'])
    
    return df_vgps

def get_alpha95(df_vgps):
    
    df_vgps['alpha95'] = np.where(df_vgps['alpha95'].isna(), 140.0/np.sqrt(df_vgps['n'] * df_vgps['k']), df_vgps['alpha95'])
    return df_vgps

def get_k(df_vgps):
    
    df_vgps['k'] = np.where(df_vgps['k'].isna(), ((140./df_vgps['alpha95'])**2)/(df_vgps['n']), df_vgps['k'])
    return df_vgps

def get_ages(df_vgps, round_ceil = True):
    '''
    Given a DF with max_age and min_age filled field, this function calculates the mean age along with the age uncertainty
    The mean age can be rounded to the ceiling if <round_ceil = True>
    '''
    
    df_vgps['age_uncertainty'] = df_vgps['max_age'] - df_vgps['min_age']
    if round_ceil == True:
        df_vgps['mean_age'] = np.ceil((df_vgps['max_age'] + df_vgps['min_age']) / 2)
    else:
        df_vgps['mean_age'] = (df_vgps['max_age'] + df_vgps['min_age']) / 2
    return df_vgps

# def dfs_vgps_recomputed_poles(data_path_VGP, by_study = True):
    
#     '''
#     From the path where all the datasheets live, this function returns two DataFrames, one for
#     all the VGPs that were considered by the original author and other DataFrame with the poles.
#     The later has to flavours, one approach considers a region or datasheet, as a unit, from which
#     a paleopole may be computed, the other approach considers if there are two paleopoles from
#     the same area as two single paleopoles within the same area/study.
    
#     input: Pass a path where all the datasheets live. If by_study == True there will by one pole 
#     for each "Study",otherwise, for each study we evaluate if there is more than one pole for study 
#     and we return in accordance.
#     output: Selected entries that were considered by their authors to the calculation of PPs.
    
#     Note: This is achieved taking advantage on the values in column df_unfiltered[`in_study_pole`]. Zero
#     valued entries were discarder and integers labels different poles within the same study.
#     '''    

#     files_names = get_files_in_directory(data_path_VGP)
#     xlsx_file_names = [os.path.splitext(os.path.basename(open(file,'r').name))[0] for file in files_names if file.endswith('.xlsx')]
#     paths = [file for file in files_names if file.endswith('.xlsx')]
#     df_files = pd.DataFrame({'path': paths,  'name_xlsx': xlsx_file_names})

#     df_vgp_unfiltered, df_poles_original = merge_files(df_files)

#     if by_study == True: 
#         df_filtered, df_pole_compilation = original_selection(df_vgp_unfiltered, df_poles_original, by_study = True)
#     else:
#         df_filtered, df_pole_compilation = original_selection(df_vgp_unfiltered, df_poles_original, by_study = False)
    
#     df_pole_compilation = df_pole_compilation.astype({"N":int})
    
#     return df_filtered, df_pole_compilation    
    

# def merge_files(df_files):
    
#     for i in df_files.index:   # cycle over each file in database

#         # import data and assign to dataframes
#         df_poles_temp, df_vgps_temp = split_datasheet(df_files, i)

#         df_vgps_temp['rej_crit'] = [[int(float(j)) for j in str(i.rej_crit).split(';') if ~np.isnan(float(j))] for _,i in df_vgps_temp.iterrows()]
#         df_vgps_temp['Study'] = df_files.name_xlsx[i]
#         df_poles_temp['Study'] = df_files.name_xlsx[i]
#         # 
#         if not df_vgps_temp.empty:

#             df_vgps_temp = recalc_vgps(df_vgps_temp)       
#             df_vgps_temp = go_reverse(df_vgps_temp)        
#             df_vgps_temp = get_alpha95(df_vgps_temp)        
#             df_vgps_temp = get_k(df_vgps_temp)
#             df_vgps_temp = get_ages(df_vgps_temp, round_ceil = False)

#             df_vgps_temp['age_uncertainty'] = df_vgps_temp['max_age'] - df_vgps_temp['min_age']

#             if i == 0 : df_vgp_unfiltered = pd.DataFrame(data=None, columns=df_vgps_temp.columns); df_poles_original = pd.DataFrame(data=None, columns=df_poles_temp.columns)

#             # parse data
#             df_vgp_unfiltered = df_vgp_unfiltered.append(df_vgps_temp, ignore_index=True)
#             df_poles_original = df_poles_original.append(df_poles_temp, ignore_index=True)
    
#     return df_vgp_unfiltered, df_poles_original


# def original_selection(df_unfiltered, df_poles_original, by_study = True):
    
#     '''
#     input: Pass an unfiltered DF 
#     output: A selection of entries that were considered by their authors for the calculation of PPs.
    
#     Note: It is achieved by utilizing the values in column df_unfiltered[`in_study_pole`]. If the value 
#     in this column was zero (0), he authors discarded entries. If entries are considered for different 
#     studies, they are grouped according to a different value.
#     '''    
    
#     df_unfiltered['keep'] = np.nan    
#     df_unfiltered['keep'] = df_unfiltered.apply(lambda row: True if row.in_study_pole != 0 else row.keep, axis = 1)
#     df_filtered = df_unfiltered.loc[df_unfiltered['keep'] == True]
      
#     # iterate through each study in order to recompute and store the paleomagnetic poles
#     df_pole_compilation = pd.DataFrame(data = None, columns = df_poles_original.columns)    
    
    
#     if by_study == True:
        
#         # iterate through each study
#         for study, df_study in df_filtered.groupby('Study'):

#             ppole = ipmag.fisher_mean(dec = df_study['vgp_lon_SH'].tolist(), inc = df_study['vgp_lat_SH'].tolist()) # final paleopole
#             mean_site = ipmag.fisher_mean(dec = df_study['slon'].tolist(), inc = df_study['slat'].tolist())

#             if len(df_study) == 1: 
#                 mean_site['inc'] = df_study['slat'].values[0]; mean_site['dec'] = df_study['slon'].values[0] 
#                 ppole['inc'] = df_study['vgp_lat_SH'].values[0]; ppole['dec'] = df_study['vgp_lon_SH'].values[0]
#                 ppole['n'] = df_study['n'].values[0]; ppole['k'] = df_study['k'].values[0]; ppole['alpha95'] = df_study['alpha95'].values[0]


#             df_pole_compilation = df_pole_compilation.append({'Study': study, 
#                                                               'slat': mean_site['inc'], 'slon': mean_site['dec'],
#                                                               'Plat': ppole['inc'], 'Plon': ppole['dec'],
#                                                               'N': ppole['n'], 'K': ppole['k'], 'A95': ppole['alpha95'],
#                                                               'min_age': df_study.min_age.min(), 'max_age': df_study.max_age.max(), 
#                                                               'mean_age': (df_study.max_age.max() + df_study.min_age.min()) / 2 },
#                                                               ignore_index = True)
#     else:
        
#         for study, df_study in df_filtered.groupby('Study'):

#             # iterate through each study
#             for pole, df_pole in df_study.groupby('in_study_pole'):

#                 ppole = ipmag.fisher_mean(dec = df_pole['vgp_lon_SH'].tolist(), inc = df_pole['vgp_lat_SH'].tolist()) # final paleopole
#                 mean_site = ipmag.fisher_mean(dec = df_study['slon'].tolist(), inc = df_study['slat'].tolist())

#                 if len(df_pole) == 1: 
#                     mean_site['inc'] = df_pole['slat'].values[0]; mean_site['dec'] = df_pole['slon'].values[0] 
#                     ppole['inc'] = df_pole['vgp_lat_SH'].values[0]; ppole['dec'] = df_pole['vgp_lon_SH'].values[0]
#                     ppole['n'] = df_pole['n'].values[0]; ppole['k'] = df_pole['k'].values[0]; ppole['alpha95'] = df_pole['alpha95'].values[0]


#                 df_pole_compilation = df_pole_compilation.append({'Study': study, 'pole': pole, 
#                                                                   'slat': mean_site['inc'], 'slon': mean_site['dec'],
#                                                                   'Plat': ppole['inc'], 'Plon': ppole['dec'],
#                                                                   'N': ppole['n'], 'K': ppole['k'], 'A95': ppole['alpha95'],
#                                                                   'min_age': df_pole.min_age.min(), 'max_age': df_pole.max_age.max(), 
#                                                                   'mean_age': (df_pole.max_age.max() + df_pole.min_age.min()) / 2 }, 
#                                                                   ignore_index = True)

    
            
            
            
#     df_pole_compilation = df_pole_compilation[['Study','pole','N','K','A95','slat','slon','Plat','Plon','min_age','max_age','mean_age']]
#     return df_filtered, df_pole_compilation




# ####################### BEFORE MEETING 

# def merge_all_files(current_path):
#     '''
#     Given the path in which all the files live, this function generates a dataframe that gathers all the vgps and another
#     different DF for the reported poles.
#     '''
    
#     def get_df_files(current_path):
#         '''
#         from the current path, this function generates a dataframe with the name and complete path to the file within the folder
#         '''
#         data_path_VGP = current_path + '/vgp_database'
#         files_names = get_files_in_directory(data_path_VGP)
#         xlsx_file_names = [os.path.splitext(os.path.basename(open(file,'r').name))[0] for file in files_names if file.endswith('.xlsx')]
#         paths = [file for file in files_names if file.endswith('.xlsx')]
#         df_files = pd.DataFrame({'path': paths,  'name_xlsx': xlsx_file_names})
#         return df_files
    
#     def merge_unfiltered_datasheets(df_files):    
#         '''
#         from the current path, this function generates a dataframe that contains all the entries included in all the files in the current directory
#         input: dataframe with the name of the files and the complete to path to access.
#         output: Two different DataFrames. One containing all the vgp entries and other for the reported poles.
#         '''
#         for i in df_files.index:   # cycle over each file in database
#             print (f'processing file {i}')

#             # import data and assign to dataframes
#             df_poles_temp, df_vgps_temp = split_datasheet(df_files, i)
#             df_poles_temp['Study'] = df_files.name_xlsx[i]
            
#             df_vgps_temp['rej_crit'] = [[int(float(j)) for j in str(i.rej_crit).split(';') if ~np.isnan(float(j))] for _,i in df_vgps_temp.iterrows()]
#             df_vgps_temp['Study'] = df_files.name_xlsx[i]
#             df_poles_temp['Study'] = df_files.name_xlsx[i]
            
#             if not df_vgps_temp.empty:

#                 df_vgps_temp = recalc_vgps(df_vgps_temp)       
#                 df_vgps_temp = go_reverse(df_vgps_temp)        
#                 df_vgps_temp = get_alpha95(df_vgps_temp)        
#                 df_vgps_temp = get_k(df_vgps_temp)
#                 df_vgps_temp = get_ages(df_vgps_temp, round_ceil = False)

#                 df_vgps_temp['age_uncertainty'] = df_vgps_temp['max_age'] - df_vgps_temp['min_age']

#                 if i == 0 : df_vgp_unfiltered = pd.DataFrame(data=None, columns=df_vgps_temp.columns); df_poles_original = pd.DataFrame(data=None, columns=df_poles_temp.columns)

#                 # parse data
#                 df_vgp_unfiltered = df_vgp_unfiltered.append(df_vgps_temp, ignore_index=True)
#             df_poles_original = df_poles_original.append(df_poles_temp, ignore_index=True)    
               
#         return df_vgp_unfiltered, df_poles_original    
    
#     df_files = get_df_files(current_path)
#     df_vgp_unfiltered, df_poles_original = merge_unfiltered_datasheets(df_files)
    
#     return df_vgp_unfiltered, df_poles_original


# def generates_compilation(df_vgp_unfiltered, df_poles_original, incl_criteria):
#     '''
#     This function generates a dataframe that contains all the entries that pass the given including cirterias
#     '''


        
#     def separate_synch_unit_means(df_vgp_unfiltered):
#         '''
#         separate to a different DataFrame the means of the sites that were considered as synchronous by the original author.
#         '''
#         # Set a DataFrame for the means of the synchronous units
#         df_synch_unit_means = df_vgp_unfiltered[df_vgp_unfiltered['synch_unit'].str.contains("M", na=False)]
#         # Ignore the means of the synch_units and cast the synch_unit to int
#         # df_vgp_unfiltered = df_vgp_unfiltered[~df_vgp_unfiltered['synch_unit'].str.contains("M", na=False)]
#         # df_vgp_unfiltered['synch_unit'] = df_vgp_unfiltered['synch_unit'].astype(int)

#         return df_vgp_unfiltered, df_synch_unit_means
        
#     def parse_including_criteria(df_vgp_unfiltered):
#         '''
#         Send to columns the rejection criterias applied in the original reference. 
#         '''
#         # Boolean type criterias
#         df_vgp_unfiltered['author_selection'] = df_vgp_unfiltered.apply(lambda row: True if row['in_study_pole'] != 0 else False, axis=1)
#         df_vgp_unfiltered['undemagnetized'] = [True if (1 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['sample_count'] = [True if (2 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['alpha_95'] = [True if (3 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['overprints'] = [True if (4 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['remagnetizations'] = [True if (5 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['uncertain_struct'] = [True if (6 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['rotated'] = [True if (7 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['shallowed'] = [True if (8 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['anomalous_dir'] = [True if (9 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['uncertain_age'] = [True if (10 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['distinct_age'] = [True if (11 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['sub-time_units'] = [True if (12 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]
#         df_vgp_unfiltered['otherwise_rej'] = [True if (14 in i.rej_crit) else False for _,i in df_vgp_unfiltered.iterrows()]

#         return df_vgp_unfiltered
 
#     def selection_boolean_conditionals(df_unfiltered, incl_criteria):
#         '''
#         Pass an unfiltered DF of vgps and makes a subselection from a dictionary of criterias.
#         '''    

#         df_unfiltered['keep'] = np.nan

#         if incl_criteria['author_selection'] == True:        
#             df_unfiltered['keep'] = df_unfiltered.apply(lambda row: True if row.in_study_pole != 0 else row.keep, axis = 1)

#         for i in incl_criteria:                       
#             if incl_criteria[i] == True:  # boolean cases          
#                 df_unfiltered['keep'] = df_unfiltered.apply(lambda row: True if row[i] == True else row.keep, axis = 1)

#         if incl_criteria['sample_count']: df_unfiltered['keep'] = df_unfiltered.apply(lambda row: True if row['n'] >= incl_criteria['sample_count'] else False, axis = 1)
#         if incl_criteria['alpha_95']: df_unfiltered['keep'] = df_unfiltered.apply(lambda row: True if row['alpha95'] <= incl_criteria['alpha_95'] else False, axis = 1)
#         if incl_criteria['uncertain_age']: df_unfiltered['keep'] = df_unfiltered.apply(lambda row: True if row['age_uncertainty'] <= incl_criteria['uncertain_age'] else False, axis = 1)
#         if incl_criteria['rock_type']: df_unfiltered['keep'] = df_unfiltered.apply(lambda row: True if row['rock_typ_1'] == incl_criteria['rock_type'] else False, axis = 1)

#         # discard those discarded by other criterias
#         for i in incl_criteria:        
#             if incl_criteria[i] or i == 'rock_type': continue # boolean cases                    
#             df_unfiltered['keep'] = df_unfiltered.apply(lambda row: False if row[i] == True else row.keep, axis = 1)

#         df_filtered = df_unfiltered.loc[df_unfiltered['keep'] == True]
#         del df_unfiltered['keep']
#         return df_filtered
        
#     def selection_recursive_conditionals(df_selection, incl_criteria):
#         '''
#         '''
        
#         # Final DataFrame with the entries that follow all the specified criteria
#         df_vgp_compilation = pd.DataFrame(data = None, columns = df_vgp_unfiltered.columns)
#         df_pole_compilation = pd.DataFrame(data = None, columns = df_poles_original.columns)

#         # We have to iterate through the different studies to evaluate nested conditionals
#         for study, df_study in df_selection.groupby('Study'):

#             # list wtih ALL the synch_units
#             synch_units = [i for i in df_study['synch_unit'].unique() if isinstance(i, (float, int))] 
#             # means of the synch_units reported by authors (M entries)
#             synch_units_means = [int(i.replace('M', '')) for i in df_study['synch_unit'].unique() if isinstance(i, str)]
#             # entries with no reported mean by authors (but considered to belong to the same cooling unit)
#             synch_without_mean = list(set(synch_units) - set(synch_units_means)) 

#             # dim a DF to fill-in with those entries that pass the criterias
#             df_temp = pd.DataFrame(data = None, columns = df_selection.columns)

#             # Strip age distinct (and send to the final DF), then, continue working with the subselection (this is independet of whether or not the criteria was applied)
#             if (incl_criteria['distinct_age']) and (not df_study[df_study['distinct_age'] == True].empty): 
#                 df_vgp_compilation = df_vgp_compilation.append(df_study[df_study['distinct_age'] == True]) # send the vgps with distinct directly to the final compilation since there is no reason to discard them from it.
#                 df_study = df_study.loc[~df_study['distinct_age'] == True]

#             # Group by synchronous units to evaluate if they represent the same spot reading or not (the same spot reading should have high concentration parameter -- to be setted e.g. higher than 150)
#             for synch_unit, df_synch_unit in df_study.groupby('synch_unit'):

#                 if incl_criteria['sub-time_units']: 

#                     # if the synch_unit equals zero, it is considered as a sopt-reading, and thus, it goes directly to the compilation      
#                     if synch_unit == 0: df_temp = df_temp.append(df_synch_unit); continue

#                     # evaluate the concebtration of individual entries that were considered all together as spot reading and averaged accordingly (M)
#                     if isinstance(synch_unit, str):

#                         df_synch_unit = df_study[df_study['synch_unit'] == float(synch_unit.replace('M', ''))] #subselect entries 
#                         df_synch_unit['in_study_pole'] = df_synch_unit_means.loc[(df_synch_unit_means['synch_unit'] == synch_unit) & 
#                                                                                  (df_synch_unit_means['Study'] == study)]['in_study_pole'].to_list()[0]

#                         if df_synch_unit.shape[0] > 1: 

#                             synch_unit_mean = ipmag.fisher_mean(dec = df_synch_unit['vgp_lon_SH'].tolist(), inc = df_synch_unit['vgp_lat_SH'].tolist())

#                             if synch_unit_mean['k'] < 150: # append all the entries otherwise consider as a single mean
#                                 df_temp = df_temp.append(df_synch_unit)                                         
#                             else:
#                                 df_temp = df_temp.append(df_synch_unit_means.loc[(df_synch_unit_means['synch_unit'] == synch_unit) & 
#                                                                                  (df_synch_unit_means['Study'] == study)])        
#                     if synch_unit in synch_without_mean:                 
#                         if int(synch_unit) == 0: continue
#                         synch_unit_mean = ipmag.fisher_mean(dec = df_synch_unit['vgp_lon_SH'].tolist(), inc = df_synch_unit['vgp_lat_SH'].tolist())                        
#                         site = ipmag.fisher_mean(dec = df_synch_unit['slon'].tolist(), inc = df_synch_unit['slat'].tolist()) 

#                         df_temp = df_temp.append({'Study': study, 'slat': site['inc'], 'slon': site['dec'],                                          
#                                                   'vgp_lat_SH': synch_unit_mean['inc'], 'vgp_lon_SH': synch_unit_mean['dec'], 'k': synch_unit_mean['k'],'alpha95': synch_unit_mean['alpha95'],
#                                                   'mean_age': df_synch_unit['mean_age'].mean(), 'min_age': df_synch_unit['min_age'].mean(), 'max_age': df_synch_unit['max_age'].mean(),
#                                                   'synch_unit' : synch_unit},                   
#                                                    ignore_index = True)            

#                 # if we not evaluate the entries within the cooling units, we take the M entries as face value (and we average the entries with no reported mean -i.e. M entry)                               
#                 else:
#                     # spot readings (0) goes straight to the final compilation
#                     if synch_unit == 0: df_temp = df_temp.append(df_synch_unit); continue

#                     #print(df_temp[['Study','slat','synch_study']])


#                     # reported means (M's) goes sraight to the final compilation
#                     if isinstance(synch_unit, str):
#                         #print('-'+str(synch_unit))
#                         df_temp = df_temp.append(df_synch_unit)

#                     # if there are 
#                     if synch_unit in synch_without_mean: 
#                         if int(synch_unit) == 0: continue
#                         synch_unit_mean = ipmag.fisher_mean(dec = df_synch_unit['vgp_lon_SH'].tolist(), inc = df_synch_unit['vgp_lat_SH'].tolist())                        
#                         site = ipmag.fisher_mean(dec = df_synch_unit['slon'].tolist(), inc = df_synch_unit['slat'].tolist()) 

#                         df_temp = df_temp.append({'Study': study, 'slat': site['inc'], 'slon': site['dec'],                                          
#                                                   'vgp_lat_SH': synch_unit_mean['inc'], 'vgp_lon_SH': synch_unit_mean['dec'], 'k': synch_unit_mean['k'],'alpha95': synch_unit_mean['alpha95'],
#                                                   'mean_age': df_synch_unit['mean_age'].mean(), 'min_age': df_synch_unit['min_age'].mean(), 'max_age': df_synch_unit['max_age'].mean(),
#                                                   'synch_unit' : synch_unit},                   
#                                                    ignore_index = True)

#             if incl_criteria['anomalous_dir']: 
#                 before = len(df_temp)
#                 df_temp = discard_vgps_recursively(df_temp, incl_criteria['anomalous_dir'])
#                 print('discarded : ', before - len(df_temp))
#             ppole = ipmag.fisher_mean(dec = df_temp['vgp_lon_SH'].tolist(), inc = df_temp['vgp_lat_SH'].tolist()) # final paleopole
#             mean_site = ipmag.fisher_mean(dec = df_temp['slon'].tolist(), inc = df_temp['slat'].tolist())

#             df_vgp_compilation = df_vgp_compilation.append(df_temp.copy(deep=True), ignore_index = True)

#             if len(ppole) == 0: continue 
#             df_pole_compilation = df_pole_compilation.append({'Study': study, 
#                                                               'slat': mean_site['inc'], 'slon': mean_site['dec'],
#                                                               'Plat': ppole['inc'], 'Plon': ppole['dec'],
#                                                               'N': ppole['n'], 'K': ppole['k'], 'A95': ppole['alpha95'],
#                                                               'min_age': df_study.min_age.min(), 'max_age': df_study.max_age.max(), 
#                                                               'mean_age': (df_study.max_age.max() + df_study.min_age.min()) / 2 }, 
#                                                               ignore_index = True)           
    
#         return df_vgp_compilation, df_pole_compilation
               
 
#     df_vgp_unfiltered, df_synch_unit_means = separate_synch_unit_means(df_vgp_unfiltered)  
#     df_vgp_unfiltered = parse_including_criteria(df_vgp_unfiltered)
#     df_pre_selection = selection_boolean_conditionals(df_vgp_unfiltered, incl_criteria)
#     df_vgp_compilation, df_pole_compilation = selection_recursive_conditionals(df_pre_selection, incl_criteria)
    
#     return df_vgp_compilation, df_pole_compilation


# def get_files_in_directory(path): 
#     """
#     Retrieves file names from a directory \
#     \n\nInput: path = directory \
#     \n\nOutput: list of subdirectories
#     """

#     # The last conditional here is in order to ignore the /DS_store file in macs 
#     return [os.path.join(path, name) for name in os.listdir(path)
#             if (os.path.isfile(os.path.join(path, name)) and (not name.startswith('.')))  ]







# def discard_vgps_recursively(df_vgps, max_angle):
#     '''
#     takes a group of vgps and discards iteratively those entries that lie "angle" degrees away from the recomputed mean
#     '''
#     if len(df_vgps) < 4:
#         return df_vgps
    
#     ppole = ipmag.fisher_mean(dec = df_vgps['vgp_lon_SH'].tolist(), inc = df_vgps['vgp_lat_SH'].tolist())       
#     df_vgps['dist_2_mean'] = df_vgps.apply(lambda row: pmag.angle([row.vgp_lon_SH, row.vgp_lat_SH], [ppole['dec'], ppole['inc']]), axis=1)    
#     max_angle_temp = df_vgps['dist_2_mean'].max()
#     if max_angle_temp < max_angle: return df_vgps

#     while max_angle_temp > max_angle:
#         if len(df_vgps) < 4: return df_vgps
#         df_vgps = df_vgps[df_vgps['dist_2_mean'] < max_angle]
#         ppole = ipmag.fisher_mean(dec = df_vgps['vgp_lon_SH'].tolist(), inc = df_vgps['vgp_lat_SH'].tolist())       
#         df_vgps['dist_2_mean'] = df_vgps.apply(lambda row: pmag.angle([row.vgp_lon_SH, row.vgp_lat_SH], [ppole['dec'], ppole['inc']]), axis=1)    
#         max_angle_temp = df_vgps['dist_2_mean'].max()
#         if max_angle_temp < max_angle: return df_vgps
