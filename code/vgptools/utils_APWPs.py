import numpy as np
import pandas as pd
from pmagpy import pmag, ipmag

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.geodesic import Geodesic
from shapely.geometry import Polygon

from vgptools.utils import spherical2cartesian, shape, eigen_decomposition, cartesian2spherical, GCD_cartesian, get_angle, PD


    
def running_mean_APWP_shape(data, plon_label, plat_label, age_label, window_length, time_step, max_age, min_age):
    """
    function to generate running mean APWP..
    """
    
    mean_pole_ages = np.arange(min_age, max_age + time_step, time_step)
    
    running_means = pd.DataFrame(columns=['age','N','n_studies','k','A95','csd','plon','plat', 'foliation','lineation','collinearity','coplanarity','elong_dir'])
    
    for age in mean_pole_ages:
        window_min = age - (window_length / 2.)
        window_max = age + (window_length / 2.)
        poles = data.loc[(data[age_label] >= window_min) & (data[age_label] <= window_max)]

        if poles.empty: continue
        
        number_studies = len(poles['Study'].unique())
        mean = ipmag.fisher_mean(dec=poles[plon_label].tolist(), inc=poles[plat_label].tolist())
        
        ArrayXYZ = np.array([spherical2cartesian([np.radians(i[plat_label]), np.radians(i[plon_label])]) for _,i in poles.iterrows()])        
        if len(ArrayXYZ) > 3:
            shapes = shape(ArrayXYZ)
            PrinComp=PD(ArrayXYZ)
            eVal, eVec = eigen_decomposition(ArrayXYZ)
            elong_dir = np.degrees(cartesian2spherical(eVec[:,1]))[1] # from T&K2004 (declination od the intermediate Evec)
            # mean['inc']=np.degrees(cartesian2spherical(PrinComp))[0]
            # mean['dec']=np.degrees(cartesian2spherical(PrinComp))[1]
        else:
            shapes = [np.nan,np.nan,np.nan,np.nan]
        
        if len(poles)>2: #ensures that dict isn't empty
            running_means.loc[age] = [age, mean['n'], number_studies, mean['k'],mean['alpha95'], mean['csd'], mean['dec'], mean['inc'], 
                                      shapes[0], shapes[1], shapes[2], shapes[3], elong_dir]
    # Set longitudes in [-180, 180]
    running_means['plon'] = running_means.apply(lambda row: row.plon - 360 if row.plon > 180 else row.plon, axis =1)   
    
    # The following block calculates rate of polar wander (degrees per million years) 
    running_means['PPcartesian'] = running_means.apply(lambda row: spherical2cartesian([np.radians(row['plat']),np.radians(row['plon'])]), axis = 1)
    running_means['PP_prev'] = running_means['PPcartesian'].shift(periods = 1)
    running_means['PP_next'] =  running_means['PPcartesian'].shift(periods = -1)
    running_means['GCD'] = running_means.apply(lambda row: np.degrees(GCD_cartesian(row['PP_prev'], row['PPcartesian'])), axis = 1)
    running_means['APW_rate'] = running_means['GCD']/running_means['age'].diff()
    # Calculate a 'kink' angle for each position of the path
    running_means['angle'] = running_means.apply(lambda row: get_angle(row['PP_prev'], row['PPcartesian'], row['PP_next']), axis = 1)

    running_means = running_means.drop(['PPcartesian', 'PP_prev', 'PP_next'], axis=1)      
    running_means.reset_index(drop=1, inplace=True)
    
    #set the present day field for the present
    running_means['plat'] = np.where(running_means['age']==0, -90, running_means['plat'])
    running_means['plon'] = np.where(running_means['age']==0, 0, running_means['plon'])
    
    return running_means

def get_pseudo_vgps(df):  
    '''
    takes a DF with paleomagnetic poles and respective statistics, it draws N randomly generated VGPs
    following the pole location and kappa concentration parameter. In the present formulation we follow
    a very conservative apporach for the assignaiton of ages to each VGP, it is taken at random between
    the lower and upper bounds of the distribution of reported VGPs.
    Note: column labels are presently hard-coded into this, if relevant.
    '''
    
    data = {'Study': [], 'Plat': [], 'Plon': [], 'mean_age': []}

    for index, row in df.iterrows():
        directions_temp = ipmag.fishrot(k = row.K, n = row.N, dec = row.Plon, inc = row.Plat, di_block = False)
        
        vgp_lon_bst = directions_temp[0]
        vgp_lat_bst = directions_temp[1]
    
        if row.uncer_dist == 'uniform':
            ages = [np.random.randint(np.floor(row.min_age),np.ceil(row.max_age)) for _ in range(row.N)]    
        elif row.uncer_dist == 'normal':
            ages = [np.random.normal(row.mean_age,(row.max_age - row.min_age) / 2) for _ in range(row.N)]
            
        studies = [row.Study for _ in range(row.N)]
 
        data['Study'] += studies
        data['Plat'] += vgp_lat_bst
        data['Plon'] += vgp_lon_bst
        data['mean_age'] += ages
    
    pseudo_vgps = pd.DataFrame(data)

    return pseudo_vgps



def get_vgps_sampling_from_direction(df, study_label= 'Study',
                                     slat_label='slat', slon_label='slon', 
                                     dec_label='dec_reverse', inc_label='inc_reverse', k_label='k',
                                     mean_age_lab='mean_age', min_age_lab='min_age', max_age_lab='max_age'):
    
    '''
    Input:
    
    DF with the following site-level information: 
    - Study, site coordinates, mean direction, concentration parameter, mean age and error distribution. 
    
    Steps:
    1. generate a pseudo-sample from the original dataset using nonparametric random sampling 
    with replacement (bootstrap sample).
    2. For each row (site-level entry) in the bootstrap sample draws a random direction following 
    the kappa concentration parameter and mean direction. 
    3. Assing an age to the entry samopling from the corresponding error distribution
    3. For each row, calculates the correspongin VGP. 
    
    Note: input directions must be all in the same mode (sensu pmagPy) so that we get coherent vgps. 
    Given that we want southern hemisphere vgps, we will work with reversed (<0) inclinations (step already
    done in 0_Preprocessing.ipynb)    
    
    Output:
    
    - A DataFrame with the same size than the original dataset, with randomized parameters.  
    '''    
    Study, age_bst, decs, incs, slat, slon, indexes = [], [], [], [], [], [], [] # parameters of the pseudo-sample to be filled
    
    k_mean = df[k_label].mean() # if site-level data has no reported kappa, we take the mean of the population instead.
    
    df_pseudo = df.sample(frac = 1, replace = True) # generates a bootstrapped sample of the dataframe by randomly sampling with replacement
    
    for index, row in df_pseudo.iterrows():        
        
        # we first generate one random direction from the original entry.
        kappa = k_mean if np.isnan(row[k_label]) else row[k_label] # if we don't have kappa, we take the mean of the reported ones       
        
        directions_temp = ipmag.fishrot(k = kappa, n = 1, dec = row[dec_label], inc = row[inc_label], di_block = False)
        
        decs.append(directions_temp[0][0])
        incs.append(directions_temp[1][0])
        slat.append(row[slat_label])
        slon.append(row[slon_label])
        indexes.append(index)
        Study.append(row[study_label])
        
        # Assessing the uncertianty distribution (uniform or normal)
        if row.uncer_dist == 'uniform':
            age_bst.append(np.random.randint(np.floor(row[min_age_lab]),np.ceil(row[max_age_lab])))
        else:            
            age_bst.append(np.random.normal(row[mean_age_lab], (row[max_age_lab] - row[mean_age_lab]) / 2)) 
            # the files were completed in such a way to have min and max ages, so we take sigma as a half of that range
    
    dictionary = {
                  'Study': Study,
                  'age': age_bst,
                  'dec': decs,    
                  'inc': incs,
                  'slat': slat,
                  'slon': slon 
                  }    
    new_df = pd.DataFrame(dictionary)        
       
    #new_df['plon'] = new_df.apply(lambda row: pmag.dia_vgp(row.dec, row.inc, 1, row.slat, row.slon)[0], axis =1)
    new_df['plon'] = pmag.dia_vgp(new_df['dec'], new_df['inc'], 1, new_df['slat'], new_df['slon'])[0]
    
    #new_df['plat'] = new_df.apply(lambda row: pmag.dia_vgp(row.dec, row.inc, 1, row.slat, row.slon)[1], axis =1)
    new_df['plat'] = pmag.dia_vgp(new_df['dec'], new_df['inc'], 1, new_df['slat'], new_df['slon'])[1]
    
    # set longitude in [-180,180]
    #new_df['plon'] = new_df.apply(lambda row: row.plon - 360 if row.plon > 180 else row.plon, axis =1)
    new_df['plon'] = new_df['plon'].where(new_df['plon'] <= 180, new_df['plon'] - 360)

    new_df.index = indexes

    return new_df



def MC_error_prop_ensemble_results(df_vgps_original, n_sims = 100, 
                                   study_label= 'Study', slat_label='slat', slon_label='slon', 
                                   dec_label='dec_reverse', inc_label='inc_reverse', k_label='k',
                                   mean_age_lab='mean_age', min_age_lab='min_age', max_age_lab='max_age',
                                   plon_label = 'plon', plat_label='plat', age_label = 'age',
                                   window_length=20, time_step=1, max_age=65, min_age=0):
    
    '''
    This function does everyting 
    '''
    
    running_means_global = pd.DataFrame(columns=['run','N','k','A95','csd','foliation','lineation','collinearity','coplanarity'])
    
    for i in range(n_sims):
    
        # Generate a pseudo-sample of the original dataset in which every entry is a pseudo-sample taken for the error PDF
        pseudo_df = get_vgps_sampling_from_direction(df_vgps_original, study_label= study_label,
                                             slat_label=slat_label, slon_label=slon_label, 
                                             dec_label=dec_label, inc_label=inc_label, k_label=k_label,
                                             mean_age_lab=mean_age_lab, min_age_lab=min_age_lab, max_age_lab=max_age_lab)

        # Construct a Moving Average on the former data-set
        RM = running_mean_APWP_shape(pseudo_df, plon_label = 'plon', plat_label='plat', age_label = 'age', 
                            window_length=window_length, time_step=time_step, max_age=max_age, min_age=min_age)
 
        RM['run'] = float(i)
        running_means_global = running_means_global.append(RM, ignore_index=True)

    running_means_global['plon'] = running_means_global['plon'].where(running_means_global['plon'] <= 180, running_means_global['plon'] - 360)
    return running_means_global

