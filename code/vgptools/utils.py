import os
import numpy as np
import pandas as pd
from pmagpy import pmag, ipmag
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
import matplotlib.image as mpimg
import cartopy.crs as ccrs
import cartopy.feature as cfeature 
# from cartopy.io.img_tiles import StamenTerrain
import cartopy.io.img_tiles as cimgt
import seaborn as sns

def get_files_in_directory(path): 
    """
    Retrieves file names from a directory \
    \n\nInput: path = directory \
    \n\nOutput: list of subdirectories
    """

    # The last conditional here is in order to ignore the /DS_store file in macs 
    return [os.path.join(path, name) for name in os.listdir(path)
            if (os.path.isfile(os.path.join(path, name)) and (not name.startswith('.')))  ]


def cartesian2spherical(v):
    """
    Take an array of length corresponding to a 3-dimensional vector and returns a array of length 2
    with latitade and longitude
    """
    theta = np.arcsin(v[2])         #facu     theta = np.arccos(v[2]) 
    phi = np.arctan2(v[1], v[0])
        
    return [theta, phi]


def spherical2cartesian(v):
    """
    calculates the cartesian coordinates of an array (or list) representing a vector in spherical coordimates
    v[0] = theta - Latitude
    """
    
    x = np.cos(v[0]) * np.cos(v[1])  
    y = np.cos(v[0]) * np.sin(v[1])  
    z = np.sin(v[0])                 
    
    return [x,y,z]


def GCD_cartesian(cartesian1, cartesian2):
    '''
    Computes the great circle distance from the dot product of two vectors in cartesian coordinates 
    '''
    try:
        dot = np.dot(cartesian1, cartesian2)
    except :
        # print(f" Inside except block - GCD_cartesian: {cartesian1} {cartesian2}")
        return np.nan
    
    if np.isnan(dot).any(): return np.nan
    
    # Check if the dot product is not between -1 and 1
    if not (-1 <= dot <= 1):
        return np.nan
    gcd =  np.arccos(dot)    
    return gcd


def PD(ArrayXYZ):
    '''
    Calculates the principal component of an array of vectors in cartesian coordinates
    '''
    
    
    Eval, Evec = eigen_decomposition(ArrayXYZ)
    mean = fisher_mean(ArrayXYZ)
    
    if GCD_cartesian(Evec[:,0], mean) < np.pi /3:
        return Evec[:,0]
    elif GCD_cartesian(Evec[:,0]*-1, mean)< np.pi /3:
        return Evec[:,0] * -1
    elif GCD_cartesian(Evec[:,1], mean)< np.pi /3:
        return Evec[:,1]
    elif GCD_cartesian(Evec[:,1]*-1, mean)< np.pi /3:
        return Evec[:,1] *-1               
    elif GCD_cartesian(Evec[:,2], mean)< np.pi /3:
        return Evec[:,2]
    elif GCD_cartesian(Evec[:,2]*-1, mean)< np.pi /3:
        return Evec[:,2]*-1
    
def fisher_mean(ArrayXYZ):
    '''
    Numpy based calculation of the Fisher mean, computationally more efficient than pmagpy
    '''
    sumx = ArrayXYZ[:,0].sum()
    sumy = ArrayXYZ[:,1].sum()
    sumz = ArrayXYZ[:,2].sum()
    
    R = np.sqrt(sumx**2+sumy**2+sumz**2)

    meanxyz = [sumx/R, sumy/R, sumz/R]
    
    return meanxyz


def get_k(df_vgps): 
    '''
    fills rows in which there's missing the kappa value and it is possible to compute from alpha and n
    '''
    df_vgps['k'] = np.where(df_vgps['k'].isna(), ((140./df_vgps['alpha95'])**2)/(df_vgps['n']), df_vgps['k'])
    return df_vgps


def get_alpha95(df_vgps): 
    '''
    fills rows in which there's missing the a95 value and it is possible to compute from kappa and n
    '''
    df_vgps['alpha95'] = np.where(df_vgps['alpha95'].isna(), 140.0/np.sqrt(df_vgps['n'] * df_vgps['k']),df_vgps['alpha95'])
    return df_vgps

def get_angle(a, b, c):
    
    '''
    calculates the angle between three points in the sphere (b is pivot) using the law of cosines.
    arguments: a,b and c are cartesian points in the form [x,y,z]
    returns: an angle
    '''
    
    s = GCD_cartesian(a, c)  # an angular distance between the first and last pole of a track
    p1 = GCD_cartesian(a, b)     # is a distance between the rotation pole and the first pole of a track 
    p2 = GCD_cartesian(c, b)    # is an angular distance between the rotation pole and the last pole of a track. 
    
    angle = np.arccos((np.cos(s) - np.cos(p1) * np.cos(p2)) / (np.sin(p1) * np.sin(p2)))
    
    return np.degrees(angle)


def print_pole_statistics(reported_pole, vgp_mean, vgp_mean_recomputed):
    
    pole_summary = pd.DataFrame(columns = ['N', 'Plat', 'Plon', 'A95'],
                                 index=['Reported mean pole', 
                                        'Mean pole (calculated from VGPs)',
                                        'Mean pole (calculated from transformed directions)'])
    
    if len(reported_pole) > 0:
        pole_summary['Plon']['Reported mean pole'] = reported_pole.iloc[0]['Plon']
        pole_summary['Plat']['Reported mean pole'] = reported_pole.iloc[0]['Plat']
        pole_summary['A95']['Reported mean pole'] = reported_pole.iloc[0]['A95']
        pole_summary['N']['Reported mean pole'] = reported_pole.iloc[0]['N']
        
    if len(vgp_mean) > 2:
        pole_summary['Plon']['Mean pole (calculated from VGPs)'] = round(vgp_mean['dec'],1)
        pole_summary['Plat']['Mean pole (calculated from VGPs)'] = round(vgp_mean['inc'],1)
        pole_summary['A95']['Mean pole (calculated from VGPs)'] = round(vgp_mean['alpha95'],1)
        pole_summary['N']['Mean pole (calculated from VGPs)'] = round(vgp_mean['n'],1)
        
    if len(vgp_mean_recomputed) > 2:
        pole_summary['Plon']['Mean pole (calculated from transformed directions)'] = round(vgp_mean_recomputed['dec'],1)
        pole_summary['Plat']['Mean pole (calculated from transformed directions)'] = round(vgp_mean_recomputed['inc'],1)
        pole_summary['A95']['Mean pole (calculated from transformed directions)'] = round(vgp_mean_recomputed['alpha95'],1)
        pole_summary['N']['Mean pole (calculated from transformed directions)'] = round(vgp_mean_recomputed['n'],1)
        
    display(pole_summary)
        
    return pole_summary
       
    
def test_fishqq(merged):
    print('')
    print('Fisher distribution Q-Q test on VGPs')
    print('------------------------------------')
    if len(merged) <= 10: print ('Not enough sites to conduct quantile-quantile test')
    else:                              
        plt.figure(1,figsize=(7,3))
        ipmag.fishqq(di_block=merged)               
        plt.show()

        
def statistical_tests(di_mode1, di_mode2, vgps, 
                      study_number, pole_number, 
                      save_folder):
    """
    Conduct reversal tests on data of two polarities that have been split using the
    pmag.separate_directions() function. This function implements the bootstrap 
    reversal mean test (Tauxe et al., 1991), the Watson common mean reversal test 
    with classification (McFadden and McElhinny, 1990), and the Bayesian reversal 
    test (Heslop & Roberts, 2018). It generates a Dataframe that summarizes the test
    results.
    
    Parameters
    ----------
    di_mode1 : a di_block of directions [dec, inc] in one polarity
    di_mode2 : a di_block of directions [dec, inc] in the opposite polarity
    vgps : a di_block of vgp [lon, lat] group in single polarity
    
    """
    
    stat_test_results = pd.DataFrame(columns = ['result'],
                                         index=['Bootstrap reversal test', 
                                                'Parametric reversal test',
                                                'Bayesian reversal test',
                                                'Fisher Q-Q test'])
    
    if len(di_mode1) == 0 or len(di_mode2) == 0: 
        print (' ')
        print ('Only one polarity; cannot conduct reversal tests')
        stat_test_results['result']['Bootstrap reversal test'] = 'Only one polarity'
        stat_test_results['result']['Parametric reversal test'] = 'Only one polarity'
        stat_test_results['result']['Bayesian reversal test'] = 'Only one polarity'
        
    elif len(di_mode1) < 3 or len(di_mode2) < 3: 
        print ('')
        print ('Not enough sites from one (or both) polarity populations to conduct Bootstrap and parametric reversal tests')
        print ('')
        stat_test_results['result']['Bootstrap reversal test'] = 'Too few sites for test'
        stat_test_results['result']['Parametric reversal test'] = 'Too few sites for test'
        
        print('Bayesian reversal test (Heslop & Roberts, 2018)')
        print('------------------------------------------------')
        B, P, support  = ipmag.common_mean_bayes(np.array(di_mode1), np.array(di_mode2))
        print('')
        stat_test_results['result']['Bayesian reversal test'] = support
    else: 
        flipped2 = ipmag.do_flip(di_block=di_mode2, unit_vector=False)
        
        print('Bootstrap reversal test (Tauxe et al., 1991)')
        print('------------------------------------------------')
        bootstrap_result = ipmag.common_mean_bootstrap(di_mode1, flipped2, save=True,
                                                       save_folder=save_folder,
                                                       fmt='png') 
        print('')
        if bootstrap_result==1:
            stat_test_results['result']['Bootstrap reversal test'] = 'Pass'
        if bootstrap_result==0:
            stat_test_results['result']['Bootstrap reversal test'] = 'Fail'
        
        print('Parametric reversal test with classification (McFadden and McElhinny, 1990)')
        print('------------------------------------------------')
        result, angle, critical_angle, classification = ipmag.common_mean_watson(di_mode1, flipped2)
        print('')
        if result==0:
            stat_test_results['result']['Parametric reversal test'] = 'Fail (angle ' + str(round(angle,1)) + 'º above ' + str(round(critical_angle,1)) + 'º critical angle)'
        if result==1:
            stat_test_results['result']['Parametric reversal test'] = 'Pass (angle ' + str(round(angle,1)) + 'º below ' + str(round(critical_angle,1)) + 'º critical angle); ' + classification + ' classification' 
    
        print('Bayesian reversal test (Heslop & Roberts, 2018)')
        print('------------------------------------------------')
        B, P, support  = ipmag.common_mean_bayes(np.array(di_mode1), np.array(di_mode2), reversal_test=True)
        print('')
        stat_test_results['result']['Bayesian reversal test'] = support
           
    print('Fisher distribution Q-Q test (Fisher et al., 1987)')
    print('------------------------------------------------')
    if len(vgps) <= 10: 
        print('Not enough sites to conduct Q-Q test')
        stat_test_results['result']['Fisher Q-Q test'] = 'Too few sites'
    else:
        QQ_dict = ipmag.fishqq(di_block = vgps,save=True,
                               save_folder=save_folder,
                               fmt='png')
        plt.show()
        print(QQ_dict['Test_result'])
        stat_test_results['result']['Fisher Q-Q test'] = QQ_dict['Test_result']
    
    display(stat_test_results)

    return stat_test_results 

        
def invert_polarity(mode1, mode2):
    '''
    Direcctions are flipped to a single polarity and pooled into a single variable (merged)
    If mode2 is empty, the pooled group is equal to mode1
    '''
       
    if mode2.size == 0: 
        merged = mode1
    elif mode1.size == 0: 
        merged = np.array(ipmag.do_flip(di_block=mode2))[:,:2]
    else:
        flipped2 = np.array(ipmag.do_flip(di_block=mode2))[:,:2]
        #flipped2 = np.delete(np.array(ipmag.do_flip(di_block=mode2)), -1, axis=1)
        merged = np.concatenate((mode1, flipped2))
    
    return merged


def scale_bar(ax, length=None, location=(0.5, 0.05), 
              linewidth=3, zorder=100):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth,zorder=zorder)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom',zorder=zorder)
    

def summary_figure(dataframe, vgp_mode1, vgp_mode2,
                   reported_pole, vgp_mean, pole_folder):
    
    fig = plt.figure(figsize=(8, 12))

    NA_extent = [-120, -65, 15, 65]
    central_lon = np.mean(NA_extent[:2])
    central_lat = np.mean(NA_extent[2:]) 
    NA_proj = ccrs.AlbersEqualArea(central_lon, central_lat)

    lat_buffer = 0.05 * (max(dataframe['slat'])-min(dataframe['slat']))
    study_extent = [min(dataframe['slon']),max(dataframe['slon']),
                    min(dataframe['slat'])-lat_buffer,
                    max(dataframe['slat'])+lat_buffer]
    study_proj = ccrs.Mercator()

    vgp_proj = ccrs.Orthographic(
            central_longitude=0, central_latitude=90)

    gs = GridSpec(5, 2, 
                  width_ratios=[1, 1], 
                  height_ratios=[1, 1, 0.8, 1.1,0.8])

    ax1 = fig.add_subplot(gs[0], projection=NA_proj)
    ax1.set_extent(NA_extent)
    ax1.stock_img()
    ax1.add_feature(cfeature.BORDERS)
    ax1.gridlines(draw_labels=False)
    ax1.add_feature(cfeature.LAND)
    ax1.add_feature(cfeature.COASTLINE)
    ax1.scatter(x = dataframe['slon'], y = dataframe['slat'], color='r', s=20, 
                marker='*', transform=ccrs.PlateCarree(),zorder=100,label='site locations')
    plt.legend()
    
    ax2 = fig.add_subplot(gs[1], projection=study_proj)
    ax2.set_extent(study_extent)
    
    #tiler = StamenTerrain()
    # tiler = cimgt.Stamen(style = 'terrain')
    tiler = cimgt.GoogleTiles(style='terrain')
    
    mercator = tiler.crs
    ax2.add_image(tiler, 6)

    gl = ax2.gridlines(draw_labels=True)
    gl.xlabels_bottom = False
    gl.ylabels_left = False

    ax2.scatter(x = dataframe['slon'], y = dataframe['slat'], color='r', s=20, 
                marker='*', transform=ccrs.PlateCarree(),zorder=100,label='site locations')
    scale_bar(ax2, 10)

    ax3 = fig.add_subplot(gs[2])
    ipmag.plot_net()
    ipmag.plot_di(dataframe['dec'].tolist(), dataframe['inc'].tolist(), 
                        label = "directions")
    ax3.set_title('site mean directions', 
                  y=-0.1, backgroundcolor= 'white')

    lat_grid=[-60.0, -30.0, 0.0, 30.0, 60.0]
    lon_grid=[-180.0, -150.0, -120.0, -90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0]
    ax4 = fig.add_subplot(gs[3], projection=vgp_proj)
    ax4.add_feature(cfeature.OCEAN, zorder=0, facecolor='lightblue')
    ax4.add_feature(cfeature.LAND, zorder=0, facecolor='bisque')
    ax4.gridlines(xlocs=lon_grid, ylocs=lat_grid, linewidth=1,
                 color='black', linestyle='dotted')

    mode2_flipped = ipmag.do_flip(di_block=vgp_mode2)

    if len(vgp_mode1)>0:
        if vgp_mode1[0][1]>0:
            ipmag.plot_vgp(ax4,di_block=vgp_mode1.tolist(), 
                           label = "VGPs (normal)", 
                           marker='s',color='blue')
            ipmag.plot_vgp(ax4,di_block=mode2_flipped, 
                           label = "VGPs (reverse)",
                           marker='d',color='red')
            plt.legend()
        else:
            ipmag.plot_vgp(ax4,di_block=mode2_flipped, 
                           label = "VGPs (normal)", 
                           marker='s',color='blue')
            ipmag.plot_vgp(ax4,di_block=vgp_mode1.tolist(), 
                           label = "VGPs (reverse)",
                           marker='d',color='red')
    elif len(vgp_mode2)>0:
        if vgp_mode2[0][1]>0: 
            ipmag.plot_vgp(ax4,di_block=mode2_flipped, 
                           label = "VGPs (normal)", 
                           marker='s',color='blue')
        elif vgp_mode2[0][1]<0:
            ipmag.plot_vgp(ax4,di_block=mode2_flipped, 
                           label = "VGPs (reverse)",
                           marker='d',color='red')
            
            
            
    if (np.isnan(reported_pole['Plon'].values[0]) == False) & (np.isnan(reported_pole['Plat'].values[0]) == False) & (np.isnan(reported_pole['A95'].values[0]) == False):
        pole_lat = reported_pole['Plat'].values[0]
        pole_lon = reported_pole['Plon'].values[0]
        pole_A95 = reported_pole['A95'].values[0]
        ipmag.plot_pole(ax4, pole_lon, pole_lat, pole_A95, 
                        label="reported pole", color='y')
    # if not len(vgp_mean) == 0:
    if len(vgp_mean) >2:
        ipmag.plot_pole(ax4, vgp_mean['dec'], vgp_mean['inc'], vgp_mean['alpha95'], 
                             label="recalculated pole", color='red')
    plt.legend(loc='lower center')    
    ax4.set_title('virtual geomagnetic poles and mean poles', 
                  y=-0.1, backgroundcolor= 'white')

    reversal_test_img_path = pole_folder + '/common_mean_bootstrap.png'
    if os.path.exists(reversal_test_img_path):
        reversal_test_img = mpimg.imread(reversal_test_img_path)
        ax5 = fig.add_subplot(gs[2, :])
        ax5.imshow(reversal_test_img)
        ax5.axis("off")

    QQ_test_img_path =  pole_folder + '/QQ_mode1.png'
    if os.path.exists(QQ_test_img_path):
        QQ_test_img = mpimg.imread(QQ_test_img_path)
        ax6 = fig.add_subplot(gs[7])
        ax6.imshow(QQ_test_img)
        ax6.axis("off")
    
    ax5 = fig.add_subplot(gs[6])
    ax5.set_title('Age distribution of VGPs')
    ax5.set_xlabel(r'Age (Ma)')
    
    if dataframe['mean_age'].any() == False:
        dataframe['mean_age']=(dataframe['min_age'] + dataframe['max_age'])/2
    
    if (dataframe['mean_age'].max() - dataframe['mean_age'].min())<20:
        ax5.set_xlim(dataframe['mean_age'].mean()-10, dataframe['mean_age'].mean()+10)
        
    sns.histplot(dataframe, x = 'mean_age', ax=ax5)

    plt.tight_layout()
    plt.savefig(pole_folder + '/pole_summary.png',dpi=450)


# def Plot_VgpsAndSites(df, pole, vgp_mean, mode1, mode2):
#     '''
#     Function to plot Site coordinates and directions in a stereonet along with the recomputed and reported pole  
#     '''    
    
#     fig = plt.figure(figsize=(12, 5))
    
#     extent = [-135, -55,5, 90]
#     central_lon = np.mean(extent[:2])
#     central_lat = np.mean(extent[2:])    
    
#     proj = ccrs.AlbersEqualArea(central_lon, central_lat)

#     #ax = plt.axes(projection=proj)
#     ax = fig.add_subplot(1,2,1, projection=proj)
#     ax.set_extent(extent)

#     # Put a background image on for nice sea rendering.
#     ax.stock_img()

#     ax.add_feature(cfeature.BORDERS)
#     ax.gridlines()
#     ax.add_feature(cfeature.LAND)
#     ax.add_feature(cfeature.COASTLINE)

        
#     #if not df[df['slat'].isna()].empty: 
#     ax.scatter(x = df['slon'], y = df['slat'], color='r', s=20, marker='*', transform=ccrs.PlateCarree())
    
#     ax2 = fig.add_subplot(1,2,2)    
#     ax2.set_title('VGPs and Paleomagnetic poles')
    
#     ax2 = ipmag.plot_net()
#     ax2 = ipmag.plot_di([x[0] for x in mode1], [x[1] for x in mode1], 
#                         label = "VGPs")
    
#     if len(mode2) != 0:
#         ax2 = ipmag.plot_di([x[0] for x in mode2], [x[1] for x in mode2],
#                             label = "VGPs")
            
#     if not pole.empty:
#         pole_lat = pole['Plat'].values[0]
#         pole_lon = pole['Plon'].values[0]
#         pole_A95 = pole['A95'].values[0]
#         ax2 = ipmag.plot_di_mean(pole_lon, pole_lat, pole_A95, 
#                                  label="Reported Pole", color='y')    
#     if not len(vgp_mean) == 0:
#         ax2 = ipmag.plot_di_mean(vgp_mean['dec'], vgp_mean['inc'], vgp_mean['alpha95'], 
#                              label="Recalculated Pole", color='red')
#     plt.legend(loc=3, fontsize=12)
#     plt.show()


def get_site_coordinates(D, I, Plat, Plon):
    '''
    The following function retrives Site coordinates from Dec/Inc and Plat/Plom
    NOTE! there are always two possible solutions so the outcome is a list in the form [[Slat1, Slon1],[Slat2, Slon2]]
    '''
    paleolat = np.degrees(np.arctan(0.5 * np.tan(np.radians(I))))
    colatitude = 90 - paleolat
    beta = np.degrees(np.arcsin((np.sin(np.radians(colatitude)) * np.sin(np.radians(D))) / (np.cos(np.radians(Plat)))))
    
    def guess(Slat):
        guess = np.arcsin(np.sin(np.radians(Slat)) * np.cos(np.radians(colatitude)) +
                     np.cos(np.radians(Slat)) * np.sin(np.radians(colatitude)) * np.cos(np.radians(D))) - np.radians(Plat)
        return np.degrees(guess)

    sol = optimize.root(guess, x0 = [-90,90],  method='hybr')
    
    res = []
    
    for i in sol.x:
                
        if np.cos(np.radians(colatitude)) > np.sin(np.radians(i)) * np.sin(np.radians(Plat)):            
            Slon = (Plon - beta) % 360.
        else:
            Slon = (Plon - 180 + beta) % 360.
        
        res.append([i,Slon])
    
    return res

def orientation_matrix(ArrayXYZ):
    '''input : np.array of x,y,z coordinates
    returns the orientation matrix
    '''
    X =[
        [sum(ArrayXYZ[:,0]*ArrayXYZ[:,0]),sum(ArrayXYZ[:,0]*ArrayXYZ[:,1]),sum(ArrayXYZ[:,0]*ArrayXYZ[:,2])],
        [sum(ArrayXYZ[:,0]*ArrayXYZ[:,1]),sum(ArrayXYZ[:,1]*ArrayXYZ[:,1]),sum(ArrayXYZ[:,1]*ArrayXYZ[:,2])],
        [sum(ArrayXYZ[:,0]*ArrayXYZ[:,2]),sum(ArrayXYZ[:,1]*ArrayXYZ[:,2]),sum(ArrayXYZ[:,2]*ArrayXYZ[:,2])]
        ]
    X = np.array(X) / len(ArrayXYZ)
    
    return X

def eigen_decomposition(ArrayXYZ):
    '''
    input: np.array of direction cosines
    outuput : watch out to the shape !!! => eigen_vectors[:,0] => represents the principal direction vector
    '''  
    Means = [ArrayXYZ[:,0].mean(),ArrayXYZ[:,1].mean(),ArrayXYZ[:,2].mean()]   

    X = orientation_matrix(ArrayXYZ)
    eigenValues, eigenVectors = np.linalg.eig(X)
    
    #the following block sorts the eigenvectros acrodnig to the eigenvalues
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
           
    return eigenValues, eigenVectors

def shape(ArrayXYZ):
    '''given an array of direction cosines [x,y,z], computes its eigen parameter and 
    returns the shape as [oblateness, prolateness,collinearity(K),coplanarity(M)]'''
    eigen_values = eigen_decomposition(ArrayXYZ)[0]
    
    O = np.log(eigen_values[1]/eigen_values[2]) # Oblateness
    P = np.log(eigen_values[0]/eigen_values[1]) # Prolateness
    K = np.log(eigen_values[0]/eigen_values[1])/np.log(eigen_values[1]/eigen_values[2]) #Collinearity
    M = np.log(eigen_values[0]/eigen_values[2]) #Coplanarity
    return [O, P, K, M] 

