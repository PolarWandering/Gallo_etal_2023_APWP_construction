import numpy as np
import pandas as pd
from pmagpy import pmag, ipmag

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.geodesic import Geodesic
from shapely.geometry import Polygon
from vgptools.utils import eigen_decomposition, spherical2cartesian, cartesian2spherical, PD, GCD_cartesian


# def RM_stats(df, title, xlabel, ylabel, size=(7,3)):
      
#     fig, ax = plt.subplots(figsize=size)
#     ax2 = ax.twinx()
#     plt.title(title, fontsize = 14)

#     df['Circular Standart Deviation']=df['csd']
#     df['Number of Studies']=df['n_studies']
    
#     dfm = df[['age', 'A95', 'Number of Studies', 'Circular Standart Deviation']].melt('age', var_name='Statistic', value_name='value')

#     sns.lineplot(data  = dfm, x = dfm['age'], y = dfm['value'], hue = dfm['Statistic'], marker="o", ax=ax)    
#     sns.lineplot(data  = df, x = df['age'], y = df['k'], marker="o",  ax=ax2, color= "#C44D52")
    
#     ax2.yaxis.label.set_color('#C44D52')
#     ax.set_ylabel("Value", fontsize = 12)
#     ax2.set_ylabel("Kappa", fontsize = 13)
#     ax.set_xlabel("Mean age (Ma)", fontsize = 13)
#     plt.xlabel(xlabel, fontsize = 13)
    
def plot_pole_A95(Plat, Plon, A95, age, min_age, max_age, ax):
    """
    Before calling this function set the figure and axes as, for instance:
    
    fig = plt.figure(figsize=(20,10))
    proj = ccrs.Orthographic(central_longitude=0, central_latitude=-55) #30, -60
    ax = plt.axes(projection=proj)    
    ax.stock_img()
    ax.set_global() 
    ax.coastlines(linewidth=1, alpha=0.5)
    ax.gridlines(linewidth=1)

    """    
    cmap = mpl.cm.get_cmap('viridis')
    norm = mpl.colors.Normalize(min_age, max_age)
        
    ax.add_geometries([Polygon(Geodesic().circle(lon=Plon, lat=Plat, radius=A95*111139, n_samples=360, endpoint=True))], 
                      crs=ccrs.PlateCarree().as_geodetic(), 
                      facecolor='none', 
                      edgecolor=cmap(norm(age)), 
                      alpha=0.8, 
                      linewidth=1.5)
    plt.scatter(x = Plon, y = Plat, color = cmap(norm(age)),
                s=50, transform = ccrs.PlateCarree(), zorder=4,edgecolors='black', alpha = 1)
    

def plot_pole(Plat, Plon, age, min_age, max_age, ax):
    """
    Before calling this function set the figure and axes as, for instance:
    
    fig = plt.figure(figsize=(20,10))
    proj = ccrs.Orthographic(central_longitude=0, central_latitude=-55) #30, -60
    ax = plt.axes(projection=proj)    
    ax.stock_img()
    ax.set_global() 
    ax.coastlines(linewidth=1, alpha=0.5)
    ax.gridlines(linewidth=1)
    """    
    cmap = mpl.cm.get_cmap('viridis')
    norm = mpl.colors.Normalize(min_age, max_age)
        
    plt.scatter(x = Plon, y = Plat, color = cmap(norm(age)),
                s=50, transform = ccrs.PlateCarree(), zorder=4,edgecolors='black', alpha = 1)    
    
    

# def plot_APWP_df (df_apwp, extent, plot_A95s=True, connect_poles=False):
#     """
#     Functions to plot an APWP in the form of DF (as spited by the running mean function)
#     size scaling is N-dependant
#     Input:  a DataFrame with columns ['age','N','n_studies','k','A95','csd','plon','plat']
#     """
    
#     fig = plt.figure(figsize=(20,10))   
    
#     proj = ccrs.Orthographic(central_longitude=0, central_latitude=-55) #30, -60

#     ax = fig.add_subplot(1,2,1, projection=proj)    
#     ax.stock_img()
#     ax.coastlines(linewidth=1, alpha=0.5)
#     ax.gridlines(linewidth=1)
    
#     cmap = mpl.cm.get_cmap('viridis')
    
#     # plot the A95s
#     if plot_A95s:
#         norm = mpl.colors.Normalize(df_apwp["age"].min(), df_apwp["age"].max())
#         df_apwp['geom'] = df_apwp.apply(lambda row: Polygon(Geodesic().circle(lon=row["plon"], lat=row["plat"], radius=row["A95"]*111139, n_samples=360, endpoint=True)), axis=1)
#         for idx, row in df_apwp.iterrows():
#             ax.add_geometries([df_apwp['geom'][idx]], crs=ccrs.PlateCarree().as_geodetic(), facecolor='none', edgecolor=cmap(norm(df_apwp["age"][idx])), 
#                               alpha=0.6, linewidth=1.5)
#         df_apwp.drop(['geom'], axis=1)

    
#     sns.scatterplot(x = df_apwp["plon"], y = df_apwp["plat"], hue = df_apwp["age"], palette=cmap, size = df_apwp["N"], sizes=(50, 200),
#                 transform = ccrs.PlateCarree(), zorder=4)
#     if connect_poles:
#         plt.plot(df_apwp["plon"], df_apwp["plat"], transform = ccrs.Geodetic(),  linewidth=1.5)
        
#     if extent != 'global':
#         ax.set_extent(extent, crs = ccrs.PlateCarree())
        
#     handles, labels = ax.get_legend_handles_labels()
    
#     ax.legend(reversed(handles), reversed(labels))
    
# def plot_poles_and_APWP(extent, df_poles, df_apwp):
#     """
#     generates a side by side graphic showing the underlying paleomagnetic poles (in the left) along
#     with the Running Means path (in the right plot)
#     """
#     proj = ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0)
#     plt.figure(figsize=(12,10))

#     ax1 = plt.subplot(221, projection = proj)
#     ax1.patch.set_visible(False)
#     ax1.add_feature(cfeature.BORDERS)
#     ax1.add_feature(cfeature.LAND)
#     ax1.add_feature(cfeature.COASTLINE)
#     ax1.stock_img()
#     ax1.set_extent(extent, crs = ccrs.PlateCarree())
#     gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
#                       linewidth=0.8, color='gray', alpha=0.5, linestyle='-')
#     gl.ylabels_left = True

#     plt.title('Overview map of Paleomagnetic Poles', fontsize = 13)
#     for _, pole in df_poles.iterrows():
#         plot_pole(pole.Plat, pole.Plon, pole.mean_age, df_poles.mean_age.min(), df_poles.mean_age.max(), ax1)
#     plt.tight_layout()

#     ax2 = plt.subplot(222, projection=proj)
#     ax2.set_title('Classic Running Means path on Paleomagnetic Poles', fontsize = 13)
#     ax2.add_feature(cfeature.BORDERS)
#     gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
#                       linewidth=0.8, color='gray', alpha=0.5, linestyle='-')
#     gl.ylabels_left = True
#     ax2.add_feature(cfeature.LAND)
#     ax2.add_feature(cfeature.COASTLINE)
#     ax2.stock_img()
#     ax2.set_extent(extent, crs = ccrs.PlateCarree())

#     for _, pole in df_apwp.iterrows():
#         plot_pole_A95(pole.plat, pole.plon, pole.A95, pole.age, df_apwp.age.min(), df_apwp.age.max(), ax2)

#     plt.plot(df_apwp["plon"], df_apwp["plat"], transform = ccrs.Geodetic(),  linewidth=1.5, color = "black")
#     plt.tight_layout()
    
#     s = plt.scatter(
#         df_apwp.plon,
#         df_apwp.plat,
#         c = df_apwp.age,
#         edgecolors= "black", marker = "o", s = 25,
#         cmap="viridis",
#         transform=ccrs.PlateCarree(),
#     )
#     plt.tight_layout()

#     plt.colorbar(s, fraction=0.035).set_label("Age (My)",fontsize=12)    
    
# def RM_APWP_lat_lon_A95 (df_apwp):
    
#     df_apwp["plon_"] = df_apwp.apply(lambda row: row.plon - 360 if row.plon > 180 else row.plon, axis =1)
    
#     fig, axes = plt.subplots(2, 1, sharex=True, figsize=(15,6))
#     fig.suptitle('Running Mean APWP', fontsize= 16)
#     axes[0].set_title('Latitude (°N)', fontsize=12)
#     axes[1].set_title('Longitude (°E)', fontsize=12)

#     # plot latitude
#     axes[0].errorbar(df_apwp["age"].to_list(), df_apwp["plat"].to_list(), yerr=df_apwp["A95"].to_list(),zorder=1) #, fmt="o")
#     axes[0].scatter(df_apwp["age"].to_list(), df_apwp["plat"].to_list(), edgecolors = "black",zorder=2)
#     axes[0].set_ylabel(r'Latitude (°N)', fontweight ='bold')
#     # axes[0].set_ylabel(-90,-75)
    
#     # plot longitude    
#     axes[1].errorbar(df_apwp["age"].to_list(), df_apwp["plon_"].to_list(), yerr=df_apwp["A95"].to_list(),zorder=1)#, fmt="o")
#     axes[1].scatter(df_apwp["age"].to_list(), df_apwp["plon_"].to_list(), edgecolors = "black",zorder=2)   
#     axes[1].set_xlabel(r'Age (Ma)', fontweight ='bold')
#     axes[1].set_ylabel(r'Longitude (°E)', fontweight ='bold')
    
    
#     del df_apwp["plon_"]

# def plot_pseudoVGPs_and_APWP(extent, df_vgps, df_apwp):
#     """
#     generates a side by side graphic showing the underlying paleomagnetic poles (in the left) along
#     with the Running Means path (in the right plot)
#     """
#     proj = ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0)
#     plt.figure(figsize=(12,10))

#     ax1 = plt.subplot(221, projection = proj)
#     ax1.patch.set_visible(False)
#     ax1.add_feature(cfeature.BORDERS)
#     ax1.add_feature(cfeature.LAND)
#     ax1.add_feature(cfeature.COASTLINE)
#     ax1.stock_img()
#     #ax1.set_extent(extent, crs = ccrs.PlateCarree())
#     gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
#                       linewidth=0.8, color='gray', alpha=0.5, linestyle='-')
#     gl.ylabels_left = True

#     plt.title('Parametrically sampled $pseudo$-VGPs (one realization)')
    
#     for _, pole in df_vgps.iterrows():
#         plot_pole(pole.Plat, pole.Plon, pole.mean_age, df_vgps.mean_age.min(), df_vgps.mean_age.max(), ax1)
#     plt.tight_layout()


#     ax2 = plt.subplot(222, projection=proj)
#     ax2.set_title('Running averages on $pseudo$-VGPs (one realization)')
#     ax2.add_feature(cfeature.BORDERS)
#     gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
#                       linewidth=0.8, color='gray', alpha=0.5, linestyle='-')
#     gl.ylabels_left = True
#     ax2.add_feature(cfeature.LAND)
#     ax2.add_feature(cfeature.COASTLINE)
#     ax2.stock_img()
#     ax2.set_extent(extent, crs = ccrs.PlateCarree())

#     for _, pole in df_apwp.iterrows():
#         plot_pole_A95(pole.plat, pole.plon, pole.A95, pole.age, df_apwp.age.min(), df_apwp.age.max(), ax2)

#     plt.plot(df_apwp["plon"], df_apwp["plat"], transform = ccrs.Geodetic(), linewidth=1.5, color = "black")
#     plt.tight_layout()
    
#     s = plt.scatter(
#         df_apwp.plon,
#         df_apwp.plat,
#         c = df_apwp.age,
#         edgecolors= "black", marker = "o", s = 25,
#         cmap="viridis",
#         transform=ccrs.PlateCarree(),
#     )
#     plt.tight_layout()

#     plt.colorbar(s, fraction=0.035).set_label("Age (My)")  
    
# def plot_VGPs_and_APWP(extent, df_vgps, df_apwp):
#     """
#     generates a side by side graphic showing the underlying VGPs along
#     with the Running Means path (in the right plot)
#     """
#     proj = ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0)
#     plt.figure(figsize=(12,10))

#     ax1 = plt.subplot(221, projection = proj)
#     ax1.patch.set_visible(False)
#     ax1.add_feature(cfeature.BORDERS)
#     ax1.add_feature(cfeature.LAND)
#     ax1.add_feature(cfeature.COASTLINE)
#     ax1.stock_img()
#     #ax1.set_extent(extent, crs = ccrs.PlateCarree())
#     gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
#                       linewidth=0.8, color='gray', alpha=0.5, linestyle='-')
#     gl.ylabels_left = True

#     plt.title('Overview map of recompiled VGPs')
    
#     for _, pole in df_vgps.iterrows():
#         plot_pole(pole.vgp_lat_SH, pole.vgp_lon_SH, pole.mean_age, df_vgps.mean_age.min(), df_vgps.mean_age.max(), ax1)
#     plt.tight_layout()


#     ax2 = plt.subplot(222, projection=proj)
#     ax2.set_title('Running Mean path on VGPs')
#     ax2.add_feature(cfeature.BORDERS)
#     gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
#                       linewidth=0.8, color='gray', alpha=0.5, linestyle='-')
#     gl.ylabels_left = True
#     ax2.add_feature(cfeature.LAND)
#     ax2.add_feature(cfeature.COASTLINE)
#     ax2.stock_img()
#     ax2.set_extent(extent, crs = ccrs.PlateCarree())

#     for _, pole in df_apwp.iterrows():
#         plot_pole_A95(pole.plat, pole.plon, pole.A95, pole.age, df_apwp.age.min(), df_apwp.age.max(), ax2)

#     plt.plot(df_apwp["plon"], df_apwp["plat"], transform = ccrs.Geodetic(), linewidth=1.5, color = "black")
#     plt.tight_layout()
    
#     s = plt.scatter(
#         df_apwp.plon,
#         df_apwp.plat,
#         c = df_apwp.age,
#         edgecolors= "black", marker = "o", s = 25,
#         cmap="viridis",
#         transform=ccrs.PlateCarree(),
#     )
#     plt.tight_layout()

#     plt.colorbar(s, fraction=0.035).set_label("Age (My)") 

# def plot_APWP_RM_ensemble(df, title, size=(10,3)):
#     '''
#     pass a df with the colection of bootstrapped means for each run.
#     '''
#     df['plon_east'] = df.apply(lambda row: row.plon - 360 if row.plon > 180 else row.plon, axis =1)
    
#     ensemble_lat = quantiles(df,"age","plat") # set quantiles of latitude groupedby age for visualization purposes
#     ensemble_lon = quantiles(df,"age","plon_east") # set quantiles of latitude groupedby age for visualization purposes
#     ensemble_PC = PC(df,"age","plat","plon_east") # set principal component for each window
    
#     fig, axes = plt.subplots(2, 1, sharex=True, figsize=size)
#     fig.suptitle(title, fontsize= 14, fontweight ='bold')
#     axes[0].set_ylabel(r'Latitude (°N)', fontweight ='bold', fontsize = 12)
#     axes[1].set_ylabel(r'Longitude (°E)', fontweight ='bold', fontsize = 12)
        
#     axes[1].set_xlabel(r'Age (Ma)', fontweight ='bold')

#     # plot latitude
    
#     for run, df_run in df.groupby('run'):
#         axes[0].plot(df_run.age.tolist(), df_run.plat.tolist(), color="#4F4F4F", zorder =0, linewidth=0.2)
#     axes[0].fill_between(ensemble_lat.X, ensemble_lat.q16,ensemble_lat.q84, color= "#C44D52", alpha=.40, zorder =1)
#     axes[0].plot(ensemble_PC.X, ensemble_lat.q16,color="#590E0E", zorder =3, linewidth=0.2)
#     axes[0].plot(ensemble_PC.X, ensemble_lat.q84,color="#590E0E", zorder =3, linewidth=0.2)
#     axes[0].scatter(ensemble_PC.X, ensemble_PC.PC()[1],color="#590E0E",edgecolors='black',zorder =2)
#     axes[0].plot(ensemble_PC.X, ensemble_PC.PC()[1],color="#590E0E", zorder =3)

#     # plot longitude
   
#     for run, df_run in df.groupby('run'):
#         axes[1].plot(df_run.age.tolist(), df_run.plon.tolist(), color="#4F4F4F", zorder =0, linewidth=0.2)
#     axes[1].plot(df_run.age.tolist(), df_run.plon.tolist(), color="#4F4F4F", zorder =0, linewidth=0.2, label = 'Moving Average realization')
#     axes[1].fill_between(ensemble_lon.X, ensemble_lon.q16,ensemble_lon.q84, color= "#C44D52", alpha=.40, zorder =1, label="Ensemble 0.16-0.84 percentiles ")
#     axes[1].plot(ensemble_PC.X, ensemble_lon.q16,color="#590E0E", zorder =3, linewidth=0.2)
#     axes[1].plot(ensemble_PC.X, ensemble_lon.q84,color="#590E0E", zorder =3, linewidth=0.2)
#     axes[1].scatter(ensemble_PC.X, ensemble_PC.PC()[0],color="#590E0E",edgecolors='black',zorder =2, label = 'Principal component of the age ensemble')
#     axes[1].plot(ensemble_PC.X, ensemble_PC.PC()[0],color="#590E0E", zorder =3)
    
    
#     plt.legend(loc="upper left")
#     df = df.drop(['plon_east'], axis=1)
    
    
# class quantiles:
#     '''
#     class to generate quantiles.
#     note: input fro longitudes should live in [-180,180]
#     '''
    
#     def __init__(self, df, xlabel, ylabel):        
#         self.X = df[xlabel].unique().transpose()        
#         self.Y = df.groupby(xlabel)[ylabel]
        
#         self.q5 = self.Y.quantile(.05).to_numpy()
#         self.q16 = self.Y.quantile(.16).to_numpy()
#         self.q25 = self.Y.quantile(.25).to_numpy()
#         self.q50 = self.Y.quantile(.50).to_numpy()
#         self.q75 = self.Y.quantile(.75).to_numpy()
#         self.q84 = self.Y.quantile(.84).to_numpy()
#         self.q95 = self.Y.quantile(.95).to_numpy()
#         self.mean = self.Y.mean().to_numpy()
    
# class PC:
    
#     '''
#     Class to calculate PCs as a function of time from a time dependant ensemble of directions
#     '''
#     def __init__(self, df, xlabel, LatLabel, LonLabel):
        
#         self.X = df[xlabel].unique().transpose()
#         self.df=df
#         self.xlabel=xlabel
#         self.LatLabel=LatLabel
#         self.LonLabel=LonLabel
#         self.groupby=df.groupby(xlabel)
        
#     def PC(self):
#         lats, lons, maxGCD95 = [], [], []
#         for age, df_age in self.groupby:
#             array = np.array([spherical2cartesian([ np.radians(i[self.LatLabel]),np.radians(i[self.LonLabel])]) for _,i in df_age.iterrows()])
            
#             PrinComp=PD(array)
#             gcds = [GCD_cartesian(array[i],PrinComp) for i,_ in enumerate(array)]
            
#             maxGCD95.append(np.degrees(np.percentile(gcds, 5)))
#             lats.append(np.degrees(cartesian2spherical(PrinComp))[0])
#             lons.append(np.degrees(cartesian2spherical(PrinComp))[1])
                       
#         lats[0]=-90 # set the present day field in -90
        
#         return [np.array(lons),np.array(lats),np.array(maxGCD95)]
    
    

