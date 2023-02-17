import numpy as np
import pandas as pd
from .utils import cartesian2spherical, spherical2cartesian, GCD_cartesian, PD


class quantiles:
    '''
    class to generate quantiles.
    note: input fro longitudes should live in [-180,180]
    '''
    
    def __init__(self, df, xlabel, ylabel):        
        self.X = np.sort(df[xlabel].unique().transpose())        
        self.Y = df.groupby(xlabel)[ylabel]
        
        self.q5 = self.Y.quantile(.05).to_numpy()
        self.q16 = self.Y.quantile(.16).to_numpy()
        self.q25 = self.Y.quantile(.25).to_numpy()
        self.q50 = self.Y.quantile(.50).to_numpy()
        self.q75 = self.Y.quantile(.75).to_numpy()
        self.q84 = self.Y.quantile(.84).to_numpy()
        self.q95 = self.Y.quantile(.95).to_numpy()
        self.mean = self.Y.mean().to_numpy()
    
class PC:
    
    '''
    Class to calculate Principal Componenents as a function of time from a time dependant ensemble of directions
    '''
    def __init__(self, df, xlabel, LatLabel, LonLabel):
        
        self.X = np.sort(df[xlabel].unique().transpose())
        self.df=df
        self.xlabel=xlabel
        self.LatLabel=LatLabel
        self.LonLabel=LonLabel
        self.groupby=df.groupby(xlabel)
        
    def PC(self, quantile): #This functions generates an array for each coordinates of the time-varying Principal Componente plus its circular empirical error
        '''
        Calculates the Prin. Comp. of the time varying ensemble of mean poles
        Then, given a quantile, it calculates the $\Theta_quantile$ metric, which is the maximum distance at a given confidence (circular)
        '''
        lats, lons, maxGCD95 = [], [], []
        for age, df_age in self.groupby:
            
            # calculate an array of directions in cartesian coordinates
            array = np.array([spherical2cartesian([ np.radians(i[self.LatLabel]),np.radians(i[self.LonLabel])]) for _,i in df_age.iterrows()])
            
            # calculate the principal component (maximum eigenvector) of the ensemble
            PrinComp=PD(array)
            
            # calculate the distance from each direction to the principal component 
            gcds = [GCD_cartesian(array[i],PrinComp) for i,_ in enumerate(array)]
            
            maxGCD95.append(np.degrees(np.percentile(gcds, quantile)))
            lats.append(np.degrees(cartesian2spherical(PrinComp))[0])
            lons.append(np.degrees(cartesian2spherical(PrinComp))[1])
                       
        lats[0]=-90 # set the present day field in -90
        
        return [np.array(lons),np.array(lats),np.array(maxGCD95)]
    
    def Lat_Lon_bounds(self, quantile): 
        '''
        Given a quantile, this function first discards the (100 - quantile) percent of most distanct vectors
        Returns a list with maximum and minimum latitudes and longitudes
        '''
        max_lat, min_lat, max_lon, min_lon = [], [], [], []
        for age, df_age in self.groupby:
            
            # calculate an array of directions in cartesian coordinates
            array = np.array([spherical2cartesian([ np.radians(i[self.LatLabel]),
                                                    np.radians(i[self.LonLabel])]) for _,i in df_age.iterrows()])        
            PrinComp=PD(array) # calculate the principal component (maximum eigenvector) of the array                      
            
            gcds = np.array([GCD_cartesian(array[i],PrinComp) for i,_ in enumerate(array)]) # calculate the distance from each direction to the principal component                                    
            
            # lats = np.array([i[self.LatLabel] for _,i in df_age.iterrows()])
            lats = df_age[self.LatLabel].values
            # lons = np.array([i[self.LonLabel] for _,i in df_age.iterrows()])  
            lons = df_age[self.LonLabel].values
            array = np.column_stack([lats, lons, gcds])  
            array_masked = array[array[:,2] <= np.percentile(array[:,2], quantile) , :] # discards quantile% of the most distanct veectors to the mean
            
            max_lat.append(array_masked[:,0].max())
            min_lat.append(array_masked[:,0].min())
            max_lon.append(array_masked[:,1].max())
            min_lon.append(array_masked[:,1].min())
            
        return [np.array(max_lon), np.array(min_lon), np.array(max_lat), np.array(min_lat)]    