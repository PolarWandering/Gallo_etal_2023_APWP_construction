## ***Jupyter Notebooks to accompany the paper:***
    
# ***Embracing uncertainty to resolve polar wander: a case study of Cenozoic North America***

### ***by L. C. Gallo<sup>1<sup>*** *(len.gallo@gmail.com)****, M. Domeier<sup>1</sup>, F. Sapienza<sup>2</sup>, ETAL***

<!-- B. Vaes<sup>3</sup>, N. Swanson-Hysell<sup>4</sup>, M. Arnaould<sup>5</sup>, A. Eyster,  D. Gürer<sup>7</sup>, A. Kiraly<sup>1</sup>, B. Robert<sup>8</sup>, T. Rolf<sup>1</sup>, G. Shephard<sup>1</sup>,  A. van der Boon<sup>1</sup>, L. Wu and Y. Zhang<sup>4</sup>***-->


*(1) Centre for Earth Evolution and Dynamics, University of Oslo, Norway.*
*(2) Department of Statistics, University of California Berkeley, United States.*
*(3) Department of Earth Sciences, Faculty of Geosciences, University of California Berkeley, United States.*
*(4) Department of Earth Sciences, Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands.*
*(5) Laboratoire de Géologie de Lyon - Terre, Planètes, Environnement, University Lyon 1 Claude Bernard, Lyon, France*
*(6) *
*(7) Research School of Earth Sciences, Australian National University, Canberra, Australian Capital Territory, Australia.*
*(8) Institut de Physique du Globe de Paris, Université de Paris, Paris, France.*

Available at: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PolarWandering/Gallo_etal_2023_APWP_construction/HEAD)

This manuscript is under review at *Geophysical Research Letters*     

## code

This folder contains Jupyter notebooks developed using a Python 3 kernel that generate the figures and conduct data analysis associated with the study. 


#### 0_Preprocessing.ipynb
    
This notebook preprocesses site level data from the VGP database compilation, checks for missing data, compares reported and calculated pole positions, conducts associated statistical tests, and makes visualizations of data at the study level.    


#### 1_Compile.ipynb
    
This notebook iterates through each paleomagnetic record (datasheet) of the vgp database, extracts data meeting user-specified criteria, and appends them to a new dataframe for later processing (to generate an APWP).    

#### MC_uncertainty_propagation.ipynb
    
This notebook illustrates the implementation of Section 2.2 in the paper.
    

#### 3_Comparisons.ipynb
    
    

#### Window_width_Optimization.ipynb    

In this following Notebook we find the optimal averaging window size, that is, an optimal balance between overfitting and underfitting of the Moving averages. We achieve this through a modified version of the leasts squares cross-validation technique