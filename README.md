## ***Jupyter Notebooks to accompany the paper:***
    
# ***Embracing uncertainty to resolve polar wander: a case study of Cenozoic North America***

### ***by L. C. Gallo<sup>1<sup>*** *(len.gallo@gmail.com)****, M. Domeier<sup>1</sup>, F. Sapienza<sup>2</sup>, N. Swanson-Hysell<sup>3</sup>, B. Vaes<sup>4</sup>, Y. Zhang<sup>3</sup>, M. Arnaould<sup>5</sup>, A. Eyster<sup>6</sup>,  D. Gürer<sup>7</sup>, A. Kiraly<sup>1</sup>, B. Robert<sup>8</sup>, T. Rolf<sup>1</sup>, G. Shephard<sup>1</sup> and A. van der Boon<sup>1</sup>***

<!--   and ***-->


*(1) Centre for Earth Evolution and Dynamics, University of Oslo, Norway.*
*(2) Department of Statistics, University of California, Berkeley, United States.*
*(3) Department of Earth and Planetary Science, University of California, Berkeley, United States.*
*(4) Department of Earth Sciences, Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands.*
*(5) Laboratoire de Géologie de Lyon - Terre, Planètes, Environnement, University Lyon 1 Claude Bernard, Lyon, France*
*(6) Department of Geoscience, University of Wisconsin-Madison, Madison, WI USA*
*(7) Research School of Earth Sciences, Australian National University, Canberra, Australian Capital Territory, Australia.*
*(8) Institut de Physique du Globe de Paris, Université de Paris, Paris, France.*

Available at: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PolarWandering/Gallo_etal_2023_APWP_construction/HEAD)

[![DOI](https://zenodo.org/badge/599666881.svg)](https://zenodo.org/badge/latestdoi/599666881)

This manuscript is accepted at *Geophysical Research Letters*     

## code

The Jupyter notebooks stored in this folder were created using a Python 3 kernel and are utilized to produce the figures and perform data analysis related to the study.


#### 0_Preprocessing.ipynb
    
This notebook preprocesses site level data from the VGP database compilation, checks for missing data, compares reported and calculated pole positions, conducts associated statistical tests, and makes visualizations of data at the study level.    


#### 1_Compile.ipynb
    
This notebook iterates through each paleomagnetic record (datasheet) of the vgp database, extracts data meeting user-specified criteria, and appends them to a new dataframe for later processing (to generate an APWP).    

#### 2_MC_uncertainty_propagation.ipynb
    
This notebook illustrates the implementation of Section 2.2 in the paper.

#### 3_Comparisons.ipynb

This notebook compares APWPs resulting from moving averages on paleopoles (e.g. Torsvik et al., 2012), moving averages on simulated VGPs (Vaes et al., 2022), and moving averages with a weighted window on the actual VGPs (this study). Plot the results in Figure 4 of the paper.

#### Window_width_Optimization.ipynb    

In this Notebook we find the optimal averaging window size, that is, an optimal balance between overfitting and underfitting of the Moving averages. We achieve this through a modified version of the leasts squares cross-validation technique
