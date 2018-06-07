# sipn-south-public
Data and scripts to process Sea Ice Prediction Network South (SIPN South) analyses.

# Structure of the project
```
data
  txt
    contains data of total and regional sea ice areas
  netcdf
    empty at initialization, but can be populated by running the script 
    ../scripts/retrieve_data.bash . The data will be downloaded from
    the following link: https://nextcloud.cism.ucl.ac.be/s/gTL53xhjp4iQMM8 

scripts
  contains scripts to retrieve observational data, forecast data, format them in 
  CMIP-like format. Also contains scripts to produce figures.

figs
  empty folder aimed at receiving the figures produced

doc
  contains the form with answer of participating groups to the SIPN South forecast
```
# Reference
F. Massonnet, P. Reid, J. L. Lieser, C. M. Bitz, J. Fyfe, W. Hobbs (2018). Assessment of February 2018 sea-ice forecasts for the Southern Ocean. Technical
Note, Université catholique de Louvain (2018), available at http://acecrc.org.au/sipnsouth/

# Primary Contact
francois.massonnet@uclouvain.be
