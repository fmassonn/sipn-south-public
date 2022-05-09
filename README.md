# SIPN South scripts and data
Welcome to the public SIPN South directory. Here, you will find the data and scripts necessary to process Sea Ice Prediction Network South (SIPN South) analyses.

<img src="./doc/Logo2.png" width="30%">

# Quick start
Assuming that Git is installed:
1. In a terminal, do
```
git clone https://github.com/fmassonn/sipn-south-public.git
```
2. Go to the scripts directory to download the SIPN South forecast data:
```
cd sipn-south-public/scripts
./retrieve_data.bash   # Will fill the folder netcdf with data, takes a few minutes. 
                       # You can specify in that script the year(s) for which you want
                       # to retrieve data
```

3. Process the Python and NCL scripts 


# Structure of the project
```
data
  txt
    # contains data of total and regional sea ice areas, per year of forecast
  netcdf
    # empty at initialization, but can be populated by running the script 
    # ../scripts/retrieve_data.bash . The data will be downloaded from
    # the following link: https://nextcloud.cism.ucl.ac.be/s/gTL53xhjp4iQMM8 

scripts
  # contains scripts to retrieve observational data, forecast data, format them in 
  # CMIP-like format. Also contains scripts to produce figures.

figs
  # empty folder aimed at receiving the figures produced by the scripts

doc
  # contains the form with answer of participating groups to the SIPN South forecast
```
# Important notes

CNRM usually submits two forecasts, one raw and one bias-corrected. For 2019, the official forecast was the raw one, so there is another submission named CNRM_corrected. For 2020-2021, it's the opposite, so there is a CNRM_raw forecast in addition to CNRM



# References
Massonnet, F., P. Reid, J. L. Lieser, C. M. Bitz, J. Fyfe, W. Hobbs (2019). Assessment of February 2018-2019 sea-ice forecasts for the Southern Ocean. Technical
Note, Université catholique de Louvain, available at http://acecrc.org.au/sipn-south/

Massonnet, F., P. Reid, J. L. Lieser, C. M. Bitz, J. Fyfe, W. Hobbs (2018). Assessment of February 2018 sea-ice forecasts for the Southern Ocean. https://eprints.utas.edu.au/27184

# Primary Contact
francois.massonnet@uclouvain.be
