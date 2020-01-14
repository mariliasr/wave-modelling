# wave-modelling
Notebooks and scripts used for wave modelling with SWAN.

The script spec2d.py return the 2d-spectrum from Hs, Tp, pdir and spread in degrees. The script spec2d_era.py returns the 2d-spectrum from Hs, Tp, pdir and cossene power.

The notebooks have the following objectives:

- buoy_csiro_correl and buoy_era_correl:
    - Read .csv files containing buoy observations
    - Read netcdf file of CSIRO reanalyses containing mean wave period
    - Correlate ERA5 wave period to buoy observations
    - Plot the correlations as scatter plots

- create_specs and create_specs_era:
    - Calculate wave spectra from wave parameters of ERA5 partitions
    - Write it to file to use as input for SWAN

- create_wind_and_timepoints:
    - Writes wind and timepoint inputs for SWAN from netcdf files
    
- get_ww3:
    - Download wave watch 3 dataset from noaa repository

- netcdf_xyz:
    - writes gebco netcdf bathymetry to a .xyz file
    
