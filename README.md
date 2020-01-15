# wave-modelling
Notebooks and functions used for wave modelling with SWAN. Includes data preparation, cleaning, munging, and visualisation.

The script spec2d.py returns the 2d-spectrum from Hs, Tp, pdir and wave spread. The script spec2d_era.py returns the 2d-spectrum from Hs, Tp, pdir and cossene power.

The notebooks have the following objectives:

- buoy_csiro_correl and buoy_era_correl:
    - Read .csv files containing buoy observations
    - Read netcdf file of CSIRO reanalyses containing mean wave period
    - Clean and prepare buoy data
    - Select ERA5 timeseries with the same coordinates and date/times of buoy observations
    - Correlate ERA5 wave period to buoy observations
    - Plot the correlations as scatter plots

- create_specs and create_specs_era:
    - Calculate wave spectra from wave parameters of ERA5 partitions
    - Write it to text file to use as input for SWAN

- create_wind_and_timepoints:
    - Writes wind and timepoint inputs for SWAN from netcdf files
    
- get_ww3:
    - Download wave watch 3 dataset from noaa repository via web scraping

- netcdf_xyz:
    - Writes gebco netcdf bathymetry to a .xyz file
- plot_specs and plot_specs_datasets:
    - plot_specs plots 2d-spectrum reconstructed from CSIRO wave partitions to a 3d plot
    - plot_specs_datasets plots 2d-spectra from three different sources: (1) downloaded directly from ERA5, reconstructed from ERA5 partitions, and reconstructed from CSIRO partitions. The plots are in polar coordinates. Also plots spectra-extracted Hs as timeseries for comparison purposes
    
- swan_out_processing:
    - Processes the output of SWAN. Extracts data from selected coordinates (wave buoy coordinates for correlation purposes)
    - Saves timeseries to numpy array
    
