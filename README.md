## Introduction
This is the repository of code and data for "Finding Common Ground: Assessing 
the Co-location Potential of California’s Blue Food and Clean Energy Sectors", 
published in Ocean and Coastal Management (2025).

This code should be used in coordination with files available on Zenodo (https://zenodo.org/records/17545767).

Inqueries can be directed to Claire Gonzales (clairegonzales AT ucsb.edu).

## Road Map

- Sector values are wrangled and processed by sector 
  - **Aquaculture** was modeled using a MATLAB file available on Zenodo. Environmental input data are required for this model, which are available as 
  netCDF files at `data/bivalve_model` and processed at `code/netCDFfile_testing`. After processing, environmental inputs were read into the MATLAB code 
  (see Zenodo) for production modeling.
  - **Fisheries** catch data was processed using `code/catch_wrangling`. The processed version can be found at `data/catchdata_ouput`.
    - Additional information about this data can be found at: https://github.com/cfree14/wcfish
  - **Offshore wind** data was contributed and is available on Zenodo as a .tiff file for processing in ArcGIS Pro.
  - **Wave energy** is pulled using MHKiT API (https://mhkit-software.github.io/MHKiT/mhkit-matlab/api.wave.html). Wave energy output can be found in 
  `random_points_wgs1984` and processed using `code/wave-power-data`.
- Sector-specific data are read into ArcGIS Pro as layers (ArcGIS model is 
availabe as a .py file on Zenodo)
- Sector data is transformed to be on a scale from 0-1
- Co-location Suitability score is evaluated by overlaying layers in ArcGIS Pro.
  - Scoring can be viewed by downloading the ArcGIS model from Zenodo. Detailed methodology is described in the full manuscript.