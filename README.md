## Introduction
This is the repository of code and data for "Finding Common Ground: Assessing 
the Co-location Potential of California’s Blue Food and Clean Energy Sectors", 
published in Ocean and Coastal Management (2025).

This code should be used in coordination with files available on Zenodo (https://zenodo.org/records/17545767).

Inqueries can be directed to Claire Gonzales (clairegonzales@ucsb.edu).

## Roadmap

- Sector values are wrangled and processed by sector (i.e., aquaculture, 
offshore wind energy, wave energy, wild capture fisheries catch)
- Sector-specific data are read into ArcGIS Pro as layers (ArcGIS model is 
availabe as a .py file on Zenodo)
- Sector data is transformed to be on a scale from 0-1
- Co-location Suitability score is evaluated by overlaying layers in ArcGIS Pro.