# Sensor Placements in Gaussian Processes Using Column Subset Selection

This repository contains the MATLAB implementation for optimal sensor placement in Gaussian Processes (GPs) using exchange-based algorithms and column subset selection. 

## Overview
The code provides tools to determine optimal sensor locations based on D-optimality criteria, with specific applications to 1D droplet dynamics and 2D Sea Surface Temperature (SST) datasets.

## Data Availability
The data used in these experiments can be found at the following link:
[https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html)

make sure to download the weekly mean SST and the corresponding masks data file

## Requirements
- MATLAB (Recommended R2021a or later)
- (Optional) H2Pack-Matlab for accelerated kernel computations (included in subdirectories)

## Getting Started

To run the Sea Surface Temperature (SST) example:
- download the data
- download H2Pack if using Random Projection Nystrom approximation (you can stick with Cholesky for now)

1. **Setup the Environment:** Run the setup script to load data and initialize parameters.
   ```matlab
   setup_sst
2. **Run the Main Script:** Execute the main analysis script.
   ```matlab
   script_sst
