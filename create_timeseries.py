#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil

import pandas as pd
import xarray as xr

import pcraster as pcr
import virtualOS as vos

# ~ Steps:
# ~ 1. Read the wmo table/dataframe
# ~ 2. Loop through all stations
# ~    - Plot the point
# ~    - Evaluate the catchment area: error must be < 15%
# ~    - If error > 15%: Search within 3x3 window find a point with the error < 15%
# ~    - Create the time series file. 


# WMO table containing station info
wmo_station_table_file = "data/streamflow_stations.csv"
wmo_station_table = pd.read_csv(wmo_station_table_file)

# add the following columns to the dataframe 
wmo_station_table["model_lon"]      = pd.Series(dtype = "float")
wmo_station_table["model_lat"]      = pd.Series(dtype = "float")
wmo_station_table["model_area_km2"] = pd.Series(dtype = "float")
wmo_station_table["area_deviation"] = pd.Series(dtype = "float")


# ldd and cell areas used in the model
ldd_map_file = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/lddsound_05min_version_20210330.map"
area_m2_file = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/cdo_gridarea_clone_global_05min_correct_lats.nc.map"

# set the pcraster clone based on the ldd
pcr.setclone(lddmap_file)
ldd_map =
area_m2 = 

# calculate catchment areas in km2
model_area_km2 = pcr.catchmenttotal(area_m2, ldd_map) / (1000.*1000.)

# lon and lat coordinates of the model
xcoord = pcr.xcoordinate(pcr.defined(model_area_km2))
ycoord = pcr.ycoordinate(pcr.defined(model_area_km2))



for irow in range(len(wmo_station_table)):
    
    # get the station ids, coordinates, and catchment areas based on the table provided by WMO
    wmo_id       = wmo_station_table_file["id"][irow]
    wmo_lon      = wmo_station_table_file["lon"][irow]
    wmo_lat      = wmo_station_table_file["lat"][irow]
    wmo_area_km2 = wmo_station_table_file["area"][irow]
    
    # assign the station on the global map
    abs_lon_diff = pcr.abs(xcoord - pcr.scalar(wmo_lon))
    abs_lat_diff = pcr.abs(ycoord - pcr.scalar(wmo_lat))    
    wmo_id_point = pcr.ifthen(abs_lon_diff == pcr.mapminimum(abs_lon_diff), abs_lat_diff == pcr.mapminimum(abs_lat_diff), pcr.boolean(1.0))
    
    need_coordinate_adjustment = False
    if wmo_area_km2 > 0:
        error_area_km2 = pcr.mapmaximum(pcr.ifthen(wmo_id_point, pcr.abs(model_area_km2 - wmo_area_km2))) / wmo_area_km2
        if error_area_km2 > 0.15: need_coordinate_adjustment = True
    else:
        need_adjustment = False
    
    if need_adjustment:

        wmo_id_window = pcr.windowmaximum(wmo_id_point, pcr.clone.cellsize() * 3.0)
        
        wmo_id_catchment_area_abs_error = pcr.ifthen(wmo_id_window, pcr.abs(model_area_km2 - wmo_area_km2))
        
        wmo_id_catchment_area_abs_error_areaorder = pcr.areaorder(wmo_id_catchment_area_abs_error, wmo_id_window)  

        wmo_id_point_adjusted = pcr.ifthen(wmo_id_catchment_area_abs_error_areaorder == 1, pcr.boolean(1.0))
        
        wmo_id_point = wmo_id_point_adjusted
    
    # get the model lon and lat coordinates, as well as the model catchment area
    model_lon      = pcr.ifthen(wmo_id_point, xcoord) 
    model_lat      = pcr.ifthen(wmo_id_point, ycoord)
    model_area_km2 = pcr.ifthen(wmo_id_point, model_area_km2)
    
    # put them in the dataframe
    wmo_station_table["model_lon"][irow]      = model_lon     
    wmo_station_table["model_lat"][irow]      = model_lat     
    wmo_station_table["model_area_km2"][irow] = model_area_km2
    wmo_station_table["area_deviation"][irow] = (model_area_km2 - wmo_area_km2) / wmo_area_km2

    
    # get the timeseries (using xarray)
    # - go to the selected clone and output directorty
    # - pick the values in the netcdf file using xarray
    
    # write it to a file, see the following for the format
    # ~ (pcrglobwb_python3_v20250207) edwindan@tcn578.local.snellius.surf.nl:/scratch-shared/edwindan/data/wmo_2024$ head -n 5 wflowsbm_2040000010_dis_1991_2024.csv
    # ~ Date,Discharge
    # ~ 1991-01-01,172.31717
    # ~ 1991-01-02,169.47862
    # ~ 1991-01-03,158.54332
    # ~ 1991-01-04,165.56178

