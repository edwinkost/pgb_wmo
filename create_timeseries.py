#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil

import numpy as np
import pandas as pd
import xarray as xr

import pcraster as pcr
import virtualOS as vos

# ~ pd.options.mode.copy_on_write = True

# ~ Steps:
# ~ 1. Read the wmo table/dataframe
# ~ 2. Loop through all stations
# ~    - Plot the point
# ~    - Evaluate the catchment area: error must be < 15%
# ~    - If error > 15%: Search within 3x3 window find a point with the error < 15%
# ~    - Create the time series file. 

# model output folder
model_output_folder = "/scratch-shared/edwin/pcrglobwb_wmo_run/v20250417"

# WMO table containing station info
wmo_station_table_file = "data/streamflow_stations.csv"
wmo_station_table_ori = pd.read_csv(wmo_station_table_file)
# - sort by upstream area, small first, NaN at the bottom
wmo_station_table = wmo_station_table_ori.sort_values("area")

# add the following columns to the dataframe 
wmo_station_table["model_lon"]      = pd.Series(dtype = "float")
wmo_station_table["model_lat"]      = pd.Series(dtype = "float")
wmo_station_table["model_area_km2"] = pd.Series(dtype = "float")
wmo_station_table["area_deviation"] = pd.Series(dtype = "float")

# output folder for the time series
csv_output_folder = "/scratch-shared/edwindan/pcrglobwb_for_wmo_timeseries_v20250417/2nd_try/"
# ~ csv_output_folder = "/scratch-shared/edwindan/pcrglobwb_for_wmo_timeseries_v20250417_test/"


# ldd and cell areas used in the model
ldd_map_file = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/lddsound_05min_version_20210330.map"
area_m2_file = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/cdo_gridarea_clone_global_05min_correct_lats.nc.map"

# set the pcraster clone based on the ldd
pcr.setclone(ldd_map_file)
ldd_map = pcr.readmap(ldd_map_file)
area_m2 = pcr.readmap(area_m2_file)

# mask_classes used during the parallelization
mask_parallel_run_file = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/cloneMaps/global_parallelization/mask5minFromTop.map"
mask = pcr.readmap(mask_parallel_run_file)

# calculate catchment areas in km2
model_area_km2 = pcr.catchmenttotal(area_m2, ldd_map) / (1000.*1000.)
# - using only cell within the mask
model_area_km2 = pcr.ifthen(pcr.defined(mask), model_area_km2)

model_area_km2 = pcr.cover(model_area_km2, 0.0)

# lon and lat coordinates of the model
xcoord = pcr.xcoordinate(pcr.defined(model_area_km2))
ycoord = pcr.ycoordinate(pcr.defined(model_area_km2))


for irow in range(len(wmo_station_table)):
# ~ for irow in range(20):

    # get the station ids, coordinates, and catchment areas based on the table provided by WMO
    wmo_id       = wmo_station_table["id"][irow]
    wmo_lon      = wmo_station_table["lon"][irow]
    wmo_lat      = wmo_station_table["lat"][irow]
    wmo_area_km2 = wmo_station_table["area"][irow]
    
    print(wmo_id)
    print(wmo_lon)
    print(wmo_lat)
    print(wmo_area_km2)    
    
    # assign the station on the global map
    abs_lon_diff = pcr.abs(xcoord - pcr.scalar(wmo_lon))
    # ~ pcr.aguila(abs_lon_diff)
    abs_lat_diff = pcr.abs(ycoord - pcr.scalar(wmo_lat))    
    # ~ pcr.aguila(abs_lat_diff)
    wmo_id_point = pcr.ifthen(abs_lon_diff == pcr.mapminimum(abs_lon_diff), pcr.ifthen(abs_lat_diff == pcr.mapminimum(abs_lat_diff), pcr.boolean(1.0)))
    # ~ pcr.aguila(wmo_id_point)
    
    # check whether coordinates must be adjusted or not, do this only for stations with both of their values of wmo_area_km2 and model_area_km2 defined
    need_adjustment = False
    if wmo_area_km2 > 0:
        
        error_area_km2, valid = pcr.cellvalue(pcr.mapminimum(pcr.ifthen(wmo_id_point, pcr.abs(model_area_km2 - wmo_area_km2))) / wmo_area_km2, 1)
        print(error_area_km2)
        if abs(error_area_km2) > 0.20: need_adjustment = True
    else:
        need_adjustment = False
    
    print(need_adjustment)
    
    # if coordinate adjustment is needed, find - within the surrounding cells - the cell with the model catchment area closest to wmo catchment area
    if need_adjustment:

        # define the window
        # ~ wmo_id_window = pcr.boolean(pcr.windowmaximum(pcr.scalar(wmo_id_point), pcr.clone().cellSize() * 10.0))
        wmo_id_window = pcr.boolean(pcr.windowmaximum(pcr.scalar(wmo_id_point), 0.75))
        # ~ pcr.aguila(wmo_id_window)
        
        # calculate the catchment area error
        wmo_id_catchment_area_abs_error = pcr.ifthen(wmo_id_window, pcr.abs(model_area_km2 - wmo_area_km2))
        
        # sort from the lowest error
        wmo_id_catchment_area_abs_error_areaorder = pcr.areaorder(wmo_id_catchment_area_abs_error, wmo_id_window)  

        # pick the cell with the lost error and use it
        wmo_id_point_adjusted = pcr.ifthen(wmo_id_catchment_area_abs_error_areaorder == 1, pcr.boolean(1.0))
        wmo_id_point = wmo_id_point_adjusted
    
    # get the model lon and lat coordinates, as well as the model catchment area
    model_lon, valid                    = pcr.cellvalue(pcr.mapminimum(pcr.ifthen(wmo_id_point, xcoord)        ), 1)
    model_lat, valid                    = pcr.cellvalue(pcr.mapminimum(pcr.ifthen(wmo_id_point, ycoord)        ), 1)
    model_area_km2_this_station, valid  = pcr.cellvalue(pcr.mapminimum(pcr.ifthen(wmo_id_point, model_area_km2)), 1)
    
    # calculate area_deviation
    area_deviation = (model_area_km2_this_station - wmo_area_km2) / wmo_area_km2
    
    # - check where the mask for this
    mask_for_this_station, valid = pcr.cellvalue(pcr.mapmaximum(pcr.scalar(pcr.ifthen(wmo_id_point, mask))), 1)

    if abs(area_deviation) < 100. and mask_for_this_station < 100.:

    # ~ use_all = True
    # ~ if use_all:

        # ~ # put them in the dataframe
        # ~ wmo_station_table["model_lon"].loc[irow]      = model_lon     
        # ~ wmo_station_table["model_lat"].loc[irow]      = model_lat     
        # ~ wmo_station_table["model_area_km2"].loc[irow] = model_area_km2
        # ~ wmo_station_table["area_deviation"].loc[irow] = (model_area_km2 - wmo_area_km2) / wmo_area_km2
        # ~ print(wmo_station_table["area_deviation"][irow])
	    
        # put them in the dataframe
        wmo_station_table.loc[irow, "model_lon"]      = model_lon     
        wmo_station_table.loc[irow, "model_lat"]      = model_lat     
        wmo_station_table.loc[irow, "model_area_km2"] = model_area_km2_this_station
        wmo_station_table.loc[irow, "area_deviation"] = area_deviation
        print(model_lon)
        print(model_lat)    
        print(model_area_km2_this_station)    
        print(area_deviation)


        # - go to the selected clone and netcdf file
        mask_for_this_station, valid = pcr.cellvalue(pcr.mapmaximum(pcr.scalar(pcr.ifthen(wmo_id_point, mask))), 1)
        mask_code   = 'M%07d' %(mask_for_this_station)
        netcdf_file = model_output_folder + "/" + mask_code + "/netcdf/discharge_dailyTot_output.nc"     
        print(netcdf_file)
        
        # - pick the timeseris in the netcdf file using xarray
        discharge_xr          = xr.open_dataset(netcdf_file)
        discharge_time_series = discharge_xr.sel(lon = model_lon, lat = model_lat, method = 'nearest')
        
        # - using 1991-2024 only
        discharge_time_series = discharge_time_series.sel(time=slice("1991-01-01", "2024-12-31"))
	    
        # - create the data frame
        df = pd.DataFrame({\
        'date': np.asarray(discharge_time_series["time"]),\
        'discharge_m3persecond': np.asarray(discharge_time_series["discharge"])})
        print(df)
        
	    
        # ~ # write it to a file, see the following for the format
        # ~ (pcrglobwb_python3_v20250207) edwindan@tcn578.local.snellius.surf.nl:/scratch-shared/edwindan/data/wmo_2024$ head -n 5 wflowsbm_2040000010_dis_1991_2024.csv
        # ~ Date,Discharge
        # ~ 1991-01-01,172.31717
        # ~ 1991-01-02,169.47862
        # ~ 1991-01-03,158.54332
        # ~ 1991-01-04,165.56178
	    
	    
        # - write the data frame to csv
        csv_filename = csv_output_folder + "/" + "pcrglobwb_" + str(wmo_id) + "_discharge_1991_2024.csv"  
        df.to_csv(csv_filename, index = False)
    
    
# write the station list to a csv file
csv_station_list = csv_output_folder + "/" + "_station_list.csv"  
wmo_station_table.to_csv(csv_station_list, index = False)


