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
model_output_folder = "/scratch/depfg/sutan101/wmo_2025/pcrglobwb/source_from_bram/"

# WMO table containing station info
wmo_station_table_file = "data/streamflow_stations.csv"
wmo_station_table_ori = pd.read_csv(wmo_station_table_file)
# - sort by upstream area
# ~ wmo_station_table = wmo_station_table_ori.sort_values("area", ignore_index = True, na_position = "first")
wmo_station_table = wmo_station_table_ori.sort_values("area", ascending = False, ignore_index = True, na_position = "last")

# add the following columns to the dataframe 
wmo_station_table["model_lon"]      = pd.Series(dtype = "float")
wmo_station_table["model_lat"]      = pd.Series(dtype = "float")
wmo_station_table["model_area_km2"] = pd.Series(dtype = "float")
wmo_station_table["area_deviation_in_percent"] = pd.Series(dtype = "float")

# output folder for the time series
#~ csv_output_folder = "/scratch-shared/edwindan/pcrglobwb_for_wmo_timeseries_v20250417/2nd_try/"
# ~ csv_output_folder = "/scratch-shared/edwindan/pcrglobwb_for_wmo_timeseries_v20250417_test/"
csv_output_folder = "/scratch/depfg/sutan101/wmo_2025/pcrglobwb/pcrglobwb_for_wmo_timeseries_v20260409/test/"

#~ sutan101@node044.cluster:/scratch/depfg/sutan101/wmo_2025/pcrglobwb$ ls -lah
#~ total 37K
#~ drwxr-xr-x 2 sutan101 depfg   6 Apr  9 14:46 .
#~ drwxr-xr-x 3 sutan101 depfg   3 Apr  9 13:54 ..
#~ lrwxrwxrwx 1 sutan101 depfg 128 Apr  9 14:24 cdo_gridarea_clone_global_05min_correct_lats.nc -> /scratch/depfg/sutan101/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/cdo_gridarea_clone_global_05min_correct_lats.nc
#~ lrwxrwxrwx 1 sutan101 depfg 116 Apr  9 14:24 lddsound_05min_version_20210330.map -> /scratch/depfg/sutan101/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/lddsound_05min_version_20210330.map
#~ lrwxrwxrwx 1 sutan101 depfg 177 Apr  9 14:46 mask5minFromTop.map -> /scratch/depfg/hydrowld/data/hydroworld/pcrglobwb2_input_release/version_2019_11_beta_extended/pcrglobwb2_input/global_05min/cloneMaps/global_parallelization/mask5minFromTop.map
#~ lrwxrwxrwx 1 sutan101 depfg  59 Apr  9 14:24 source_from_bram -> /depfg/dropp003/pcrglobwb_reanalysis/output/develop/domains

# ldd and cell areas used in the model
ldd_map_file = "/scratch/depfg/sutan101/wmo_2025/pcrglobwb/lddsound_05min_version_20210330.map"
area_m2_file = "/scratch/depfg/sutan101/wmo_2025/pcrglobwb/cdo_gridarea_clone_global_05min_correct_lats.nc.map"

# set the pcraster clone based on the ldd
pcr.setclone(ldd_map_file)
ldd_map = pcr.readmap(ldd_map_file)
area_m2 = pcr.readmap(area_m2_file)

# mask_classes used during the parallelization
mask_parallel_run_file = "/scratch/depfg/sutan101/wmo_2025/pcrglobwb/mask5minFromTop.map"
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
    
    # calculate area_deviation, do not use it if the deviation is large
    use_this_station = True
    area_deviation = (model_area_km2_this_station - wmo_area_km2) / wmo_area_km2
    if abs(area_deviation) > 5.: use_this_station = False
    # - also do not use if if area_deviation cannot be calculated
    if np.isnan(np.array(float(area_deviation))): use_this_station = False
    
    # - check where the mask for this
    mask_for_this_station, valid = pcr.cellvalue(pcr.mapmaximum(pcr.scalar(pcr.ifthen(wmo_id_point, mask))), 1)
    if mask_for_this_station < 0.: use_this_station = False

    print(use_this_station)
    if use_this_station:

    # ~ use_all = True
    # ~ if use_all:

        # put them in the dataframe
        wmo_station_table.loc[irow, "model_lon"]      = model_lon     
        wmo_station_table.loc[irow, "model_lat"]      = model_lat     
        wmo_station_table.loc[irow, "model_area_km2"] = model_area_km2_this_station
        wmo_station_table.loc[irow, "area_deviation_in_percent"] = area_deviation * 100.
        print(model_lon)
        print(model_lat)    
        print(model_area_km2_this_station)    
        print(area_deviation)


        # - go to the selected clone
        mask_for_this_station, valid = pcr.cellvalue(pcr.mapmaximum(pcr.scalar(pcr.ifthen(wmo_id_point, mask))), 1)
        #~ mask_code   = 'M%07d' %(mask_for_this_station)
        mask_code   = 'M%02d' %(mask_for_this_station)
        
        strt_year = 1991
        last_year = 2025
        
        for year in range(strt_year, last_year+1): 
        
            # netcdf filename (example: /scratch/depfg/sutan101/wmo_2025/pcrglobwb/source_from_bram/M01/periods/20250101-20251231/netcdf/discharge_dailyTot_output.nc)
            netcdf_file = model_output_folder + "/" + mask_code + "/periods/" + str(year) + "0101-" + str(year) + "1231/" + "/netcdf/discharge_dailyTot_output.nc"     
            print(netcdf_file)
            
            # - pick the timeseris in the netcdf file using xarray
            discharge_xr          = xr.open_dataset(netcdf_file)
            discharge_time_series = discharge_xr.sel(lon = model_lon, lat = model_lat, method = 'nearest')
            
            strt_date = str(year) + "-01-01"
            last_date = str(year) + "-12-31"
            discharge_time_series = discharge_time_series.sel(time=slice(strt_date, last_date))
	        
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
            csv_filename = csv_output_folder + "/" + "pcrglobwb_" + str(wmo_id) + "_discharge_1991_2025.csv"  
            # - append if not the first year
            mode = "w" ; header = True
            if year > strt_year: mode = "a" ; header = False
            df.to_csv(csv_filename, index = False, mode = mode, header = header)
        
        print(csv_filename)    

    # write the station list to a csv file
    csv_station_list = csv_output_folder + "/" + "_pcrglobwb_wmo_station_list.csv"  
    wmo_station_table.to_csv(csv_station_list, index = False)


