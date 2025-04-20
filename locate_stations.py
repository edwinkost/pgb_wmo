#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil

import pcraster as pcr
import virtualOS as vos

# ~ # table/column input file - this should be sorted based on reservoir capacities (the largest the lower ids)
# ~ column_input_file = "../aha_hydropowers/existing_reservoirs_v2024-04-10_selected_no-headers.csv"

column_input_file = "../aha_hydropowers/corrected_dams_2024-04-26_no-headers.csv"

# ~ # - for testing with 74 and 85 (multiple pixel problem)
# ~ column_input_file = "../aha_hydropowers/existing_reservoirs_v2024-04-10_selected_no-headers-test74and85.csv"


# set the clone (your map grid system and extent)
clone_map_file = "../pcrglobwb_maps/lddsound_05min_version_20210330.map"
pcr.setclone(clone_map_file)

# ldd map (river network)
ldd_map_file = "../pcrglobwb_maps/lddsound_05min_version_20210330.map"
ldd_map = pcr.readmap(ldd_map_file)

# cell area (unit: m2)
cell_area_file = "../pcrglobwb_maps/cdo_gridarea_clone_global_05min_correct_lats.nc.map"
cell_area = pcr.readmap(cell_area_file)

# hydrolakes
hydrolakes_file = "../pcrglobwb_maps/areaIDs.map"
hydrolakes_ids = pcr.readmap(hydrolakes_file)

# calculate catchment area (in km2) based on pcrglobwb ldd
catchment_area_km2 = pcr.catchmenttotal(cell_area, ldd_map) / (1000.*1000.)

# calculate rough estimation of hydrolakes surface area based on pcrglobwb cell area (unit: m2)
hydrolakes_pcrglobwb_area_m2 = pcr.areatotal(cell_area, hydrolakes_ids)

# convert table/column to a pcraster map of dam ids
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 1 -S -s ',' " + column_input_file + " dam_ids.map" 
print(cmd); os.system(cmd)
# - read dam ids as a variable
dam_ids = pcr.readmap("dam_ids.map") 

# convert table/column to a pcraster map of catchment areas based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 2 -S -s ',' " + column_input_file + " aha_catchment_area_km2.map"
print(cmd); os.system(cmd)
# - read aha_catchment_area_km2 as a variable
aha_catchment_area_km2 = pcr.readmap("aha_catchment_area_km2.map") 

# convert table/column to a pcraster map of the latitude coordinates based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 3 -S -s ',' " + column_input_file + " aha_latitudes.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_latitudes = pcr.readmap("aha_latitudes.map") 

# convert table/column to a pcraster map of the latitude coordinates based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 4 -S -s ',' " + column_input_file + " aha_longitudes.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_longitudes = pcr.readmap("aha_longitudes.map") 

# convert table/column to a pcraster map of the reservoir surface area based on AHA (unit: km2)
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 6 -S -s ',' " + column_input_file + " aha_surface_area_km2.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_surface_area_km2 = pcr.readmap("aha_surface_area_km2.map") 

# get the pcrglobwb catchment area for every dam id
dam_ids_pcrglobwb_catchment_area_km2 = pcr.ifthen(pcr.defined(dam_ids), catchment_area_km2)

# calculate the relative difference between two catchment areas
rel_dif_catchment_area = pcr.abs(aha_catchment_area_km2 - dam_ids_pcrglobwb_catchment_area_km2) / dam_ids_pcrglobwb_catchment_area_km2

# loop through all dams (from the largest to the smallest), if rel_dif_catchment_area > threshold, we have to reposition it
threshold = 0.1
number_of_dams = 131

# ~ # for testing
# ~ number_of_dams = 10
# ~ number_of_dams = 2


def get_distances_to_a_reference_latlon_coodinate_in_meter(lon_pcrmap, lat_pcrmap, lon_ref, lat_ref):
    pcr_distance_map = pcr.scalar(pcr.acos(pcr.sin(lat_pcrmap)*pcr.sin(lat_ref) + pcr.cos(lat_pcrmap)*pcr.cos(lat_ref)*pcr.cos(lon_ref - lon_pcrmap)))/360.*6371000. 
    return pcr_distance_map


for dam_id in range(1, number_of_dams + 1):
    
    print(dam_id)
    
    # make a point map of this dam
    this_dam_point = pcr.ifthen(dam_ids == dam_id, pcr.boolean(1.0))
    
    # evaluate relative difference for this particular dam
    rel_dif_catchment_area_this_dam = pcr.ifthen(this_dam_point, rel_dif_catchment_area)
    # - get its cell value
    rel_dif_catchment_area_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(rel_dif_catchment_area_this_dam),1)[0]
    
    # correcting the location of dam
    if rel_dif_catchment_area_this_dam_cell_value < threshold:

        location_corrected_dam_id = pcr.ifthen(this_dam_point, pcr.nominal(dam_id))
        
    else:
    
        # expanding the point to its neighbours
        search_window = pcr.windowmajority(this_dam_point, 5./60. * 3.)
        # - note using the window_length = 2 to avoid 'too large' window size
        
        # get the pcrglobwb catchment area within this search_window
        catchment_area_within_search_window = pcr.ifthen(pcr.defined(search_window), catchment_area_km2)
        
        # compare the above to aha catchment_area
        aha_catchment_area_km2_this_dam = pcr.ifthen(this_dam_point, aha_catchment_area_km2)
        aha_catchment_area_km2_this_dam = pcr.windowmaximum(aha_catchment_area_km2_this_dam, 5./60 * 3.)

        # calculate the absolute difference
        difference_catch_area = pcr.abs(catchment_area_within_search_window - aha_catchment_area_km2_this_dam)
        
        location_corrected_dam_id = pcr.ifthen(difference_catch_area == pcr.mapminimum(difference_catch_area), pcr.nominal(dam_id))
        

        # make sure that there is only one pixel in location_corrected_dam_id; if this is NOT the case, we just find the closest one to the original coordinates 

        # - get the latitude and longitude coordinates based on AHA, and transfer them to single values maps (using mapmaximum)
        lat_aha_coordinate = pcr.mapmaximum(pcr.ifthen(dam_ids == dam_id, aha_latitudes))
        lon_aha_coordinate = pcr.mapmaximum(pcr.ifthen(dam_ids == dam_id, aha_longitudes))
        
        # - calculate the distance of every candidate point to the AHA coordinate
        distance_map = get_distances_to_a_reference_latlon_coodinate_in_meter(lon_pcrmap = pcr.xcoordinate(pcr.defined(location_corrected_dam_id)),\
                                                                              lat_pcrmap = pcr.ycoordinate(pcr.defined(location_corrected_dam_id)),\
                                                                              lon_ref = lon_aha_coordinate,\
                                                                              lat_ref = lat_aha_coordinate)
        area_order = pcr.areaorder(distance_map, location_corrected_dam_id)
        # - choose the one with shortest distance (order/rank = 1)
        location_corrected_dam_id = pcr.ifthen(area_order == 1, pcr.nominal(dam_id))
        
    if dam_id == 1:    
        all_location_corrected_dam_ids = location_corrected_dam_id
    else:
        all_location_corrected_dam_ids = pcr.cover(all_location_corrected_dam_ids, location_corrected_dam_id)
    
    # get the reservoir surface area based on AHA (unit: m2)
    this_dam_point = pcr.defined(location_corrected_dam_id)
    aha_surface_area_m2 = pcr.ifthen(this_dam_point, aha_surface_area_km2) * 1000.*1000.
    aha_surface_area_m2_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(aha_surface_area_m2),1)[0]
    
    # get the cell area for this dam (unit: m2)
    cell_area_m2_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(pcr.ifthen(this_dam_point, cell_area)))
    
    # obtaining the reservoir extent (including estimating surface area)
    if aha_surface_area_m2_this_dam_cell_value < cell_area_m2_this_dam_cell_value:
        
        # this means the reservoir extent will be within a cell
        reservoir_surface_area           = pcr.ifthen(this_dam_point, aha_surface_area_m2_this_dam_cell_value)
        reservoir_surface_area_per_cell  = reservoir_surface_area
        
        # fraction of surface water within a cell
        reservoir_fraction_water = reservoir_surface_area_per_cell/cell_area

    else:

        # this is for the case if reservoirs covering multiple cells
        
        # make a search window to find the nearest hydrolakes - expanding the point to its neighbours
        search_window = pcr.windowmajority(this_dam_point, 5./60. * 3.)
        # - note using the window_length = 2 to avoid 'too large' window size
        
        # find the hydrolakes ids within this search window 
        hydrolakes_ids_within_search_window = pcr.ifthen(pcr.defined(search_window), hydrolakes_ids)
        
        # check whether there are more than one hydrolakes_ids_within_search_window
        number_of_hydrolakes_ids_within_search_window = pcr.cellvalue(pcr.maximum(pcr.scalar(pcr.clump(hydrolakes_ids_within_search_window))))[0]

        if number_of_hydrolakes_ids_within_search_window > 0:
            
            # find the one that has the most similar surface area to the estimate on hydrolakes_pcrglobwb_area
            absolute_difference_surface_area = pcr.abs(hydrolakes_pcrglobwb_area_m2 - aha_surface_area_m2_this_dam_cell_value)
            absolute_difference_surface_area = pcr.ifthen(pcr.defined(search_window), absolute_difference_surface_area)
            
            class_map_boolean = pcr.defined(absolute_difference_surface_area)
            class_map_nominal = pcr.ifthen(class_map_boolean, pcr.nominal(1.0))
            
            area_order = pcr.areaorder(absolute_difference_surface_area, class_map_nominal)
            
            # - choose the one with order/rank = 1 (minimum difference in surface area)
            hydrolakes_id_selected = pcr.mapmaximum(pcr.scalar(pcr.ifthen(area_order == 1, pcr.nominal(hydrolakes_ids_within_search_window))))
            
            # - the chosen hydro lakes id and assign in to the dam point
            hydrolakes_id_for_this_dam_point = pcr.ifthen(this_dam_point, pcr.nominal(hydrolakes_id_selected))
            
            # make the sub catchment until the dam point and its ldd
            sub_catchment = pcr.subcatchment(ldd_map, this_dam_point)
            ldd_above_the_dam_point = pcr.lddmask(ldd_map, sub_catchment)
            
            # - expand it until the entire reservoirs
            hydrolakes_id_for_this_dam_point = pcr.path(ldd_above_the_dam_point, hydrolakes_ids)
            hydrolakes_id_for_this_dam_point = pcr.ifthen(hydrolakes_id_for_this_dam_point == pcr.nominal(hydrolakes_id_selected), hydrolakes_id_for_this_dam_point)
            
            # identify the number of cells based on the hydrolakes
            number_of_cells_according_to_hydrolakes = pcr.areatotal(pcr.scalar(1.0), hydrolakes_id_for_this_dam_point)
		    
            # obtaining reservoir surface area
            reservoir_surface_area_per_cell = aha_surface_area_m2_this_dam_cell_value / number_of_cells_according_to_hydrolakes
            
            # fraction of surface water within a cell
            reservoir_fraction_water = reservoir_surface_area_per_cell / cell_area

        else:

            # if not identified in the hydrolakes

            # start with assuming everything from zero
            reservoir_surface_area             = 0.0
            reservoir_surface_area_per_cell    = 0.0
            reservoir_fraction_water           = 0.0
            number_of_cells_for_this_reservoir = 0.0

            # assign the current cell as the reservoir 
            reservoir_extent                   = this_dam_point
            reservoir_surface_area             = pcr.ifthen(reservoir_extent, cell_area)
            reservoir_surface_area_per_cell    = pcr.ifthen(reservoir_extent, cell_area)
            reservoir_fraction_water           = pcr.ifthen(reservoir_extent, 1.0)
            number_of_cells_for_this_reservoir = 1.0
            
            # calculate the remaining surface area that needs to be covered
            remaining_area = aha_surface_area_m2_this_dam_cell_value - pcr.cellvalue(pcr.maptotal(reservoir_surface_area_per_cell))
            
            while remaining_area > 0:
            
                # identify the upstream cells
                upstream_cells_of_reservoirs = pcr.boolean(pcr.downstream(ldd_map, pcr.cover(pcr.scalar(reservoir_cells), 0.0)))

                # number of upstream cells
                num_of_upstream_cells = pcr.maptotal(pcr.scalar(upstream_cells_of_reservoirs))
                
                # reservoir surface area at the upstream cells
                additional_reservoir_surface_area_per_cell = pcr.ifthen(upstream_cells_of_reservoirs, pcr.min(remaining_area / num_of_upstream_cells, cell_area))
                
                # the updated reservoir surface area per cell
                reservoir_surface_area_per_cell = pcr.cover(reservoir_surface_area_per_cell, additional_reservoir_surface_area_per_cell)
                
                # calculate the remaining surface area that needs to be covered
                remaining_area = aha_surface_area_m2_this_dam_cell_value - pcr.cellvalue(pcr.maptotal(reservoir_surface_area_per_cell))
            
            # update the reservoir with upstream cells
            reservoir_extent         = pcr.cover(reservoir_extent, pcr.ifthen(reservoir_surface_area_per_cell > 0.0, pcr.boolean(1.0))
            reservoir_surface_area   = pcr.ifthen(reservoir_extent, pcr.maptotal(reservoir_surface_area_per_cell))
            reservoir_fraction_water = reservoir_surface_area_per_cell / cell_area
        

    # we summarize all reservoirs into single variables 
    if dam_id == 1:    
        all_reservoir_extent_ids     = pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(all_location_corrected_dam_ids)))
        all_reservoir_surface_area   = pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(reservoir_surface_area)))
        all_reservoir_fraction_water = pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(all_reservoir_fraction_water))) 
    else:
        all_reservoir_extent_ids     = pcr.cover(all_reservoir_extent_ids,     pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(all_location_corrected_dam_ids))))
        all_reservoir_surface_area   = pcr.cover(all_reservoir_surface_area,   pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(reservoir_surface_area))))
        all_reservoir_fraction_water = pcr.cover(all_reservoir_fraction_water, pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(all_reservoir_fraction_water)))) 


# save all_location_corrected_dam_ids to a pcraster map - this will be the location of all dams/outlets
pcr.report(all_location_corrected_dam_ids, "corrected_dam_ids.map")

# obtain the catchment areas of all_location_corrected_dam_ids
corrected_dam_catchment_area_km2 = pcr.ifthen(pcr.defined(all_location_corrected_dam_ids), catchment_area_km2)
pcr.report(corrected_dam_catchment_area_km2, "corrected_dam_catchment_area_km2.map")

# obtain a table/column format for all_location_corrected_dam_ids and corrected_dam_catchment_area_km2
cmd = "map2col corrected_dam_ids.map corrected_dam_catchment_area_km2.map corrected_dams.txt"
print(cmd); os.system(cmd)        

# save also the following variables for pcrglobwb input:
pcr.report(all_reservoir_extent_ids,     "reservoir_extent_ids.map")
pcr.report(all_reservoir_surface_area,   "reservoir_surface_area_ids.map")
pcr.report(all_reservoir_fraction_water, "reservoir_fraction_water_ids.map")
