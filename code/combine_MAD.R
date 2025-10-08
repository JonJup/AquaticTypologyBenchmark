################################################################################
# Script Name:        combine_MAD.R
# Description:        combine all the MAD values from different environmental variables. 
#                     Evaluate which ones are high enough to merit inclusion        
# Author:             Your Name
# Date Created:       2025-09-11
# Last Modified:      2025-09-11
#
# R Version:          R 4.5.1
# Required Packages:  
#
# Notes:              
################################################################################


# setuo -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())

library(data.table)

# var                   <- readRDS("D://Arbeit/data/river_network/hydrography90m/processedData/slope_hydroDem_mad.rds")
# var                 <- readRDS("D://Arbeit/data/soil/new soil/processedData/soilOC_MAD.rds")
# var                 <- readRDS("D://Arbeit/data/soil/soil_ph_europe/processedData/mad_soil_ph_hydro_dem.rds")
# var                 <- readRDS("D://Arbeit/data/river_network/flood_area/processedData/inundationDepth_euHydro_mad.rds")
# var                 <- readRDS("D://Arbeit/data/river_network/hydrography90m/processedData/spi_hydroDem_mad.rds")
# var                  <- readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim05_euHydro_mad.rds")
# var                  <- readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim06_euHydro_mad.rds")
# var                  <- readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim13_euHydro_mad.rds")
# var                  <- readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim14_euHydro_mad.rds")
# var                 <- readRDS("D://Arbeit/data/misc/Gloreda/processedData/gloreda_max_euHydro_mad.rds")
# var                 <- readRDS("D://Arbeit/data/misc/Gloreda/processedData/gloreda_min_euHydro_mad.rds")
# var                 <- readRDS("D://Arbeit/data/misc/Gloreda/processedData/gloreda_mean_euHydro_mad.rds")
# var               <- readRDS("D://Arbeit/data/DEM/DTM_Europe/processedData/roughness_hydroDEM_mad.rds")
# var               <- readRDS("D://Arbeit/data/DEM/DTM_Europe/processedData/elevation_hydroDEM_mad.rds")
# var               <- readRDS("D://Arbeit/data/CERRA/processedData/snow_depth_hydoDem_mad.rds")
# var               <- readRDS("D://Arbeit/data/river_network/discharge copernicus/processedData/discharge_max_hydroDEM_mad.rds")
# var               <- readRDS("D://Arbeit/data/river_network/discharge copernicus/processedData/discharge_mean_hydroDEM_mad.rds")
# var               <- readRDS("D://Arbeit/data/river_network/discharge copernicus/processedData/discharge_min_hydroDEM_mad.rds")
# var             <- readRDS("D://Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/groundwater_hydroDEM_mad.rds")
# var           <- readRDS("D://Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/upstream_hydroDEM_mad.rds")
# var <- readRDS("D://Arbeit/data/river_network/soilHydraulicPropertiesEurope/processedData/SHP_hydroDEM_mad.rds")
# var  <- readRDS("D://Arbeit/data/river_network/valley_bottom_flattness/processedData/VBM_hydroDEM_mad.rds")

# Create list with all file reads
data_list <- list(
        slope = readRDS("D://Arbeit/data/river_network/hydrography90m/processedData/slope_hydroDem_mad.rds"),
        soilOC = readRDS("D://Arbeit/data/soil/new soil/processedData/soilOC_MAD.rds"),
        soil_ph = readRDS("D://Arbeit/data/soil/soil_ph_europe/processedData/mad_soil_ph_hydro_dem.rds"),
        inundationDepth = readRDS("D://Arbeit/data/river_network/flood_area/processedData/inundationDepth_euHydro_mad.rds"),
        spi = readRDS("D://Arbeit/data/river_network/hydrography90m/processedData/spi_hydroDem_mad.rds"),
        chelsaBioclim05 = readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim05_euHydro_mad.rds"),
        chelsaBioclim06 = readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim06_euHydro_mad.rds"),
        chelsaBioclim13 = readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim13_euHydro_mad.rds"),
        chelsaBioclim14 = readRDS("D://Arbeit/data/climate/bioclim_chelsa/processedData/chelsaBioclim14_euHydro_mad.rds"),
        gloreda_max = readRDS("D://Arbeit/data/misc/Gloreda/processedData/gloreda_max_euHydro_mad.rds"),
        gloreda_min = readRDS("D://Arbeit/data/misc/Gloreda/processedData/gloreda_min_euHydro_mad.rds"),
        gloreda_mean = readRDS("D://Arbeit/data/misc/Gloreda/processedData/gloreda_mean_euHydro_mad.rds"),
        roughness = readRDS("D://Arbeit/data/DEM/DTM_Europe/processedData/roughness_hydroDEM_mad.rds"),
        elevation = readRDS("D://Arbeit/data/DEM/DTM_Europe/processedData/elevation_hydroDEM_mad.rds"),
        snow_depth = readRDS("D://Arbeit/data/CERRA/processedData/snow_depth_hydoDem_mad.rds"),
        discharge_max = readRDS("D://Arbeit/data/river_network/discharge copernicus/processedData/discharge_max_hydroDEM_mad.rds"),
        discharge_mean = readRDS("D://Arbeit/data/river_network/discharge copernicus/processedData/discharge_mean_hydroDEM_mad.rds"),
        discharge_min = readRDS("D://Arbeit/data/river_network/discharge copernicus/processedData/discharge_min_hydroDEM_mad.rds"),
        groundwater = readRDS("D://Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/groundwater_hydroDEM_mad.rds"),
        upstream = readRDS("D://Arbeit/data/river_network/BasinATLAS_Data_v10/processedData/upstream_hydroDEM_mad.rds"),
        SHP = readRDS("D://Arbeit/data/river_network/soilHydraulicPropertiesEurope/processedData/SHP_hydroDEM_mad.rds"),
        VBM = readRDS("D://Arbeit/data/river_network/valley_bottom_flattness/processedData/VBM_hydroDEM_mad.rds")
)

name_vec <- c()
for (i in 1:length(data_list)) name_vec <- append(name_vec, rep(names(data_list)[i], length(data_list[[i]])))


dt <- data.table(
        variable = name_vec,
        value = unlist(data_list)
)
dt2 <- dt[, median(value, na.rm = T), by = "variable"]
setorderv(dt2, "V1")
dt2
