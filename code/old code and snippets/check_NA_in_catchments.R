library(fs)
library(arrow)
library(data.table)
library(sf)
library(mapview)
library(tidyverse)
library(magrittr)
files <- dir_ls("data/catchments/")
files2 <- lapply(files, read_parquet)
files3 <- rbindlist(files2, fill = TRUE)
#dropv <- readRDS("data/dropv_catchments.rds")
files3 <- files3[!ID %in% dropv]
look_for_na <- sort(colSums(is.na(files3)))

miss <- files3[is.na(flooded_area), ID]
geo <- readRDS("E://Arbeit/Data/river_network/eu_hydro_dem/all_eu_hydro_dem_catchments.rds")

any(miss %in% geo$ID)
geo2 <- filter(geo, ID %in% miss) 
mapview(geo2)
geo2 <- filter(geo, str_detect(ID, "Mesima"))
files3[str_detect(ID, "Mesima")]


dropv <- miss_soil
saveRDS(dropv, "001_Uni/001_projects/PULSE/WP1/AquaticTypologyBenchmark/data/dropv_catchments.rds")

files4 <- copy(files3)
files4 <- files4[!ID %in% miss_soil]
# files4[is.na(upstream_catchment_area) & !is.na(upstream_catchment_area.x), upstream_catchment_area := upstream_catchment_area.x]
# files4[is.na(upstream_catchment_area) & !is.na(upstream_catchment_area.y), upstream_catchment_area := upstream_catchment_area.y]
# files4[, upstream_catchment_area.x := NULL]
# files4[, upstream_catchment_area.y := NULL]
(look_for_na <- sort(colSums(is.na(files4))))

miss_flood <- files4[is.na(Rfactor_max), ID]
geo2 <- filter(geo, ID %in% miss_flood) 
mapview(geo2)
dropv <- readRDS("data/dropv_catchments.rds")
dropv <- append(dropv, miss_flood)
saveRDS(dropv, "data/dropv_catchments.rds")
