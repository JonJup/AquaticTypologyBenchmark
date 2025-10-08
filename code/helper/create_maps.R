################################################################################
# Script Name:        generate_maps.R
# Description:        create maps of samples
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-23
# Last Modified:      2025-09-23
#
# R Version:          R 4.5.1
# Required Packages:  package1, package2
#
# Notes:              Any notes or assumptions
################################################################################

pacman::p_load(
        sf,
        tmap,
        data.table,
        fs,
        lubridate,
        stringr,
        dplyr,
        terra,
        mapview,
        maptiles,
        magrittr
)
tmap_mode("view")

files <- c(
        "D://data/biota/diatom data/00_combine_data/PULSE/pulse_diatoms.rds",
        "D://data/biota/fish data/00_combine_data/PULSE/pulse_fish.rds",
        "D://data/biota/invertebrate data/00_combine_data/PULSE/pulse_invertebrates.rds",
        "D://data/biota/macrophyte data/00_combine_data/PULSE/pulse_macrophytes.rds"
)


files <- lapply(files, readRDS)
files2 <- lapply(files, function(x) unique(x, by = "siteID"))

files2[[1]]$taxon <- "diatoms"
files2[[2]]$taxon <- "fish"
files2[[3]]$taxon <- "invertebrates"
files2[[4]]$taxon <- "macrophytes"

sites <- rbindlist(files2, fill = T)

#- extract sites from data for each season
sites     <- st_as_sf(sites, coords = c("x.coord", "y.coord"), crs = 3035)
#- define custom color palette
custom.color.palette <- c("#EC6B4F", "#65F78D")

# mapview(sites)

tmap_mode("plot")

basemap.tile <-
        get_tiles(
                x = st_union(sites),
                provider = "OpenTopoMap",
                zoom = 1, crop = TRUE)


# tm_shape(basemap.tile) +
#         tm_rgb ()
#writeRaster(basemap.tile, "fig/basemap.tif", overwrite=TRUE)
## -- create map 
map.static <-
        tm_shape(basemap.tile) +
        tm_rgb () +
        tm_shape(sites) +
        tm_symbols(col = "#CD5C5C",
                   #col = "data.set.id",
                   shape = 21,
                   size = .1)  + 
        tm_facets(by = "taxon", free.coords = FALSE, nrow = 1) 

## -- save map to file 
tmap_save(tm = map.static, 
          filename = "output/figures/map_all.png")
## -- save map as tiff for publication 
tmap_save(tm = map.static, 
          filename = "output/figures/map_all.tiff", 
          dpi = 300)
