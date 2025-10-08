################################################################################
# Script Name:        generate_summary.R
# Description:        This script generates summary statistics for the data sets. 
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-23
# Last Modified:      2025-09-23
#
# R Version:          R 4.5.1
# Required Packages:  package1, package2
#
# Notes:              
################################################################################

setwd(rstudioapi::getActiveProject())

library(data.table)

files <- c(
        "D://data/biota/diatom data/00_combine_data/PULSE/pulse_diatoms.rds",
        "D://data/biota/fish data/00_combine_data/PULSE/pulse_fish.rds",
        "D://data/biota/invertebrate data/00_combine_data/PULSE/pulse_invertebrates.rds",
        "D://data/biota/macrophyte data/00_combine_data/PULSE/pulse_macrophytes.rds"
)



files <- lapply(files, readRDS)

# one data set has a placeholder data. this is removed 
files[[4]] <- files[[4]][data.set != "macrophytes_sweden_leo"]

files <- lapply(files, function(x) x[, MD := min(eventDate), by = "data.set"])
files <- lapply(files, function(x) x[, MD2 := max(eventDate), by = "data.set"])
files <- lapply(files, function(x) x[, date.range := MD2-MD])
date_range <- lapply(files, function(x) x[, MD2-MD])
date_range <- sapply(date_range, unique)
date_range <- unlist(date_range)
shortest_time <- min(date_range, na.rm = T)
longest_time <- max(date_range, na.rm = T)
shortest_time <- paste(shortest_time, "days")
longest_time <- paste(floor(longest_time/365), "years")
average_time <- floor(median(date_range, na.rm = T))
average_time <- paste(average_time, "days")

files <- lapply(files, function(x) x[, DSS  := uniqueN(eventID), by = "data.set"])
uf <- lapply(files, unique, by = "data.set")
uf <- rbindlist(uf, fill =T )
uf <- summary(uf$DSS)
out <- list(
        shortest_time = shortest_time,
        longest_time  = longest_time ,
        average_time  = average_time, 
        small_dataset = paste(uf[1], "samples"),
        large_dataset = paste(uf[6], "samples"),
        median_dataset = paste(uf[3], "samples")
)
saveRDS(out, "output/summary_statistics.rds")
