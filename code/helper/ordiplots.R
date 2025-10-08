################################################################################
# Script Name:        Ordination uncertainy.R
# Description:        Short description of the script
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-26
# Last Modified:      2025-09-26
#
# R Version:          R x.x.x
# Required Packages:  package1, package2
#
# Notes:              Any notes or assumptions
################################################################################

setwd(rstudioapi::getActiveProject())

library(gllvm)
library(data.table)
library(lubridate)

schemes <- list.files("data/biota", pattern = "03", full.names=T)

for (i in 1:4){
        
        biota <- readRDS(list.files("data/biota", pattern = "02", full.names = T)[i])        
        scheme <- readRDS(schemes[i])
        
        for (ii in 1:nrow(scheme)){
                
                schem <- scheme[ii, ]
                bio <- biota[data.set == schem$data.set & eventYear == schem$eventYear]
                bio <- bio[month(eventDate) %in% as.numeric(unlist(schem$focal_months))]        
                bio2 <- dcast(bio, formula = eventID ~ working.taxon , value.var = "PA", fill = 0)
                bio2[, eventID := NULL]
                test <- gllvm(y = bio2, family = "binomial")
                title <- paste(schem$data.set, "-", schem$eventYear, "-", schem$samples)
                png(paste0("output/figures/ordiplots/",title, ".png"))
                ordiplot(test, predict.region = T, main = title)
                dev.off()
        }
        
        
}
