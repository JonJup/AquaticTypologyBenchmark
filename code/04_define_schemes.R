################################################################################
# Script Name:        define_schemes.R
# Description:        Each original data set is decomposed into different schemes. 
#                     Different schemes are created for each year.
#                     Within each year only the three consecutive month with the most samples are used.
#                     Within these three month further subsets are created based on the number of samples.
#                     One scheme is run with the complete data from the focal month but 
#                     further schemes are created in decreasing steps of 100. For example, 
#                     if the three focal month contain 563 distinct samples, one scheme uses 563 samples and further schemes use
#                     500, 400, 300, 200, and 100 randomly selected samples. 
#                     Additionally, if abundance data is available, different schemes are run for abundance and for presence/absence.
#                     To reuse the previous example, the currently six schemes would double to 12 schemes if abundance data were available.
# 
# Author:             Jonathan Jupke
# Date Created:       2025-09-18
# Last Modified:      2025-10-01
#
# R Version:          R 4.5.1
# Required Packages:  package1, package2
#
# Notes:              
################################################################################

# setup -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())
library(data.table)
source("code/functions/find_max_consequtive_sum.R")
# load and prepare data -------------------------------------------------------------------------
files      <- list.files("data/biota", full.names = T, pattern = "02_")
bio.names <- sapply(files,     function (x) sub("data/biota/02_", "", x)) 
bio.names <- sapply(bio.names, function (x) sub("_w_environment.rds", "", x)) 

bio.list   <- lapply(files, readRDS)
n_bio_datasets <- length(bio.list)

id.to.enz <- readRDS("data/eu_hydro_dem_w_enz.rds")

# set loop independent parameters  --------------------------------------------------
n.data.sets      <- lapply(bio.list, function(x) uniqueN(x$data.set))
ud               <- lapply(bio.list, function(x) sort(unique(x$data.set)))

# 
# setnames(bio.list[[2]], old = "year", new = "eventYear")
# setnames(bio.list[[2]], old = "lowest.taxon", new = "working.taxon")
# setnames(bio.list[[4]], old = "lowest.taxon", new = "working.taxon")

# prepare data ----------------------------------------------------------------------
result_list <- vector(mode = "list")
for (b in 1:n_bio_datasets){
        for (d in 1:n.data.sets[[b]]) {
                print(d)
                ds.result.list <- list()
                #- select data set
                ds.data.set <- ud[[b]][d]
                ds.bio      <- bio.list[[b]][data.set == ds.data.set]
                
                if (ds.data.set == "diatoms_germany_bavaria_federal_monitoring"){
                        ds.bio <- ds.bio[organismQuantityType == "individuals per m2"]
                }
                # TODO quick fix for Brandenburg diatoms remove eventually 
                ds.check.vec <- c("Abundanzen (alle) ohne cf", "Zähldaten Taxa ohne cf", "Abundanzen (alle)", "Abundanzen (alle) mit cf")
                if (any(ds.bio$organismQuantityType %in% ds.check.vec)){
                        ds.bio[organismQuantityType %in% ds.check.vec, organismQuantityType := "indiviuals"]
                }
                #- how many years?
                if (class(ds.bio$eventYear) == "list") {
                        ds.bio[, eventYear := unlist(eventYear)]
                }
                ds.n.year <- uniqueN(ds.bio$eventYear)
                #- how many samples per year?
                ds.n.samples <- ds.bio[, uniqueN(eventID), by = "eventYear"]
                
                if (all(ds.n.samples$V1 < 50)) {
                        print(paste("data set", ds.data.set, "contains less than 50 samples in each year"))
                        rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
                        next()
                }
                #- select years with more than 50 samples
                ds.n.samples <- ds.n.samples[V1 >= 50]
                #- remove NA row, if one of the entries in the year table is NA
                if (!all(is.na(ds.n.samples$eventYear)) & any(is.na(ds.n.samples$eventYear))){
                        ds.n.samples <- ds.n.samples[-which(is.na(ds.n.samples$eventYear))]
                }
                
                # START LOOP i OVER years in data set
                for (i in 1:nrow(ds.n.samples)) {
                        
                        if (all(is.na(ds.bio$eventYear))) {
                                i.data    <- copy(ds.bio)
                        } else if (all(is.na(ds.bio$eventDate))){
                                i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]
                        } else {
                                #- Create subset of focal year.
                                i.data <- ds.bio[eventYear == ds.n.samples$eventYear[i]]
                                #--- Check seasons. To prevent strong seasonal changes from influencing
                                #--- the community composition, we identify the three consecutive month
                                #--- with the most samples.
                                #- create month variable
                                i.data[, month := month(eventDate)]
                                #- count samples per month
                                i.month_table <- unique(i.data, by = "eventID")
                                i.month_table <- i.month_table$month
                                i.month_table <- table(i.month_table)
                                #- Which three consecutive months have the most samples?
                                #- This function is loaded in the beginning of the script
                                if (length(i.month_table) > 2) {
                                        i.max_month    <- find_max_consecutive_sum(i.month_table)
                                        i.focal.months <- names(i.max_month$values)
                                        
                                        #- create subset of i.data only containing the focal months
                                        i.data <- i.data[month %in% i.focal.months]
                                        
                                } else {
                                        #- Are the two months consecutive?
                                        i.monthdiff <- diff(as.numeric(names(
                                                i.month_table
                                        )))
                                        #- No or only one month
                                        if (length(i.monthdiff) == 0){
                                                # Does one of the month have more than 100 samples?
                                                if (any(i.month_table > 100)) {
                                                        i.focal.months <- names(which.max(
                                                                i.month_table
                                                        ))
                                                        i.data <- i.data[month == i.focal.months]
                                                } else {
                                                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                                                        next()
                                                }
                                        } else if (i.monthdiff > 2) {
                                                # Does one of the month have more than 100 samples?
                                                if (any(i.month_table > 100)) {
                                                        i.focal.months <- names(which.max(
                                                                i.month_table
                                                        ))
                                                        i.data <- i.data[month == i.focal.months]
                                                } else {
                                                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                                                        next()
                                                }
                                        } else {
                                                #- yes
                                                i.focal.months <- names(i.month_table)
                                                i.data <- i.data[month %in% i.focal.months]
                                        }
                                }
                        }
                        
                        #- how many samples are left?
                        i.samples <- uniqueN(i.data$eventID)
                        if (i.samples < 50) {
                                rm(list = ls()[grepl("^i\\.", x = ls())])
                                next()
                        }
                        
                        # ———— What sorts of random combinations can be created from that?
                        i.hundreds    <- floor(i.samples / 100)
                        #if (i.hundreds > 10) i.hundreds <- 10
                        # In cases where we have less than 100 samples, we do not create further subsamples. 
                        if (i.hundreds == 0){
                                i.sample_scheme <- i.samples
                        } else {
                                i.sample_scheme <- c()
                                for (crs in 1:i.hundreds) {
                                        i.sample_scheme <-  
                                                append(i.sample_scheme, 
                                                       rep(crs * 100, times = 1)
                                                       )
                                }
                                rm(crs)
                                i.sample_scheme <- append(i.sample_scheme, i.samples)
                        }

                        # if ((i.samples / 100 - i.hundreds) > 0.5) {
                        #         i.sample_scheme <- append(i.sample_scheme, max(i.sample_scheme))
                        # }
                        rm(i.hundreds)
                        # check what environmental zones are included
                        i.enz <- id.to.enz[ID %in% i.data$ID]
                        i.enz <- unique(i.enz$EnZ_name)
                        i.enz <- sort(i.enz)
                        
                        #i.sample_scheme
                        i.out <- data.table(
                                taxon    = bio.names[b],
                                data.set = ud[[b]][d],
                                eventYear     = unique(i.data$eventYear),
                                catchments = list(i.enz),
                                samples  = i.sample_scheme,
                                sample_type = c(rep("sub", length(i.sample_scheme)-1), "full"),
                                n_taxa = uniqueN(i.data$working.taxon)
                        )
                        if (all(is.na(ds.bio$eventYear))| all(is.na(ds.bio$eventDate))) {
                                i.out[, focal_months := list(list(NA))]
                        } else {
                                i.out[, focal_months := list(list(i.focal.months))]
                        }
                        # check abundance type
                        if (all(i.data$organismQuantityType == "individuals")){
                                
                                i.out2 <- copy(i.out)
                                i.out3 <- copy(i.out)
                                
                                i.out2[, organismQuantityType := "presence"]
                                i.out3[, organismQuantityType := "individuals"]
                                i.out <- rbindlist(list(i.out2, i.out3))
                                      
                        } else {
                                i.out[, organismQuantityType := "presence"]
                        }
                        
                        ds.result.list[[length(ds.result.list) + 1]] <- i.out
                        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
                } # END of loop i over years in data.set d
                
                # combine results of all years. 
                ds.results.list <- rbindlist(ds.result.list)
                # limit maximal number of schemes per data set to 10.
                # Pick the 10 schemes with the largest samples size 
                if (nrow(ds.results.list) > 10){
                        setorderv(ds.results.list, "samples")
                        ds.results.list <- tail(ds.results.list, 10)
                }
                
                result_list[[length(result_list) + 1]] <- ds.results.list
                
                rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
                
        } # END d of loop over data sets in taxonomic group b
} # END of loop b over taxonomic groups

result.data <- rbindlist(result_list, fill = T)
result.data$taxon <- factor(result.data$taxon) 
result.data <- split(result.data, f = result.data$taxon)
result.data <- lapply(result.data, function(x) x[, number := .GRP, by =  c("data.set", "eventYear", "samples", "organismQuantityType", "sample_type")])
result.data <- rbindlist(result.data)

#result.data[, number := .GRP, by = c("taxon", "data.set", "eventYear", "samples", "organismQuantityType", "sample_type")]
result.data[, number := as.character(number)]
result.data[as.numeric(number) < 10, number := paste0("000", number)]
result.data[as.numeric(number) >= 10 & as.numeric(number) < 100, number := paste0("00", number)]
result.data[as.numeric(number) >= 100, number := paste0("0", number)]
result.data[, scheme_id := paste0(taxon, "_",number)]
result.data[, number := NULL]
MCMCParameterMatrix <- data.table(
        
        scheme_id = result.data$scheme_id,
        nChains   = 3,
        nSamples  = 2000, 
        thin      = 2,
        transient = 4000,
        converged = FALSE,
        smallAUC  = FALSE,
        envelope  = FALSE,
        c_score   = FALSE,
        spec_env  = FALSE,
        occ_dist  = FALSE,
        taxa_dropped = character(1),
        variables_dropped = character(1)
)
split.data <- split(result.data, by = "taxon")


# data output -------------------------------------------------------------
lapply(1:n_bio_datasets,
       function(x)
               saveRDS(
                       split.data[[x]], 
                       paste0("data/biota/03_",bio.names[x],"_scheme.rds")
               )
)
fwrite(MCMCParameterMatrix, "data/mcmc_parameters_null.csv")
fwrite(MCMCParameterMatrix, "data/mcmc_parameters.csv")
