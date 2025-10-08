################################################################################
# Script Name:        define_schemes.R
# Description:        Each original data set is decomposed into different schemes. 
#                     Different schemes are created for each year.
#                     Within each year only the three consecutive month with the most samples are used.
#                     Within these three month further subsets are created based on the number of samples.
#                     
#                     SAMPLING STRATEGY:
#                     Uses logarithmic spacing (geometric progression) to sample more densely 
#                     at low N where metrics stabilize rapidly, and more sparsely at high N 
#                     where changes plateau. Multiplier of ~1.5 generates sequence:
#                     50, 75, 113, 169, 254, 381, 571, 857...
#                     This captures critical transitions better than linear 100-step intervals.
#                     
#                     Additionally, if abundance data is available, different schemes are run 
#                     for abundance and for presence/absence.
# 
# Author:             Jonathan Jupke
# Date Created:       2025-09-18
# Last Modified:      2025-10-06
#
# R Version:          R 4.5.1
# Required Packages:  data.table
#
# Notes:              Modified sampling strategy from 100-step to log-scale based on
#                     red team analysis identifying non-linear metric stabilization patterns.
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
                        
                        # ———— Create log-scale sample scheme ————
                        # Rationale: Ecological metrics typically stabilize non-linearly.
                        # Log-spacing samples more densely where change is rapid (50-200 range)
                        # and more efficiently where change plateaus (>500 range).
                        # This provides better characterization of minimum viable N and 
                        # more accurate RF model fitting in critical transition zones.
                        
                        # Generate geometric sequence: 50 × 1.5^(0, 1, 2, 3...)
                        i.min_samples <- 50
                        i.multiplier <- 1.5
                        i.max_iter <- 20  # Sufficient to reach very large N
                        
                        # Generate candidate sample sizes
                        i.log_sequence <- i.min_samples * (i.multiplier ^ (0:(i.max_iter)))
                        
                        # Round to whole numbers and remove duplicates after rounding
                        i.log_sequence <- unique(round(i.log_sequence))
                        
                        # Keep only values ≤ actual sample size
                        i.sample_scheme <- i.log_sequence[i.log_sequence <= i.samples]
                        
                        # Always include the full sample size if not already present
                        if (max(i.sample_scheme) < i.samples) {
                                i.sample_scheme <- c(i.sample_scheme, i.samples)
                        }
                        
                        # If somehow we have no valid samples (shouldn't happen given i.samples >= 50)
                        if (length(i.sample_scheme) == 0) {
                                i.sample_scheme <- i.samples
                        }
                        
                        # Limit to maximum 15 schemes to prevent excessive computation
                        # (typically won't be reached unless i.samples is very large)
                        if (length(i.sample_scheme) > 15) {
                                # Keep the smallest, largest, and evenly spaced intermediate values
                                i.keep_indices <- unique(c(
                                        1,  # Always keep smallest
                                        round(seq(2, length(i.sample_scheme) - 1, length.out = 13)),
                                        length(i.sample_scheme)  # Always keep largest (full sample)
                                ))
                                i.sample_scheme <- i.sample_scheme[i.keep_indices]
                        }
                        
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
                # limit maximal number of schemes per data set to 15
                # Pick the schemes with the largest samples sizes
                if (nrow(ds.results.list) > 15){
                        setorderv(ds.results.list, "samples")
                        ds.results.list <- tail(ds.results.list, 15)
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