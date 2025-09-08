# ——————————————————————————————— #
# ——— Define Sampling Schemes ——— # 
# ——————————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 29.01.2025

# load and prepare data -------------------------------------------------------------------------

library(data.table)
source("code/functions/find_max_consequtive_sum.R")
files      <- list.files("data/biota", full.names = T, pattern = "02_")
bio.names <- sapply(files, function (x) sub("data/biota/02_", "", x)) 
bio.names <- sapply(bio.names, function (x) sub("_w_environment.rds", "", x)) 
bio.list   <- lapply(files, readRDS)
n_bio_datasets <- length(bio.list)
# out_list <- vector("list", n_bio_datasets)
# names(out_list) <- paste0("out", 1:n_bio_datasets)

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
                #- how many years?
                if (class(ds.bio$year) == "list") {
                        ds.bio[, year := unlist(year)]
                }
                ds.n.year <- uniqueN(ds.bio$year)
                #- how many samples per year?
                ds.n.samples <- ds.bio[, uniqueN(comb_sample_id), by = "year"]
                
                if (all(ds.n.samples$V1 < 100)) {
                        print(paste("data set", ds.data.set, "contains less than 100 samples in each year"))
                        rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
                        next()
                }
                #- select years with more than 100 samples
                ds.n.samples <- ds.n.samples[V1 >= 100]
                #- remove NA row, if one of the entries in the year table is NA
                if (!all(is.na(ds.n.samples)) & any(is.na(ds.n.samples$year))){
                        ds.n.samples <- ds.n.samples[-which(is.na(ds.n.samples$year))]
                }
                
                # START LOOP i OVER years in data set
                for (i in 1:nrow(ds.n.samples)) {
                        
                        if (all(is.na(ds.bio$year))) {
                                i.data    <- copy(ds.bio)
                        } else if (all(is.na(ds.bio$date))){
                                i.data <- ds.bio[year == ds.n.samples$year[i]]
                        } else {
                                #- Create subset of focal year.
                                i.data <- ds.bio[year == ds.n.samples$year[i]]
                                #--- Check seasons. To prevent strong seasonal changes from influencing
                                #--- the community composition, we identify the three consecutive month
                                #--- with the most samples.
                                #- create month variable
                                i.data[, month := month(date)]
                                #- count samples per month
                                i.month_table <- unique(i.data, by = "comb_site_id")
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
                        i.samples <- uniqueN(i.data$comb_site_id)
                        if (i.samples < 100) {
                                rm(list = ls()[grepl("^i\\.", x = ls())])
                                next()
                        }
                        
                        # ———— What sorts of random combinations can be created from that?
                        i.hundreds    <- floor(i.samples / 100)
                        if (i.hundreds > 10) i.hundreds <-10
                        i.sample_scheme <- c()
                        for (crs in 1:i.hundreds) {
                                i.sample_scheme <-  append(i.sample_scheme, rep(crs * 100, times = min(
                                        i.hundreds - (crs - 1), 3
                                )))
                        }
                        rm(crs)
                        if ((i.samples / 100 - i.hundreds) > 0.5) {
                                i.sample_scheme <- append(i.sample_scheme, max(i.sample_scheme))
                        }
                        rm(i.hundreds)
                        #i.sample_scheme
                        i.out <- data.table(
                                taxon    = bio.names[b],
                                data.set = ud[[b]][d],
                                year     = unique(i.data$year),
                                samples  = i.sample_scheme
                        )
                        if (all(is.na(ds.bio$year))| all(is.na(ds.bio$date))) {
                                i.out[, focal_months := list(list(NA))]
                        } else {
                                i.out[, focal_months := list(list(i.focal.months))]
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

result.data <- rbindlist(result_list)
split.data <- split(result.data, by = "taxon")
lapply(1:n_bio_datasets,
       function(x)
               saveRDS(
                       split.data[[x]], 
                       paste0("data/biota/03_",bio.names[x],"_scheme.rds")
               )
)

