# ——————————————————————————————— #
# ——— Define Sampling Schemes ——— # 
# ——————————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 29.01.2025

# load data -------------------------------------------------------------------------
bio <- read_parquet("data/fish_w_env.parquet")

# set loop independent parameters  --------------------------------------------------
n.data.sets      <- uniqueN(bio$data.set)
ud               <- sort(unique(bio$data.set))
# prepare data ----------------------------------------------------------------------
result_list <- vector(mode = "list")
for (d in 1:n.data.sets){
#for (d in 1:3){
        
        #if (data.set == 1) result_list <- vector(mode = "list")
        
        #- select data set 
        ds.data.set <- ud[d]
        ds.bio      <- bio[data.set == ds.data.set]
        #- how many years? 
        ds.n.year <- uniqueN(ds.bio$year)
        #- how many samples per year? 
        ds.n.samples <- ds.bio[,uniqueN(comb_sample_id), by = "year"]
        
        if (all(ds.n.samples$V1 < 100)){
                print(paste("data set", ds.data.set, "contains less than 100 samples in each year")) 
                rm(list = ls()[grepl(pattern="^ds\\.", x = ls())])
                next()
        } 
        #- select years with more than 100 samples 
        ds.n.samples <- ds.n.samples[V1 >= 100]
        for (i in 1:nrow(ds.n.samples)){ # START LOOP i OVER years in data set
                #- updater 
                #print(paste("data set: ", data.set, "; year: ", i, "of", nrow(ds.n.samples)))
                
                #- Create subset of focal year. 
                i.data <- ds.bio[year == ds.n.samples$year[i]]
                
                #--- Check seasons. To prevent strong seasonal changes from influencing  
                #--- the community composition, we identify the three consecutive month  
                #--- with the most samples.
                
                #- create month variable
                i.data[, month := month(date)]
                #- count samples per month 
                i.month_table <- 
                        unique(i.data, by = "comb_site_id") %>% 
                        pull(month) %>% 
                        table
                
                
                #- Which three consecutive months have the most samples? 
                #- This function is defined in 00_setup.R
                if (length(i.month_table) > 2){
                        
                        i.max_month    <- find_max_consecutive_sum(i.month_table) 
                        i.focal.months <- names(i.max_month$values)
                        
                        #- create subset of i.data only containing the focal months
                        i.data <- i.data[month %in% i.focal.months]
                        
                } else {
                        #- Are the two months consecutive? 
                        i.monthdiff <- diff(as.numeric(names(i.month_table)))
                        #- No or only one month
                        if (length(i.monthdiff) == 0|i.monthdiff > 2){
                                # Does one of the month have more than 100 samples? 
                                if (any(i.month_table > 100)){
                                        i.focal.months <- names(which.max(i.month_table))
                                        i.data <- i.data[month == i.focal.months]
                                } else {
                                        rm(list = ls()[pattern="^i\\.", x = ls()])
                                        next()
                                }
                        } else { #- yes 
                                i.focal.months <- names(i.month_table)
                                i.data <- i.data[month %in% i.focal.months]
                        }
                } 
                
                #- how many samples are left? 
                i.samples <- uniqueN(i.data$comb_site_id)
                if (i.samples < 100) {
                        rm(list = ls()[grepl("^i\\.", x = ls())])
                        next()
                }
                # ———— What sorts of random combinations can be created from that? 
                i.hundreds    <- floor(i.samples/ 100)
                i.sample_scheme <- c()
                for (crs in 1:i.hundreds){
                        i.sample_scheme %<>% append(rep(crs * 100, times = min(i.hundreds - (crs-1),3)))
                }
                rm(crs)
                if ((i.samples/100 - i.hundreds) > 0.5){
                        i.sample_scheme %<>% append(max(i.sample_scheme))
                }
                rm(i.hundreds)
                #i.sample_scheme
                i.out <- data.table(data.set = ud[d], 
                           year     = unique(i.data$year),
                           samples  = i.sample_scheme)
                i.out[, focal_months := list(list(i.focal.months))]
                result_list[[length(result_list) + 1]] <- i.out
                rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
        }
        rm(list = ls()[grepl(pattern = "^ds\\.", x = ls())])
}
result.data <- rbindlist(result_list)

saveRDS(result.data, "data/schemes_fish.rds")
