# load data -------------------------------------------------------------------------
bio <- read_parquet("data/fish_w_env.parquet")

# LOAD PARAMETERS -------------------------------------------------------------------
source("code/parameters/parameters_HMSC.R")

# set loop independent parameters  --------------------------------------------------
n.data.sets      <- uniqueN(bio$data.set)
#sampled.data.set <- sample(1:n.data.sets, size = 30000, replace = TRUE)

# prepare data ----------------------------------------------------------------------
data.set.sizes <- c()
for (data.set in 1:n.data.sets){
        
        #if (data.set < 4) next()
        
        #- select data set 
        ds.data.set <- unique(bio$data.set)[data.set]
        ds.bio      <- bio[data.set == ds.data.set]
        #- how many years? 
        ds.n.year <- uniqueN(ds.bio$year)
        #- how many samples? 
        ds.n.samples <- ds.bio[,uniqueN(comb_sample_id), by = "year"]
        
        if(all(ds.n.samples$V1 < 100)){
                print(paste("data set", data.set, "contains less than 100 samples in each year")) 
                next()
        } 
        #- select years with more than 100 samples 
        ds.n.samples <- ds.n.samples[V1 >= 100]
        for (i in 1:nrow(ds.n.samples)){ # START LOOP i OVER years in data set
                #- updater 
                #print(paste("data set: ", data.set, "; year: ", i, "of", nrow(ds.n.samples)))
                
                #- Create subset of focal year. 
                i.data <- ds.bio[year == ds.n.samples$year[i]]
        
                #- check seasons. To prevent strong seasonal changes from influencing the 
                #- we identify the three consecutive month with the most samples.
                #- Create month variable
                i.data[, month := month(date)]
                i.month_table <- unique(i.data, by = "comb_site_id")
                #- how often does each month occur
                i.month_table <- table(i.month_table$month)
                
                #- Which three consecutive month have the most samples? 
                #- This function is defined in 00_setup.R
                if (length(i.month_table) > 2){
                        i.max_month <- find_max_consecutive_sum(i.month_table) 
                        #- Create subset only containing these months
                        i.data <- i.data[month %in% names(i.max_month$values)]
                } else {
                        #- Are the two months consecutive? 
                        i.monthdiff <- diff(as.numeric(names(i.month_table)))
                        #- No or only one month
                        if (length(i.monthdiff) == 0|i.monthdiff > 2){
                                # Does one of the month have more than 100 samples? 
                                if (any(i.month_table > 100)){
                                      i.data <- i.data[month == names(which.max(i.month_table))]
                                } else {
                                        next()
                                }
                        } else { #- yes 
                                i.data <- i.data[month %in% names(i.month_table)]
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
                        i.sample_scheme %<>% append(rep(crs * 100, times = i.hundreds - (crs-1)))
                }
                rm(crs)
                if ((i.samples/100 - i.hundreds) > 0.5){
                        i.sample_scheme %<>% append(max(i.sample_scheme))
                }
                rm(i.hundreds)
                #- loop over different random samples 
                for (o in 1:length(i.sample_scheme)){
                        #if (i == 1 & o < 28) next()
                        print(paste("data set: ", data.set, "; year: ", i, "of", nrow(ds.n.samples), "scheme", o, "of",length(i.sample_scheme)))
                        o.sample.ids <- sample(
                                x = unique(i.data$comb_sample_id),
                                size = i.sample_scheme[o],
                                replace = F
                        )
                        o.data <- i.data[comb_sample_id %in% o.sample.ids]
                        #- establish spatial scale 
                        o.sf <- unique(o.data, by = "original_site_name")
                        o.sf <- st_as_sf(o.sf, coords = c("x.coord", "y.coord"), crs = 3035)
                        o.sf <- st_distance(o.sf)
                        o.sf <- drop_units(o.sf)
                        o.sf <- o.sf[lower.tri(o.sf)]
                        o.sf <- list(min = min(o.sf), max = max(o.sf), mean = mean(o.sf), median = median(o.sf))
                        
                        # prepare HMSC ----------------------------------------------------------------------
                        #- HMSC needs a Y response matrix. 
                        #- That matrix has one sample per row and one column per taxon
                        #- We need to reshape the data for this. 
                        o.data2 <- dcast(
                                o.data, 
                                formula = comb_sample_id ~ lowest.taxon, 
                                value.var = "PA", 
                                fill = 0
                                )
                        
                        #- Check for rare taxa (less than 5% of sites)
                        o.data3 <- apply(
                                o.data2[,-1], 2, 
                                FUN = function(x) sum(x!=0)
                        )
                        o.data3 <- names(
                                which(o.data3 > round(i.sample_scheme[o]/100 * 5))
                                )
                        o.data3 <- c("comb_sample_id", 
                                     o.data3)   
                        o.data4 <- o.data2[, o.data3, with = FALSE]
                        o.data4 <- as.matrix(o.data4[,-1])
                        
                        #- prepare environmental variables 
                        o.env.samples <- unique(
                                o.data, 
                                by = "comb_sample_id"
                        )
                        o.env.samples[, `:=` (comb_sample_id = NULL, 
                                              counter = NULL)
                        ]
                        #- If glacial area is zero everywhere, this leads to issues later on. 
                        #- Therefore, this variable is removed in such cases. 
                        if (all(o.env.samples$glacial_area == 0)){
                                o.env.samples[, glacial_area := NULL]
                        }
                        if (all(o.env.samples$Rfactor_min == 0)){
                                o.env.samples[, Rfactor_min := NULL]
                        }
                        
                        o.env.samples <- o.env.samples[, x.coord:valley_bottom_flattness]
                        # TODO cheap imputation for saturated soil water content. Needs to improve! 
                        if (any(is.na(o.env.samples$saturated_soil_water_content))){
                                o.mean <- mean(o.env.samples$saturated_soil_water_content, na.rm = T)
                                o.env.samples[is.na(saturated_soil_water_content), saturated_soil_water_content := o.mean]
                        }
                        o.env.samples.xy <- o.env.samples[, c("x.coord", "y.coord")]
                        o.env.samples[, c("x.coord", "y.coord") := NULL]
                        o.env.samples %<>% scale()
                        attributes(o.env.samples)$`scaled:scale` <- NULL
                        attributes(o.env.samples)$`scaled:center` <- NULL
                        o.env.samples <- cbind(o.env.samples.xy, o.env.samples)
                        #- MEMs 
                        o.xy.mat <- matrix(
                                c(o.env.samples$x.coord, 
                                  o.env.samples$y.coord), 
                                ncol = 2, 
                                byrow = F
                        )
                        colnames(o.xy.mat) <- c("x", "y")
                        o.x       <- dbmem(
                                xyORdist    = o.xy.mat,
                                MEM.autocor = "non-null",
                                store.listw = T)
                        o.signi   <- moran.randtest(o.x, nrepet = 999)
                        o.test.id <- which(o.signi$pvalue== min(o.signi$pvalue))
                        o.MEM     <- o.x[, o.test.id]
                        o.MEM     <- data.frame(o.MEM)
                        setDT(o.MEM)
                        o.env.samples <- cbind(o.env.samples, o.MEM)
                        o.n.env       <- ncol(o.env.samples) - 2 - length(o.test.id)
                        
                        #- prepare HMSC model 
                        #- establish a site level random factor 
                        o.studyDesign <- data.frame(sample = as.factor(1:i.sample_scheme[o]) )
                        o.rL          <- HmscRandomLevel(units = o.studyDesign$sample)
                        
                        #- create model formula
                        o.env.samples <- o.env.samples[,-c(1,2)]
                        o.predictors  <- names(o.env.samples)
                        o.formulas    <- paste("~", paste(o.predictors, collapse = "+"))
                        o.formulas    <- as.formula(o.formulas)
                        
                        #- define HMSC model 
                        o.mod1 <- Hmsc(
                                Y           = o.data4,
                                XData       = o.env.samples,
                                XFormula    = o.formulas,
                                studyDesign = o.studyDesign,
                                ranLevels   = list("sample" = o.rL),
                                distr       = "probit"
                                
                        )
                        # ———— fit model —————
                        o.m1.mcmc = sampleMcmc(
                                o.mod1,
                                thin      = thin,
                                samples   = samples,
                                transient = transient,
                                nChains   = nChains,
                                nParallel = nChains
                        )
                        
                        # ————— evaluate model fit ————
                        #- check model fit with the potential scale reduction factor
                        o.mpost <- convertToCodaObject(o.m1.mcmc)
                        o.ge    <- gelman.diag(x = o.mpost$Beta)
                        #- species names
                        o.species_names  <- o.ge$psrf %>% 
                                rownames %>%
                                str_extract(pattern = ",\\ .*\\(S") %>% 
                                str_remove(",") %>% 
                                str_remove("\\(S") %>%
                                str_trim %>% 
                                unique
                        #- Create empty vectors to store results
                        o.means <- numeric(length(o.species_names))
                        o.maxes <- numeric(length(o.species_names))
                        #- Calculate statistics for each species
                        for(j in seq_along(o.species_names)) {
                                # Get rows corresponding to current species
                                j.species_rows <- grep(o.species_names[j], rownames(o.ge$psrf))
                                # Calculate mean and max of point estimates for these rows
                                o.means[j] <- mean(o.ge$psrf[j.species_rows, "Point est."])
                                o.maxes[j] <- max( o.ge$psrf[j.species_rows, "Point est."])
                                rm(list = ls()[grepl("^j\\.", x = ls())])
                        }
                        rm(j)
                        # ——— Variation Partitioning ————
                        # -- for environmental vs space
                        o.preds <- computePredictedValues(o.m1.mcmc)
                        o.MF    <- evaluateModelFit(hM = o.m1.mcmc, predY = o.preds)
                        o.VP    <- computeVariancePartitioning(o.m1.mcmc,
                                                               group = c(
                                                                       rep(1, length(o.predictors) - length(o.test.id)),
                                                                       rep(2, length(o.test.id))
                                                               ),
                                                               groupnames = c("env", "space")
                        )
                        o.VP2   <- data.table(
                                taxon = rep(colnames(o.VP$vals), each = 3),
                                driver = rep(c("env", "space", "bio"), times = ncol(o.VP$vals)),
                                value  = c(o.VP$vals),
                                r2     = rep(o.MF$TjurR2, each = 3)
                        )
                        o.VP2[r2 < 0, r2 := 0]
                        o.VP2[,scaled_values := value * r2]
                        #- add stochasticity = 1-r2
                        o.VP3 <- data.table(taxon = colnames(o.VP$vals), 
                                            driver = "stochastic",
                                            value = 0,
                                            r2 = o.MF$TjurR2)
                        o.VP3[r2<0, r2 := 0]
                        o.VP3[,value := 1-r2]
                        o.VP3[,scaled_values := 1-r2]
                        o.VP4 <- rbindlist(list(o.VP2, o.VP3))
                        o.VP4$run <- i
                        o.VP4$psrf_mean <- rep(o.means, each = 4)
                        o.VP4$psrf_max  <- rep(o.maxes, each = 4)
                        
                        #- new VP to determine relative predictor importance 
                        o.VP5    <- computeVariancePartitioning(o.m1.mcmc)
                        
                        o.VP6 <- rowSums(o.VP5$vals)
                        #- drop morans eigenvectors 
                        o.VP6 <- o.VP6[str_detect(names(o.VP6), "MEM", negate =T )]
                        #- drop random sample 
                        o.VP6 <- o.VP6[str_detect(names(o.VP6), "Random", negate =T )]
                        o.VP6 <- o.VP6/sum(o.VP6)
                        #o.VP6 <- rank(o.VP6)
                        
                        #o.groupnames = c("env", "space")
                        
                        #- save fitted model 
                        o.out.list <- 
                                list(biota       = o.data4,
                                     environment = o.env.samples,
                                     model       = o.m1.mcmc,
                                     VP          = o.VP4, 
                                     coordinates = o.xy.mat,
                                     scale       = o.sf,
                                     importance  = o.VP6,
                                     data.et     = unique(bio$data.set)[data.set],
                                     year        = ds.n.samples$year[i],
                                     size        = i.sample_scheme[o]
                                )
                        # print_i <- ifelse(i <10, paste0(0,i), i)
                        saveRDS(o.out.list, 
                                paste0("data/fitted_models/",
                                       data.set,"_",
                                       i,
                                       "_",
                                       o, 
                                       ".rds")
                                )
                        rm(list = ls()[grepl("^o\\.", x = ls())])
                } # END OF o LOOP
                rm(list = ls()[grepl("^i\\.", x = ls())])
        } # END OF i LOOP
        rm(list = ls()[grepl("^ds\\.", x = ls())])
}# END OF data.set LOOP



