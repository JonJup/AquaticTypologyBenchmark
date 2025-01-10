# load data -------------------------------------------------------------------------
bio <- read_parquet("data/fish_w_env.parquet")

# LOAD PARAMETERS -------------------------------------------------------------------
source("code/parameters/parameters_HMSC.R")

#- currently there are many NAs in the slope data; needs to be fixd 
#bio[, slope := NULL]


# set loop independent parameters  --------------------------------------------------
n.data.sets      <- uniqueN(bio$data.set)
sampled.data.set <- sample(1:n.data.sets, size = 30000, replace = TRUE)

# prepare data ----------------------------------------------------------------------
for (i in 1:3){
        #- updater 
        print(i)
        # ——— create random sample ——— #
        #- which data set
        i.bio.sample <- bio[data.set == unique(data.set)[sampled.data.set[i]]]
        if(all(is.na(i.bio.sample$abundance))) next()
        i.n.samples <- min(uniqueN(i.bio.sample$comb_sample_id), n.samples)
        i.sample.ids <- sample(
                x = unique(i.bio.sample$comb_sample_id), 
                size = i.n.samples,
                replace = F
                )
        i.bio.sample <- i.bio.sample[comb_sample_id %in% i.sample.ids]
        #- establish spatial scale 
        i.sf <- unique(i.bio.sample, by = "original_site_name")
        i.sf <- st_as_sf(i.sf, coords = c("x.coord", "y.coord"), crs = 3035)
        i.sf <- st_distance(i.sf)
        i.sf <- drop_units(i.sf)
        i.sf <- i.sf[lower.tri(i.sf)]
        i.sf <- list(min = min(i.sf), max = max(i.sf), mean = mean(i.sf), median = median(i.sf))
        # prepare HMSC ----------------------------------------------------------------------
        #- HMSC needs a Y response matrix. That matrix has one sample per row and one column per taxon
        #- I need to reshape the data for this. 
        i.bio.sample2 <- dcast(i.bio.sample, 
                               formula = comb_sample_id ~ lowest.taxon, 
                               value.var = "PA", 
                               fill = 0)
        
        #- check for rare taxa (less than 5% of sites)
        i.bio.sample3 <- apply(
                i.bio.sample2[,-1], 2, 
                FUN = function(x) sum(x!=0)
                )
        i.bio.sample3 <- names(which(i.bio.sample3 > 5))
        i.bio.sample3 <- c("comb_sample_id", 
                           i.bio.sample3)   
        i.bio.sample4 <- i.bio.sample2[, i.bio.sample3, with = FALSE]
        i.bio.sample4 <- as.matrix(i.bio.sample4[,-1])
        
        #- prepare environmental variables 
        i.env.samples <- unique(
                i.bio.sample, 
                by = "comb_sample_id"
                )
        i.env.samples[, `:=` (comb_sample_id = NULL, 
                              counter = NULL)
                      ]
        i.env.samples <- i.env.samples[, x.coord:spi]
        
        
        #- MEMs 
        i.xy.mat <- matrix(
                c(i.env.samples$x.coord, 
                  i.env.samples$y.coord), 
                ncol = 2, 
                byrow = F
                )
        colnames(i.xy.mat) <- c("x", "y")
        i.x <- dbmem(xyORdist = i.xy.mat)
        i.signi <- moran.randtest(i.x, nrepet = 999)
        i.test.id <- which(i.signi$pvalue== min(i.signi$pvalue))
        #if (length(i.test.id > 5)) i.test.id <- c(1:5)
        i.MEM <- i.x[, i.test.id]

        setDT(i.MEM)
        i.env.samples <- cbind(i.env.samples, i.MEM)
        i.n.env <- ncol(i.env.samples) - 2 - length(i.test.id)
        #- prepare HMSC model 
        #- establish a site level random factor 
        i.studyDesign <- data.frame(sample = as.factor(1:i.n.samples) )
        i.rL          <- HmscRandomLevel(units = i.studyDesign$sample)
        
        #- create model formula
        i.env.samples <- i.env.samples[,-c(1,2)]
        i.predictors  <- names(i.env.samples)
        i.formulas    <- paste("~", paste(i.predictors, collapse = "+"))
        i.formulas    <- as.formula(i.formulas)
        
        #- define HMSC model 
        i.mod1 <- Hmsc(
                Y           = i.bio.sample4,
                XData       = i.env.samples,
                XFormula    = i.formulas,
                studyDesign = i.studyDesign,
                ranLevels   = list("sample" = i.rL),
                distr       = "probit"
                
        )
        # ———— fit model —————
        i.m1.mcmc = sampleMcmc(
                i.mod1,
                thin = thin,
                samples = samples,
                transient = transient,
                nChains = nChains,
                nParallel = nChains
        )
        
        # ————— evaluate model fit ————
        #- check model fit with the gelman diagnostic
        i.mpost <- convertToCodaObject(i.m1.mcmc)
        i.ge    <- gelman.diag(x = i.mpost$Beta)
        #- species names
        i.species_names  <- i.ge$psrf%>%rownames%>%stringr::str_extract(pattern = ",\\ .*\\(S") %>% str_remove(",") %>% str_remove("\\(S")%>%str_trim%>%unique
        # Create empty vectors to store results
        i.means <- numeric(length(i.species_names))
        i.maxes <- numeric(length(i.species_names))
        # Calculate statistics for each species
        for(j in seq_along(i.species_names)) {
                # Get rows corresponding to current species
                j.species_rows <- grep(i.species_names[j], rownames(i.ge$psrf))
                # Calculate mean and max of point estimates for these rows
                i.means[j] <- mean(i.ge$psrf[j.species_rows, "Point est."])
                i.maxes[j] <- max( i.ge$psrf[j.species_rows, "Point est."])
                rm(list = ls()[grepl("^j\\.", x = ls())])
        }
        rm(j)
        # ——— Variation Partitioning ————
        # -- for environmental vs space
        i.preds <- computePredictedValues(i.m1.mcmc)
        i.MF    <- evaluateModelFit(hM = i.m1.mcmc, predY = i.preds)
        i.VP    <- computeVariancePartitioning(i.m1.mcmc,
                                               group = c(
                                                       rep(1, length(i.predictors) - length(i.test.id)),
                                                       rep(2, length(i.test.id))
                                                       ),
                                               groupnames = c("env", "space")
        )
        i.VP2   <- data.table(
                taxon = rep(colnames(i.VP$vals), each = 3),
                driver = rep(c("env", "space", "bio"), times = ncol(i.VP$vals)),
                value  = c(i.VP$vals),
                r2     = rep(i.MF$TjurR2, each = 3)
        )
        i.VP2[r2 < 0, r2 := 0]
        i.VP2[,scaled_values := value * r2]
        #- add stochasticity = 1-r2
        i.VP3 <- data.table(taxon = colnames(i.VP$vals), 
                            driver = "stochastic",
                            value = 0,
                            r2 = i.MF$TjurR2)
        i.VP3[r2<0, r2 := 0]
        i.VP3[,value := 1-r2]
        i.VP3[,scaled_values := 1-r2]
        i.VP4 <- rbindlist(list(i.VP2, i.VP3))
        i.VP4$run <- i
        i.VP4$psrf_mean <- rep(i.means, each = 4)
        i.VP4$psrf_max  <- rep(i.maxes, each = 4)
        
        #- new VP to determine relative predictor importance 
        i.VP5    <- computeVariancePartitioning(i.m1.mcmc)
        
        i.VP6 <- rowSums(i.VP5$vals)
        #- drop morans eigenvectors 
        i.VP6 <- i.VP6[str_detect(names(i.VP6), "MEM", negate =T )]
        #- drop random sample 
        i.VP6 <- i.VP6[str_detect(names(i.VP6), "Random", negate =T )]
        i.VP6_total <- sum(i.VP6)
        i.VP6 <- rank(i.VP6)
        
        
        groupnames = c("env", "space")
        
        #- save fitted model 
        out.list <- 
                list(biota       = i.bio.sample4,
                     environment = i.env.samples,
                     model       = i.m1.mcmc,
                     VP          = i.VP4, 
                     coordinates = i.xy.mat,
                     scale       = i.sf,
                     importance  = i.VP6
                )
        
        print_i <- ifelse(i <10, paste0(0,i), i)
        saveRDS(out.list, paste0("data/fitted_models/", print_i, ".rds"))
        rm(list = ls()[grepl("^i\\.", x = ls())])
}
