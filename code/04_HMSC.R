# load data -------------------------------------------------------------------------
bio <- read_parquet("data/fish_w_env.parquet")
schemes <- readRDS("data/schemes_fish.rds")
# LOAD PARAMETERS -------------------------------------------------------------------
source("code/parameters/parameters_HMSC.R")
# prepare data ----------------------------------------------------------------------

#- loop over schemes
for (o in 1:nrow(schemes)) {
        set.seed(1)
        print(paste(o, "of", nrow(schemes)))
        o.scheme <- schemes[o, ]
        o.data   <- bio[data.set == o.scheme$data.set & 
                                year == o.scheme$year & 
                                month(date) %in% as.numeric(o.scheme$focal_months[[1]])]
        o.sample.ids <- sample(
                x = unique(o.data$sample_id),
                size = o.scheme$samples,
                replace = F
        )
        o.data <- o.data[sample_id %in% o.sample.ids]
        #- establish spatial scale
        o.sf <- determine_spatial_scale(o.data)
        
        # prepare data for HMSC ----------------------------------------------------------------------
        # HMSC needs a Y response matrix.
        # That matrix has one sample per row and one column per taxon
        # We need to reshape the data for this.
        o.data2 <- dcast(
                o.data,
                formula = sample_id ~ lowest.taxon,
                value.var = "PA",
                fill = 0
        )
        
        #- Check for rare taxa (less than 5% of sites)
        o.data3 <-
                o.data2[, {
                        cols_to_keep <- colnames(o.data2)[-1][colSums(o.data2[, -1] != 0) > (o.scheme$samples * 0.05)]
                        c(cols_to_keep)
                }, with = FALSE] %>%
                as.matrix
        
        
        #- prepare environmental variables
        o.env.samples <- unique(o.data, by = "sample_id")
        o.env.samples[, `:=` (sample_id = NULL, counter = NULL, comb_sample_id = NULL)]
        #- If glacial area is zero everywhere, this leads to issues later on.
        #- Therefore, this variable is removed in such cases.
        o.cols.to.keep <- o.env.samples[, names(which(lapply(.SD, uniqueN) > 1))]
        o.env.samples  <- o.env.samples[, ..o.cols.to.keep]
        rm(o.cols.to.keep)
        o.env.samples <- o.env.samples[, x.coord:valley_bottom_flattness]
        # Saturated soil water content has some missing values 
        # These are imputed here with a OLS Regression
        if (any(is.na(o.env.samples$saturated_soil_water_content))) {
                o.row <- which(is.na(o.env.samples$saturated_soil_water_content))
                o.pre <- o.env.samples[-o.row, ]
                o.pre <- o.pre[, comb_sample_id := NULL]
                o.mod <- lm(saturated_soil_water_content ~ ., data = o.pre)
                o.pre <- o.env.samples[, comb_sample_id := NULL] %>% .[o.row]
                o.pre <- predict(object = o.mod, newdata =  o.pre)
                o.env.samples[is.na(saturated_soil_water_content), saturated_soil_water_content := o.pre]
                rm(o.row, o.pre,o.mod)
        }
        #- scale predictor variables 
        o.env.samples.xy <- o.env.samples[, c("x.coord", "y.coord")]
        o.env.samples[, c("x.coord", "y.coord") := NULL]
        o.env.samples %<>% scale()
        attributes(o.env.samples)$`scaled:scale` <- NULL
        attributes(o.env.samples)$`scaled:center` <- NULL
        o.env.samples <- cbind(o.env.samples.xy, o.env.samples)
        
        # Asymmetric Eigenvector Maps -------------------------------------------------------
        # Load river network for the EU HYDRO DEM catchments in which samples are located.
        o.rivers <- load_river_networks(o.data)
        o.sites  <- o.data %>%
                unique(by = "sample_id") %>%
                st_as_sf(coords = c("x.coord", "y.coord"), crs = 3035) %>%
                st_transform(4326) 
        o.network <- prepare_directional_network(o.rivers, o.sites)
        # Calculate network costs
        o.cost_matrix <-
                st_network_cost(x       = o.network,
                                from    = o.sites,
                                to      = o.sites,
                                weights = "weight") %>%
                drop_units()
        o.cost_matrix <- create_directional_incidence(o.cost_matrix, o.sites)
        o.aem_result <- aem(binary.mat = o.cost_matrix$incidence,
                            weight     = o.cost_matrix$weights)
        #- select AEM vectors 
        o.signi.x <- c()
        for (focal.spe in 1:ncol(o.data3)){
                
                fs.env.samples <- o.env.samples[, -c(1, 2)]
                fs.env.samples <- cbind(o.data3[, focal.spe], fs.env.samples)
                fs.env.samples$V1 %<>% factor
                fs.predictors  <- names(fs.env.samples)[-1]
                fs.formulas    <- paste("V1 ~", paste(fs.predictors, collapse = "+"))
                fs.formulas    <- as.formula(fs.formulas)
                fs.model       <- glm(fs.formulas, data = fs.env.samples, family = "binomial")
                fs.residuals   <- residuals(fs.model)
                fs.p.values    <- apply(o.aem_result$vectors, 2, function(x) cor.test(x,fs.residuals)$p.value)
                fs.p.values    <- p.adjust(fs.p.values, method = "holm")
                if (any(fs.p.values < 0.05)){
                        fs.sp <- which(fs.p.values < 0.05)
                        o.signi.x <- append(o.signi.x, fs.sp)
                        o.signi.x <- unique(o.signi.x)
                }
                rm(list = ls()[grepl(pattern = "^fs\\.", x = ls())])
        }
        o.AEM     <- 
                o.aem_result$vectors[, o.signi.x] %>% 
                data.frame %>% 
                setDT
        names(o.AEM) <- paste0("AEM", 1:ncol(o.AEM))
        o.env.samples <- cbind(o.env.samples, o.AEM)
        
        #- Morans eigenvector maps 
        o.xy.mat <- matrix(
                c(o.env.samples$x.coord, o.env.samples$y.coord),
                ncol = 2,
                byrow = F
        )
        colnames(o.xy.mat) <- c("x", "y")
        o.x       <- dbmem(
                xyORdist    = o.xy.mat,
                MEM.autocor = "non-null",
                store.listw = T
        )
        #- Following Bini et al (2009) [10.1111/j.1600-0587.2009.05717.x]. In the paper this is called SEVM-3. 
        #- Loops over species 
        o.signi.x <- c()
        for (focal.spe in 1:ncol(o.data3)){
                
                fs.env.samples <- o.env.samples[, -c(1, 2)]
                fs.env.samples <- cbind(o.data3[, focal.spe], fs.env.samples)
                fs.env.samples$V1 %<>% factor
                fs.predictors  <- names(fs.env.samples)[-1]
                fs.formulas    <- paste("V1 ~", paste(fs.predictors, collapse = "+"))
                fs.formulas    <- as.formula(fs.formulas)
                fs.model       <- glm(fs.formulas, data = fs.env.samples, family = "binomial")
                fs.residuals   <- residuals(fs.model)
                fs.p.values    <- apply(o.x, 2, function(x) cor.test(x,fs.residuals)$p.value)
                fs.p.values    <- p.adjust(fs.p.values, method = "holm")
                if (any(fs.p.values < 0.05)){
                        fs.sp <- which(fs.p.values < 0.005)
                        o.signi.x <- append(o.signi.x, fs.sp)
                        o.signi.x <- unique(o.signi.x)
                }
                rm(list = ls()[grepl(pattern = "^fs\\.", x = ls())])
        }
        rm(focal.spe)
        #o.signi   <- moran.randtest(o.x, nrepet = 999)
        #o.test.id <- which(o.signi$pvalue == min(o.signi$pvalue))
        o.MEM     <- o.x[, o.signi.x] %>% 
                data.frame %>%
                setDT
        o.env.samples <- cbind(o.env.samples, o.MEM)
        o.n.env       <- ncol(o.env.samples) - 2 - length(o.signi.x) - ncol(o.AEM)
        
        #- prepare HMSC model
        #- establish a site level random factor
        o.studyDesign <- data.frame(sample = as.factor(1:o.scheme$samples))
        o.rL          <- HmscRandomLevel(units = o.studyDesign$sample)
        
        #- create model formula
        o.env.samples <- o.env.samples[, -c(1, 2)]
        o.predictors  <- names(o.env.samples)
        o.formulas    <- paste("~", paste(o.predictors, collapse = "+"))
        o.formulas    <- as.formula(o.formulas)

        #- define HMSC model
        o.mod1 <- Hmsc(
                Y           = o.data3,
                XData       = o.env.samples,
                XFormula    = o.formulas,
                studyDesign = o.studyDesign,
                ranLevels   = list("sample" = o.rL),
                distr       = "probit"
                
        )
        

        # FIT HMSC MODEL --------------------------------------------------------------------
        while.condition <- TRUE
        o.counter = 1
        while (while.condition) {
                print(paste("fitting model", o.counter))
                # ———— fit model —————
                o.m1.mcmc = sampleMcmc(
                        o.mod1,
                        thin      = thin * o.counter,
                        samples   = samples * o.counter,
                        transient = transient * o.counter,
                        nChains   = nChains,
                        nParallel = nChains
                )
                
                # ————— evaluate model fit ————
                #- check model fit with the potential scale reduction factor
                o.mpost <- convertToCodaObject(o.m1.mcmc)
                #- purpose-build function (see 00_setup.R)
                o.means      <- gelman_check(o.mpost)
                o.means.rate <- sum(o.means >= 1.1) / length(o.means) 
                if (o.means.rate >= 0.1) {
                        while.condition <- TRUE
                        o.counter <- o.counter + 1
                } else {
                        while.condition <- FALSE
                }
                if (o.counter == 6) {
                        while.condition == F
                        o.model.fit = FALSE
                }
        }        
        
        # variation partitioning ------------------------------------------------------------
        # -- for environmental vs space
        o.preds     <- computePredictedValues(o.m1.mcmc)
        o.MF        <- evaluateModelFit(hM = o.m1.mcmc, predY = o.preds)
        # how many spatial predictors are there?
        o.n.spatial <- sum(str_detect(names(o.env.samples), "AEM|MEM"))
        o.VP    <- computeVariancePartitioning(
                o.m1.mcmc,
                group = c(rep(
                        1, length(o.predictors) - o.n.spatial
                ), rep(2, o.n.spatial)),
                groupnames = c("env", "space")
        )
        o.VP2   <- data.table(
                taxon = rep(colnames(o.VP$vals), each = 3),
                driver = rep(c("env", "space", "bio"), times = ncol(o.VP$vals)),
                value  = c(o.VP$vals),
                r2     = rep(o.MF$TjurR2, each = 3)
        )
        o.VP2[r2 < 0, r2 := 0]
        o.VP2[, scaled_values := value * r2]
        #- add stochasticity = 1-r2
        o.VP3 <- data.table(
                taxon = colnames(o.VP$vals),
                driver = "stochastic",
                value = 0,
                r2 = o.MF$TjurR2
        )
        o.VP3[r2 < 0, r2 := 0]
        o.VP3[, value := 1 - r2]
        o.VP3[, scaled_values := 1 - r2]
        o.VP4 <- rbindlist(list(o.VP2, o.VP3))
        o.VP4$run <- o
        #o.VP4$psrf_mean <- rep(o.means, each = 4)
        #o.VP4$psrf_max  <- rep(o.maxes, each = 4)
        
        #- new VP to determine relative predictor importance
        o.VP5    <- computeVariancePartitioning(o.m1.mcmc)
        
        o.VP6 <- rowSums(o.VP5$vals)
        #- drop morans eigenvectors
        o.VP6 <- o.VP6[str_detect(names(o.VP6), "MEM|AEM", negate =
                                          T)]
        #- drop random sample
        o.VP6 <- o.VP6[str_detect(names(o.VP6), "Random", negate =
                                          T)]
        o.VP6 <- o.VP6 / sum(o.VP6)
        #o.VP6 <- rank(o.VP6)
        
        #o.groupnames = c("env", "space")
        
        #- save fitted model
        o.out.list <-
                list(
                        biota       = o.data3,
                        environment = o.env.samples,
                        model       = o.m1.mcmc,
                        VP          = o.VP4,
                        coordinates = o.xy.mat,
                        scale       = o.sf,
                        importance  = o.VP6,
                        scheme      = o.scheme
                )
        # print_i <- ifelse(i <10, paste0(0,i), i)
        saveRDS(o.out.list,
                paste0("data/fitted_models/", data.set, "_", i, "_", o, ".rds"))
        rm(list = ls()[grepl("^o\\.", x = ls())])
} # END OF o LOOP
rm(list = ls()[grepl("^i\\.", x = ls())])




