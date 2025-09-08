# Simulate new data
#' Explantion of variables and functions 
#' VARIABLES 
# setup ------------------------------------------------------------------
# library(groundhog)
# pkgs <- c(
#         "cluster",
#         "data.table",
#         "dplyr",
#         "fs",
#         "Hmsc",
#         "kernlab",
#         "vegclust",
#         "magrittr",
#         "purrr",
#         "stringr"
# 
# )
# groundhog.library(pkgs,'2024-12-01')
# rm(pkgs)

library(data.table)
library(kernlab)
library(cluster)
library(vegclust)
library(Hmsc)
setwd(rstudioapi::getActiveProject())
source("code/parameters/run_parameters.R")
# setup loop ------------------------------------------------------------------------
model.files <- list.files(paste0("data/fitted_hmsc_models/"), full.names = T)
vp.files    <- list.files(paste0("data/variation_partitioning/"), full.names = T)
spatial.scale <- list.files(paste0("data/spatial_scale/"), full.names = T)
source("code/functions/balance_clusters.R")
set.seed(1)
# Loop ------------------------------------------------------------------------------
# this could be parallelized
for (i in 3:length(model.files)){
        
        #- Print Status Update to console 
        print(paste("i =",i))
        
        #- extract taxon from model name 
        i.taxon <- sub(x = model.files[i], pattern = "data/fitted_hmsc_models/", replacement = "")
        i.taxon <- sub(x =i.taxon        , pattern = "_.*\\.rds", replacement = "")
        
        
        #- Prepare list to store results of loop
        i.out <- vector(mode = "list", length = length(n_types))
        
        ## Verify that the names are the same 
        if (sub(x = model.files[i], pattern = "data/fitted_hmsc_models/", replacement = "") != sub(x = vp.files[i], pattern = "data/variation_partitioning/", replacement = "")){
                stop(paste("loaded files have different names in", i))
        }
        
        #- Load model results 
        i.model <- readRDS(model.files[i])
        #- Load the results of the variation partitioning 
        i.vp   <- readRDS(vp.files[i])
        
        #- Extract name of model 
        i.model.name <- sub("data/fitted_hmsc_models/", "", model.files[i])
                
        
        #- Extract number of samples in model
        i.nrow  <- nrow(i.model[[2]])
        #- Which variables follow the typology?
        #- Here we make 5 subsets which contain: Four times random assortments and one time 
        #- all variables. 
        #- The strength of the typology is partly determined by how important the variables are. 
        #- This is measured by the fraction of VP scores captured by the included variables. 
        #- These vectors hold information on the simulations 
        i.n.variables <- 
                i.contraction.points <- 
                i.contraction.centroids <- 
                i.importance  <- 
                i.asw <- c()
        
        i.cluster.assignments <- list()
        i.fuzzy.assignments   <- list()
        
        #--- loop over different number of types 
        for (j in 1:length(n_types)){
                
                #- print iteration number to console 
                print(paste("i =",i, "j =",j))
                
                #- assign number of types for this iteration 
                j.types <- n_types[j]
                
                #  — — — NEW PREDCICTOR VALUES  —  —  —  —
                #- extract predictor names 
                j.all.vars1 <- colnames(i.model$XData)[-c(1, which(colnames(i.model$XData) == "."))]
                #- remove spatial predictors (MEM = MORANS EIGENVECTOR MAPS) 
                j.all.vars1 <- j.all.vars1[!grepl("MEM|AEM",  j.all.vars1)]
                #- create a list to store the results of different predictor variable sets
                j.predictions <- list()
                
                #- START LOOP OVER Q, varies the number of variables that 
                #- follow the typology
                for (q in 1:max.q){
                        print(paste("model =", i, "number of types =", j, "number of variables =", q))
                        
                        # How many variables are included in this subset? 
                        # This is a random draw
                        q.nvariables <- sample(
                                x = 3:(length(j.all.vars1) - 1),
                                size = 1
                                )
                        # Add number to vector for later reference
                        i.n.variables <- append(i.n.variables, q.nvariables)
                        # Which variables are included 
                        # This is a random sample from all variables 
                        q.all.vars   <- sample(
                                j.all.vars1, 
                                size = q.nvariables, 
                                replace = F
                                )
                        # vector with names of all variables that do not follow 
                        # the typology
                        q.missing  <- j.all.vars1[which(!j.all.vars1 %in% q.all.vars)]
                        
                        # At this point we can already determine the typology 
                        # strength as judged by the selection of variables.
                        q.importance <- sum(i.vp$importance[q.all.vars])
                        i.importance <- append(i.importance, q.importance) 
        
                        #- Extract original environmental variables from HMSC model 
                        q.env <- copy(i.model$XData)
                        if (":" %in% names(q.env)) {
                                q.env <- q.env[, - which(names(q.env) == ".")]
                        }
                        #- Create a classification of sites based on the selected
                        #- environmental variables   
                        q.cluster.env <- copy(q.env)
                        q.cluster.env <- as.matrix(as.data.frame(q.cluster.env[, q.all.vars]))
                        # q.cluster.env %<>%
                        #         select(all_of(q.all.vars)) %>% 
                        #         as.data.frame %>% 
                        #         as.matrix
                        #- spectral clustering was chosen as allows for a balanced
                        #- number of sites among clusters  
                        test.out <- c()
                        #difer.out <- c()
                        for (numberOfClusters in 2:10){
                                test <- kmeans(x = q.cluster.env, 
                                       centers = numberOfClusters,
                                       nstart = 100,
                                       iter.max = 100)
                                if (any(test$size < 10)) next()
                                test2 <- silhouette(test$cluster, dist(q.cluster.env))
                                test.out[length(test.out) + 1] <- mean(test2[,3])
                                names(test.out)[length(test.out)] <- paste(numberOfClusters)
                                #sizes <- test$size
                                #print(sizes)
                        }     
                        opt.nc <- which.min(rank(rank(difer.out) - rank(test.out)))

                        
                        
                        q.type <- specc(x       = q.cluster.env, 
                                        centers = j.types)
                        q.clusters <- balance_clusters(
                                data    = q.cluster.env,
                                clusters = q.type@.Data,
                                min_size = nrow(q.cluster.env) / j.types * 0.75
                        ) 

                        
                        #- store cluster assignment for later 
                        i.cluster.assignments[[length(i.cluster.assignments) + 1]] <- q.clusters
                        #j.fuzzy.assignments[[q]]   <- q.fc.b
                        setDT(q.env)
                        q.env[, type := q.type]
                        #q.env[, type.fuzzy := q.fc.b]
                        #- find the centroid for each variable that follows the typology
                        q.centroids <- q.env[, lapply(.SD, mean), .SDcols = q.all.vars, by = "type"]
                        setorderv(q.centroids, "type")
                        q.centroids[, type := NULL]
                        q.centroids <- as.matrix(as.data.frame(q.centroids))
                        
                        #- prepare environmental variables for adjustments in script below 
                        q.observations <- as.matrix(as.data.frame(q.cluster.env))
                        
                        # Based on two vectors of how strongly the centroids and sites are contracted
                        # recompute environmental varibles 
                        source("code/helper/new_predictors_for_hard_classification.R")
                        
                        # fuzzy classification 
                        q.fuzzy <- lapply(
                                1:within.q, 
                                function(x) {
                                        vegclust(q.newenv[[x]][, which(colnames(q.newenv[[x]]) %in% q.all.vars)],
                                        mobileCenters = j.types,
                                        method = "FCM",
                                        m = 1.5
                                        )
                                }        
                        )
                        i.fuzzy.assignments[[length(i.fuzzy.assignments) + 1]] <- q.fuzzy
                        #- There used to be a balancing step here. 
                        #- I removed it.
                        #- Since the spectral clustering is already balanced, and the fuzzy classification uses 
                        #- the results of the spectral classification, it is likely not necessary.
                        # Comine MCMC chains for prediction 
                        q.posterior_samples <- poolMcmcChains(i.model$postList)
                        q.selected_samples  <- q.posterior_samples[c(1,length(q.posterior_samples))]
                        
                        #  — — — predict new biotic communities  —  —  —  —
                        q.out <- lapply(q.newenv, 
                                        FUN = function(pre) {
                                                x <- predict(object = i.model, 
                                                        post = q.selected_samples, 
                                                        X = pre
                                                    )
                                                x <- x[[2]]
                                                return(x)
                                        }
                                        )
                        j.predictions [[q]] <- q.out
                        rm(list = ls()[grepl("^q\\.", x = ls())])
                }
                rm(q)
                i.out[[j]] <- j.predictions
                rm(list = ls()[grepl("^j\\.", x = ls())])
        }
        rm(j)
        i.out <- lapply(i.out, unlist, recursive = FALSE)
        i.out <- unlist(i.out, recursive = FALSE)  
        
        # read in spatial scale for additional information
        i.ss <- readRDS(spatial.scale[i])
        out_list <- 
                list(
                        data = i.out,
                        number_of_variables = i.n.variables,
                        contraction_points  = i.contraction.points,
                        contraction_centroids= i.contraction.centroids, 
                        variable_importance = i.importance,   
                        asw = i.asw, 
                        hard_cluster_assignment = i.cluster.assignments, 
                        fuzzy_cluster_assignment = i.fuzzy.assignments, 
                        info = i.ss[,c(1,2,3,4)]
                )
        # name.i <- ifelse(i < 10, paste0("0",i), as.character(i))
        i.save.name <- sub(x = model.files[i], pattern = "fitted_hmsc_models", replacement = "simulated_data")
        saveRDS(out_list, i.save.name)
        rm(out_list)
        rm(list = ls()[grepl("^i\\.", x = ls())])
} 
rm(i)


# rm(i)
# out.all <- rbindlist(out)

# out.all %>% 
#         filter(metric == "cs") %>% 
#         ggplot(aes(y = value , x = importance)) + geom_point(aes(col = contraction)) + geom_smooth()
# 
# saveRDS(out.all, "data/250115_results.rds")
# rm(out.all, out, model.files, n_types)

