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
source("code/functions/normalized_partitioning_entropy.R")
set.seed(1)
# Loop ------------------------------------------------------------------------------
# this could be parallelized
for (i in 1:length(model.files)) {
        #- Print Status Update to console
        print(paste("i =", i))
        
        #- extract taxon from model name
        i.taxon <- sub(x = model.files[i],
                       pattern = "data/fitted_hmsc_models/",
                       replacement = "")
        i.taxon <- sub(x = i.taxon        ,
                       pattern = "_.*\\.rds",
                       replacement = "")
        
        
        #- Prepare list to store results of loop
        i.out <- vector(mode = "list", length = length(n_types))
        
        ## Verify that the names are the same
        if (sub(x = model.files[i],
                pattern = "data/fitted_hmsc_models/",
                replacement = "") != sub(x = vp.files[i],
                pattern = "data/variation_partitioning/",
                replacement = "")
            ) {
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
                i.nc  <-
                i.asw <- c()
        
        i.cluster.assignments <- list()
        i.fuzzy.assignments   <- list()
        
        #  — — — NEW PREDCICTOR VALUES  —  —  —  —
        #- extract predictor names
        i.all.vars1 <- colnames(i.model$XData)[-c(1, which(colnames(i.model$XData) == "."))]
        #- remove spatial predictors (MEM = MORANS EIGENVECTOR MAPS)
        i.all.vars1 <- i.all.vars1[!grepl("MEM|AEM", i.all.vars1)]
        #- create a list to store the results of different predictor variable sets
        i.predictions <- list()
        
        #- START LOOP OVER Q, varies the number of variables that
        #- follow the typology
        for (q in 1:max.q) {
                print(paste("model ", i, "run ", q))
                
                # How many variables are included in this subset?
                # This is a random draw
                q.nvariables <- sample(x = 3:(length(i.all.vars1) - 1), size = 1)
                # Add number to vector for later reference
                i.n.variables <- append(i.n.variables, q.nvariables)
                # Which variables are included
                # This is a random sample from all variables
                q.all.vars   <- sample(i.all.vars1,
                                       size = q.nvariables,
                                       replace = F)
                # all variables that do not follow the typology
                q.missing  <- i.all.vars1[which(!i.all.vars1 %in% q.all.vars)]
                
                # At this point we can already determine the typology
                # strength as judged by the selection of variables.
                q.importance <- sum(i.vp$importance[q.all.vars])
                i.importance <- append(i.importance, q.importance)
                
                #- Extract original environmental variables from HMSC model
                q.env <- copy(i.model$XData)
                if ("." %in% names(q.env)) {
                        q.env <- q.env[, -which(names(q.env) == ".")]
                }
                #- Create a classification of sites based on the selected
                #- environmental variables
                q.cluster.env <- copy(q.env)
                q.cluster.env <- as.matrix(as.data.frame(q.cluster.env[, q.all.vars]))
                # find optimal k means clustering
                q.kmeans.asw <-  q.smallest.cluster <- q.n.small.clusters <- c()
                for (numberOfClusters in 2:10) {
                        nc.clust <- kmeans(
                                x = q.cluster.env,
                                centers = numberOfClusters,
                                nstart = 100,
                                iter.max = 100
                        )
                        q.smallest.cluster[numberOfClusters - 1] <- min(nc.clust$size)
                        q.n.small.clusters[numberOfClusters - 1] <- sum(nc.clust$size < 10)
                        if (any(nc.clust$size < 10)){
                                rm(list = ls()[grepl("^nc\\.", x = ls())])
                                next()
                        }
                                
                        nc.asw <- silhouette(nc.clust$cluster,
                                             dist(q.cluster.env))
                        q.kmeans.asw[length(q.kmeans.asw) + 1] <- mean(nc.asw[, 3])
                        names(q.kmeans.asw)[length(q.kmeans.asw)] <- paste(numberOfClusters)
                        rm(list = ls()[grepl("^nc\\.", x = ls())])
                }
                # if there is no solution where the smallest cluster has at least 10 samples we need to balance manually
                if (is.null(q.kmeans.asw)) {
                        # best available number
                        # which.max only returns one even if multiple solutions share the same number
                        names(q.smallest.cluster) <- 2:10
                        q.smallest.cluster <- q.smallest.cluster[q.n.small.clusters == 1]
                        q.nc <- as.numeric(names(q.smallest.cluster[max(which(
                                q.smallest.cluster == max(q.smallest.cluster)
                        ))]))
                        q.clusters <- kmeans(
                                x = q.cluster.env,
                                centers = q.nc,
                                nstart = 100,
                                iter.max = 100
                        )
                        q.clusters <- q.clusters$cluster
                        q.clusters <- balance_clusters(
                                data = q.cluster.env,
                                clusters = q.clusters,
                                min_size = 10
                        )
                        
                } else {
                        q.nc <- as.numeric(names(which.max(q.kmeans.asw)))
                        q.clusters <- kmeans(
                                x = q.cluster.env,
                                centers = q.nc,
                                nstart = 100,
                                iter.max = 100
                        )
                        q.clusters <- q.clusters$cluster
                }
                
                #                 q.type <- specc(x       = q.cluster.env, centers = j.types)
                #                 q.clusters <- balance_clusters(
                #                         data    = q.cluster.env,
                #                         clusters = q.type@.Data,
                #                         min_size = nrow(q.cluster.env) / j.types * 0.75
                #                 )
                
                
                #- store cluster assignment for later
                i.cluster.assignments[[length(i.cluster.assignments) + 1]] <- q.clusters
                i.nc[length(i.nc) + 1] <- q.nc
                setDT(q.env)
                q.env[, type := q.clusters]
                #- find the centroid for each variable that follows the typology
                q.centroids <- q.env[, lapply(.SD, mean), .SDcols = q.all.vars, by = "type"]
                setorderv(q.centroids, "type")
                q.centroids[, type := NULL]
                q.centroids <- as.matrix(as.data.frame(q.centroids))
                
                #- prepare environmental variables for adjustments in script below
                q.observations <- as.matrix(as.data.frame(q.cluster.env))
                
                # Based on two vectors of how strongly the centroids and sites are contracted
                # recompute environmental variables
                source("code/helper/new_predictors_for_hard_classification.R")
                
                # fuzzy classification
                q.fuzzy <- lapply(1:within.q, function(x) {
                        vegclust(
                                q.newenv[[x]][, which(colnames(q.newenv[[x]]) %in% q.all.vars)],
                                mobileCenters = q.nc,
                                method = "FCM",
                                m = 1.5
                        )
                })
                i.fuzzy.assignments[[length(i.fuzzy.assignments) + 1]] <- q.fuzzy
                # Judge quality of fuzzy classification. Analog to ASW for hard classification. 
                i.npe <- sapply(q.fuzzy, NPE)
                
                
                #- There used to be a balancing step here.
                #- I removed it.
                #- Since the spectral clustering is already balanced, and the fuzzy classification uses
                #- the results of the spectral classification, it is likely not necessary.
                # Comine MCMC chains for prediction
                q.posterior_samples <- poolMcmcChains(i.model$postList)
                q.selected_samples  <- q.posterior_samples[c(1, length(q.posterior_samples))]
                
                #  — — — predict new biotic communities  —  —  —  —
                q.out <- lapply(
                        q.newenv,
                        FUN = function(pre) {
                                x <- predict(
                                        object = i.model,
                                        post = q.selected_samples,
                                        X = pre
                                )
                                x <- x[[2]]
                                return(x)
                        }
                )
                i.out[[q]] <- q.out
                #i.predictions [[q]] <- q.out
                rm(list = ls()[grepl("^q\\.", x = ls())])
                
        }
        rm(q)
        #i.out <- lapply(i.out, unlist, recursive = FALSE)
        i.out <- unlist(i.out, recursive = FALSE)
        
        # read in spatial scale for additional information
        i.ss <- readRDS(spatial.scale[i])
        out_list <-
                list(
                        data = i.out,
                        number_of_clusters = i.nc,
                        number_of_variables = i.n.variables,
                        contraction_points  = i.contraction.points,
                        contraction_centroids = i.contraction.centroids,
                        variable_importance = i.importance,
                        asw = i.asw,
                        npe = i.npe,
                        hard_cluster_assignment = i.cluster.assignments,
                        fuzzy_cluster_assignment = i.fuzzy.assignments,
                        info = i.ss[, c(1, 2, 3, 4)]
                )
        # name.i <- ifelse(i < 10, paste0("0",i), as.character(i))
        i.save.name <- sub(x = model.files[i],
                           pattern = "fitted_hmsc_models",
                           replacement = "simulated_data")
        saveRDS(out_list, i.save.name)
        rm(out_list)
        rm(list = ls()[grepl("^i\\.", x = ls())])
}
rm(i)


