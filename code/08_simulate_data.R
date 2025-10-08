################################################################################
# Script Name:        simulate_data.R
# Description:        Short description of the script
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-15
# Last Modified:      2025-09-22
#
# R Version:          R 4.5.1
# Required Packages: data.table, kernlab, cluster, vegclust, Hmsc, isotree
#
# Notes:             
################################################################################


# setup -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())

# libraries  --------------------------------------------------------------
library(data.table)
library(kernlab)
library(cluster)
library(vegclust)
library(Hmsc)
library(isotree)

# scripts -----------------------------------------------------------------
source("code/parameters/run_parameters.R")
source("code/functions/balance_clusters.R")
source("code/functions/normalized_partitioning_entropy.R")

# load files ------------------------------------------------------------------------
model.files   <- list.files("data/converged_hmsc_models/", full.names = T)
vp.files      <- list.files("data/variation_partitioning/", full.names = T)
spatial.scale <- list.files("data/spatial_scale/", full.names = T)
schemes       <- list.files("data/biota", full.names = T)

schemes <- schemes[grepl("03", schemes)]
schemes <- lapply(schemes, readRDS)
schemes <- rbindlist(schemes)
#names(schemes) <- c("diatoms", "fish", "invertebrates", "macrophytes")

isomaha <- readRDS("data/isomaha.rds")
mcmcPar <- readRDS("data/mcmc_parameters.rds")



set.seed(1)
# Loop ------------------------------------------------------------------------------

for (i in 1) {
        

        
        # load the model file
        i.model <- readRDS(model.files[i])
        #TODO remove 
        fit.model    <- readRDS(model.files[i])[[1]]
        fit.model    <- try(from_json(fit.model))
        unfittedModel <- readRDS("data/unfitted_hmsc_models/fish_0006.rds")
        
        # extract model name 
        model.name <- gsub(pattern="data/converged_hmsc_models//","", model.files[i])
        model.name <- gsub(pattern="data/converged_hmsc_models/","", model.name)
        model.name <- gsub(pattern="\\.rds", "", model.name)
        
        #  Print Status Update to console
        print(paste(model.name))
        
        
        # select the corresponding scheme 
        i.schemes <- schemes[scheme_id == model.name]
        
        # select the corresponding parameter set 
        i.para <- mcmcPar[scheme_id == model.name]
        
        
        postList = fit.model[1:i.para$nChains]
        fitTF = importPosteriorFromHPC(
                m = unfittedModel,  # the Hmsc object containing the unfitted model
                postList = postList, 
                nSamples = i.para$nSamples, 
                thin = i.para$thin, 
                transient = i.para$transient)
        

        
        # input in spatial scale for additional information
        #i.ss <- grep(x = spatial.scale, pattern = model.name)
        #i.ss <- readRDS(spatial.scale[i])
                
        
        #  Load the results of the variation partitioning
        # i.vp   <- readRDS(vp.files[i])
        
       #  check if the model is any good (i.e. Median AUC > 0.75) 
        # if ("C.SR2" %in% names(i.vp$MF)){
        #         if (!median(i.vp$MF$O.AUC) > 0.75) next()
        # } else {
        #         if (!median(i.vp$MF$O.AUC) > 0.75) next()
        # }

        #  Prepare list to store results of loop
        i.out <- vector(mode = "list", length = length(n_types))
        
        # ## Verify that the names are the same
        # if (sub(x = model.files[i],
        #         pattern = "data/fitted_hmsc_models/",
        #         replacement = "") != sub(x = vp.files[i],
        #         pattern = "data/variation_partitioning/",
        #         replacement = "")
        #     ) {
        #         stop(paste("loaded files have different names in", i))
        # }
        
        # #  Load model results
        # i.model <- readRDS(model.files[i])
        
        # convergence_state = i.model$converged
        # i.model <- i.model$model
        #  Extract name of model
        # i.model.name <- sub("data/fitted_hmsc_models/", "", model.files[i])
        
        
        #  Extract number of samples in model
        # i.nrow  <- nrow(i.model[[2]])
        #  Which variables follow the typology?
        #  Here we make 5 subsets which contain: Four times random assortments and one time
        #  all variables.
        #  The strength of the typology is partly determined by how important the variables are.
        #  This is measured by the fraction of VP scores captured by the included variables.
        #  These vectors hold information on the simulations
        i.n.variables <-
                i.contraction.points <-
                i.contraction.centroids <-
                i.importance  <-
                i.nc  <-
                i.npe <-
                i.asw <- c()
        
        i.cluster.assignments <- 
        i.fuzzy.assignments   <- list()
        
        #  — — — NEW PREDCICTOR VALUES  —  —  —  —
        #  extract predictor names
        i.all.vars1 <- colnames(i.model$XData)
        if (any (colnames(i.model$XData) == ".")){
                i.all.vars1 <- i.all.vars1[-which(colnames(i.model$XData) == ".")]
        }
        #  remove spatial predictors (MEM = MORANS EIGENVECTOR MAP)
        i.all.vars1 <- i.all.vars1[!grepl("MEM", i.all.vars1)]
        
        #  sort i.all.vars alphabetically. This matters for assigning the correct probabilities
        #  in the sampling. 
        i.all.vars1 <- sort(i.all.vars1)

        
        #  create a list to store the results of different predictor variable sets
        i.predictions <- list()
        
        #  Create a names vector to store how often each variable has been included in the artificial typologies. 
        #  The probability that any variable is chosen is inversely proportional to the number of times it has been included
        #  In successful artificial typologies. 
        
        i.variable.counter <- rep(1, length(i.all.vars1))
        # if (!"glacial_area" %in% i.all.vars1) i.variable.counter <- i.variable.counter[-1]
        # if (!"mean_snow_equivalent" %in% i.all.vars1) i.variable.counter <- i.variable.counter[-1]
        
        names(i.variable.counter) <- i.all.vars1
        #  START LOOP OVER Q, varies the number and identity of variables that
        #  follow the typology
        for (q in 1:max.q) {
                
                if (q == 1) i.success <- 0
                if (i.success > 200) next()
                print(paste(q, i.success))
                #print(paste("model ", i, "run ", q))
                
                # set up probabilities 
                q.prob <- 1 / i.variable.counter
                
                
                # How many variables are included in this subset?
                # This is a random draw
                q.nvariables <- sample(x = 3:(length(i.all.vars1) - 1), size = 1)
                # Which variables are included? This is a random sample from all variables
                q.all.vars   <- sample(i.all.vars1,
                                       size = q.nvariables,
                                       replace = F, 
                                       prob = q.prob)
                # if any of the geology variables are included, all must be included. 
                if (("area_calcareous" %in% q.all.vars | "area_sediment" %in% q.all.vars | "area_siliceous" %in% q.all.vars & 
                     !all(c("area_calcareous", "area_sediment", "area_siliceous") %in% q.all.vars))){
                        if (! "area_calcareous" %in% q.all.vars) q.all.vars <- append(q.all.vars, "area_calcareous")
                        if (! "area_sediment" %in% q.all.vars)   q.all.vars <- append(q.all.vars, "area_sediment")
                        if (! "area_siliceous" %in% q.all.vars)  q.all.vars <- append(q.all.vars, "area_siliceous")
                }
                
                

                # all variables that do not follow the typology
                q.missing  <- i.all.vars1[which(!i.all.vars1 %in% q.all.vars)]
                
                #  Extract original environmental variables from HMSC model
                q.env <- copy(i.model$XData)
                if ("." %in% names(q.env)) {
                        q.env <- q.env[, -which(names(q.env) == ".")]
                }
                #  Create a classification of sites based on the selected
                #  environmental variables
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
                
                
               
                setDT(q.env)
                q.env[, type := q.clusters]
                #  find the centroid for each variable that follows the typology
                q.centroids <- q.env[, lapply(.SD, mean), .SDcols = q.all.vars, by = "type"]
                setorderv(q.centroids, "type")
                q.centroids[, type := NULL]
                q.centroids <- as.matrix(as.data.frame(q.centroids))
                
                #  prepare environmental variables for adjustments in script below
                q.observations <- as.matrix(as.data.frame(q.cluster.env))
                
                # Based on two vectors of how strongly the centroids and sites are contracted
                # recompute environmental variables
                source("code/helper/new_predictors_for_hard_classification.R")
                i.success <- i.success + nrow(q.rating.out2)
                i.nc <- append(i.nc, rep(q.nc, nrow(q.rating.out2)))
                #print(i.success)
                # If we reach this part of the loop we have successfully extracted realistic artificial environments 
                # Thus it makes now sense to store results for later 
                # Add number to vector for later reference
                i.n.variables <- append(i.n.variables, rep(q.nvariables, nrow(q.rating.out2)))
                # Add importance of variables
                q.importance <- sum(i.vp$importance[q.all.vars])
                i.importance <- append(i.importance, rep(q.importance, nrow(q.rating.out2)))
                #  store cluster assignment for later
                for (foo in 1:nrow(q.rating.out2)) i.cluster.assignments[[length(i.cluster.assignments) + 1]] <- q.clusters
                rm(foo)
                # update the vector that determined probabilities 
                i.variable.counter[q.all.vars] <- i.variable.counter[q.all.vars] + 1 
                
                # fuzzy classification
                q.fuzzy <- lapply(1:length(q.newenv), function(x) {
                        vegclust(
                                q.newenv[[x]][, which(colnames(q.newenv[[x]]) %in% q.all.vars)],
                                mobileCenters = q.nc,
                                method = "FCM",
                                m = 1.5
                        )
                })
                for (foo in 1:nrow(q.rating.out2)) i.fuzzy.assignments[[length(i.fuzzy.assignments) + 1]] <- q.fuzzy
                rm(foo)
                # Judge quality of fuzzy classification with normalized partitioning entropy.
                # Analog to ASW for hard classification. 
                q.npe <- sapply(q.fuzzy, NPE)
                i.npe <- append(i.npe, q.npe)
                
                
                # Combine MCMC chains for prediction
                
                # select the last fifty samples from each chain 
                chain.length <- mcmcPar[scheme_id == model.name, nSamples]
                chain.length.low <- chain.length - 50
                chain.sampels <- c(chain.length.low:chain.length, (chain.length + chain.length.low):(2*chain.length), ((2*chain.length)+chain.length.low):(3*chain.length))
                
                q.posterior_samples <- poolMcmcChains(i.model$postList)
                q.selected_samples  <- q.posterior_samples[chain.sampels]
                
                #  — — — predict new biotic communities  —  —  —  —
                q.out <- lapply(
                        q.newenv,
                        FUN = function(pre) {
                                x <- predict(
                                        object = i.model,
                                        post = q.selected_samples,
                                        X = pre
                                )
                                # x <- x[[2]]
                                return(x)
                        }
                )
                
                i.out[[length(i.out) + 1]] <- q.out
                rm(list = ls()[grepl("^q\\.", x = ls())])
                
        }
        rm(q)
        i.out2 <- unlist(i.out, recursive = FALSE)
        #  if more than 200 reduce to 200, if less than two hundred thann all available 
        i.index <- min(length(i.out2), 200)
        i.index.range <- 1:i.index
        out_list <-
                list(
                        data = i.out2[i.index.range],
                        number_of_clusters = i.nc[i.index.range],
                        number_of_variables = i.n.variables[i.index.range],
                        contraction_points  = i.contraction.points[i.index.range],
                        contraction_centroids = i.contraction.centroids[i.index.range],
                        variable_importance = i.importance[i.index.range],
                        asw = i.asw[i.index.range],
                        npe = i.npe[i.index.range],
                        hard_cluster_assignment = i.cluster.assignments[i.index.range],
                        fuzzy_cluster_assignment = i.fuzzy.assignments[i.index.range],
                        info = i.ss,
                        converged = convergence_state
                )
        i.save.name <- sub(x = model.files[i],
                           pattern = "fitted_hmsc_models",
                           replacement = "simulated_data")
        saveRDS(out_list, i.save.name)
        rm(out_list)
        rm(list = ls()[grepl("^i\\.", x = ls())])
}
rm(i)


