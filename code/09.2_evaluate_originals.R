################################################################################
# Script Name:        evaluate_originals.R
# Description:        evaluate coherence between simulated typologies and communities 
#                        
#
# Author:             Jonathan Jupke 
# Date Created:       2025-09-15
# Last Modified:      2025-09-15
#
# R Version:          R 4.5.1
# Required Packages:  data.table, mvabund, parallel, vegan, parallelDist, doParallel, zetadiv, foreach, indicspecies, philentropy, labdsv
#
# Notes:              
################################################################################


# setup -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())

# libraries ---------------------------------------------------------------
library(data.table)
library(mvabund)
library(parallel)
library(vegan)
library(parallelDist)
library(doParallel)
library(zetadiv)
library(foreach)
library(indicspecies)
library(philentropy)
library(labdsv)

# custom scripts ----------------------------------------------------------

source("code/functions/prop_sample.R")
source("code/functions/calculate_auc.R")
source("code/functions/group_sites_by_cluster.R")


# script wide parameters --------------------------------------------------

cores_to_use <- detectCores() - 2

scheme_types    <- list.files("data/scheme_types/", full.names = T)
unfitted_models <- list.files("data/unfitted_hmsc_models/", full.names = T)
ss_files        <- list.files("data/spatial_scale/", full.names = T)
#vp.files        <- list.files("data/variation_partitioning/", full.names = T)

out <- list()

#for (i in 1:length(unfitted_models)) {
for (i in 20) {
     
        if (i == 1) first_test <- list()
        # update to console 
        print(paste(i, "___", Sys.time()))
        i.out <- list()
        #- Read the simulated data sets created in script 08_simulate_data.R
        #- Each of these files is a list with the following elements: 
        # data = a predicted community data set,
        # number_of_variables = the number of variables the typology is based on
        # contraction_points  = the factor by which samples are contracted in environmental space 
        # contraction_centroids= the factor by which cluster centroids are contracted in environmental space 
        # variable_importance = cumulative importance of variables on which the typology is based   
        # asw = average silhouette width of environmental clusters
        # hard_cluster_assignment = type assignment of each sample in hard clustering
        # fuzzy_cluster_assignment =  type assignment of each sample in fuzzy clustering
        
        i.file <- readRDS(unfitted_models[i])
        #- Identify and read the corresponding spatial scale file to add scheme ID to model
        i.ss.index <- sub("unfitted_hmsc_models", "spatial_scale", x = unfitted_models[i])
        i.ss.index <- which(ss_files == i.ss.index)
        if (length(i.ss.index) == 0) break()
        i.ss <- readRDS(ss_files[i.ss.index])
        #- same for typology
        i.ss.index <- sub("unfitted_hmsc_models", "scheme_types", x = unfitted_models[i])
        i.ss.index <- which(scheme_types == i.ss.index)
        if (length(i.ss.index) == 0) break()
        i.typ <- readRDS(scheme_types[i.ss.index])
        
        ## tests 
        if (nrow(i.typ) != i.ss$samples) stop("The rows in the typology data set do not match the scheme.")
        
        # Extract biotic communities 
        i.Y <- i.file$Y 
        if (i.ss$organismQuantityType == "presence"){
                i.D <- vegdist(i.Y, method = "jaccard")
        } else if (i.ss$organismQuantityType == "individuals"){
                i.D <- vegdist(i.Y, method = "bray")
        }
       
        # Extract environmental variables 
        i.X <- i.file$XData
        # remove MEMs
        if (any(grepl("MEM", names(i.X)))) i.X <- i.X[, -which(grepl("MEM", names(i.X)))]
        
        # k-means clustering 
        i.kmeans.asw <-  i.smallest.cluster <- i.n.small.clusters <- c()
        for (numberOfClusters in 2:10) {
                nc.clust <- kmeans(
                        x = i.X,
                        centers = numberOfClusters,
                        nstart = 100,
                        iter.max = 100
                )
                i.smallest.cluster[numberOfClusters - 1] <- min(nc.clust$size)
                i.n.small.clusters[numberOfClusters - 1] <- sum(nc.clust$size < 10)
                if (any(nc.clust$size < 10)){
                        rm(list = ls()[grepl("^nc\\.", x = ls())])
                        next()
                }
                
                nc.asw <- cluster::silhouette(nc.clust$cluster,
                                     dist(i.X))
                i.kmeans.asw[length(i.kmeans.asw) + 1] <- mean(nc.asw[, 3])
                names(i.kmeans.asw)[length(i.kmeans.asw)] <- paste(numberOfClusters)
                rm(list = ls()[grepl("^nc\\.", x = ls())])
        }
        # if there is no solution where the smallest cluster has at least 10 samples we need to balance manually
        if (is.null(i.kmeans.asw)) {
                # best available number
                # which.max only returns one even if multiple solutions share the same number
                names(i.smallest.cluster) <- 2:10
                i.smallest.cluster <- i.smallest.cluster[i.n.small.clusters == 1]
                i.nc <- as.numeric(names(i.smallest.cluster[max(which(
                        i.smallest.cluster == max(i.smallest.cluster)
                ))]))
                i.clusters <- kmeans(
                        x = i.X,
                        centers = i.nc,
                        nstart = 100,
                        iter.max = 100
                )
                i.clusters <- i.clusters$cluster
                i.clusters <- balance_clusters(
                        data = i.X,
                        clusters = i.clusters,
                        min_size = 10
                )
                
        } else {
                i.nc <- as.numeric(names(which.max(i.kmeans.asw)))
                i.clusters <- kmeans(
                        x = i.X,
                        centers = i.nc,
                        nstart = 100,
                        iter.max = 100
                )
                i.clusters <- i.clusters$cluster
        }
        
        # create a data.frame of types. Contains Kmeans result and all real typologies which have multiple types 

        # i.D.list creation -------------------------------------------------------
        i.D.list <- list()
        for (ii in 1:ncol(i.typ)){
                
                i.D.list[[ii]] <- list(distance = NA, types = NA, Y = NA)
                ii.typ  <- names(i.typ)[ii]
                ii.typ  <- i.typ[, ..ii.typ]
                ii.typ  <- unlist(ii.typ)
                ii.typ  <- as.character(ii.typ)
                ii.drop <- which(is.na(ii.typ))
                ii.tab  <- table(ii.typ)
                # only one type? No comparisons possible
                if (length(ii.tab) == 1){
                        next()
                }
                if (any(ii.tab < 5)){
                        ii.small.id <- which(ii.tab < 5)
                        ii.small.id <- names(ii.small.id)
                        ii.drop <- append(ii.drop, which(ii.typ %in% ii.small.id))
                        # is there more than one type left after dropping 
                        ii.type.dropped <- ii.typ[-ii.drop]
                        ii.tab <- table(ii.type.dropped)
                        if (length(ii.tab) == 1){
                                next()
                        }
                        ii.D <- as.matrix(i.D)
                        ii.D <- ii.D[-ii.drop, -ii.drop]
                        ii.D <- as.dist(ii.D)
                        i.D.list[[ii]]$distance <- ii.D
                        i.D.list[[ii]]$types <- ii.type.dropped
                        i.D.list[[ii]]$Y <- i.Y[-ii.drop, ]
                        
                        # tests 
                        if (attr(i.D.list[[ii]]$distance, "Size") != length(i.D.list[[ii]]$types)) stop("dimensions of distance and types do not match")
                        if (nrow(i.D.list[[ii]]$Y) != length(i.D.list[[ii]]$types)) stop ("dimensions of Y and types do not match")
                        
                } else {
                        if (length(ii.drop) > 0){
                                ii.typ <- ii.typ[-ii.drop]
                                ii.D <- as.matrix(i.D)
                                ii.D <- ii.D[-ii.drop, -ii.drop]
                                ii.D <- as.dist(ii.D)   
                                i.D.list[[ii]]$Y <- i.Y[-ii.drop, ]
                                i.D.list[[ii]]$distance <- ii.D
                        } else {
                                i.D.list[[ii]]$Y <- i.Y
                                i.D.list[[ii]]$distance <- i.D
                        }

                        i.D.list[[ii]]$types <- unlist(ii.typ)
                        
                        
                        # tests 
                        if (attr(i.D.list[[ii]]$distance, "Size") != length(i.D.list[[ii]]$types)) stop("dimensions of distance and types do not match")
                        if (nrow(i.D.list[[ii]]$Y) != length(i.D.list[[ii]]$types)) stop ("dimensions of Y and types do not match")
                        
                }
                rm(list = ls()[grepl(x = ls(), pattern = "^ii\\.")])
        };rm(ii);gc()
        
        # ADD kmeans to list 
        i.D.list[[length(i.D.list) + 1]] <- list(distance = NA, types = NA, Y = NA)
        i.D.list[[length(i.D.list)]]$distance <- i.D
        i.D.list[[length(i.D.list)]]$types     <- i.clusters
        i.D.list[[length(i.D.list)]]$Y         <- i.Y
        names(i.D.list) <- c(names(i.typ), "kmeans")
        
        # drop empty list elements, i.e., typologies with only one type. 
        i.id <- sapply(i.D.list, function(x) !all(is.na(x$distance)))
        i.D.list <- i.D.list[i.id]
        
        # evaluations -------------------------------------------------------------
        # Classification Strength -------------------------------------------------
        
        i.cs <- sapply(i.D.list, function(x) {
                y <- meandist(dist = x$distance, grouping = x$types)
                y <- unlist(summary(y))["CS"]
                y <- round(y, 3)
                y
        })
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.cs, types = names(i.D.list), metric = "cs")
        # ANOSIM -----------------------------------------------------------------
        print(paste(i, "___ ANOSIM START"))
        i.N = length(i.D.list)
        # Preallocate objects 
        i.out.r.min <- 
                i.out.r.men <- 
                i.out.r.max <-  numeric(length = i.N)
        
        i.combinations <- 
                lapply(
                        i.D.list, 
                        function(x) {
                                x$types |>
                                        unique() |>
                                        combn(m = 2) |>
                                        t()
                        }
                )  
        # Pre-compute cluster indices
        cluster_indices <- list()
        # Loop through each unique cluster set
        for (ci in 1:length(i.D.list)) {
                # Get the cluster assignment for this set
                cluster_assignments <- i.D.list[[ci]]$types
                # Pre-compute all possible pair indices for this cluster assignment
                # A pair index returns the indices of all samples that are in a given
                # Combination of types. 
                pair_indices <- list()
                for (pa in 1:nrow(i.combinations[[ci]])) {
                        pair_indices[[pa]] <- which(cluster_assignments %in% i.combinations[[ci]][pa, ])
                }
                cluster_indices[[ci]] <- pair_indices
        }
        rm(ci)
        # Setup parllel backend
        cl <- makeCluster(cores_to_use)
        registerDoParallel(cl)
        
        
        results <- foreach(evaluation = 1:i.N, .combine = "rbind", .packages = c('vegan', 'dplyr')) %dopar% {
       
                        # Get cluster set index
                combination_set <- i.combinations[[evaluation]]
                
                pa_results <- 
                        t(
                                sapply(
                                        1:nrow(combination_set), 
                                        function(pa) {
                                                pa_id <- cluster_indices[[evaluation]][[pa]]
                                                dist_matrix <- as.matrix(i.D.list[[evaluation]]$distance)
                                                pa_dist <- as.dist(dist_matrix[pa_id, pa_id])
                                                
                                                pa_ano <- anosim(
                                                        x = pa_dist,
                                                        grouping = i.D.list[[evaluation]]$types[pa_id],
                                                        permutations = 5
                                                )
                                                
                                                return(c(pa_ano$statistic, evaluation))
                                        }
                                )
                        )
        
                return(pa_results)
        }
        stopCluster(cl)
        for (evaluation in 1:i.N) {
                eval_results <- results[which(results[,2] == evaluation), 1:2]
                
                # Get min, mean, max for this evaluation
                if (is.null(nrow(eval_results))){
                        i.out.r.min[evaluation] <- 
                        i.out.r.men[evaluation] <- 
                        i.out.r.max[evaluation] <-eval_results[1]
                } else {
                        i.out.r.min[evaluation] <- min(eval_results[,1])
                        i.out.r.men[evaluation] <- mean(eval_results[,1])
                        i.out.r.max[evaluation] <- max(eval_results[,1])
                }
        rm(eval_results)
        };rm(evaluation)
        
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.out.r.min, types = names(i.D.list), metric = "ANOSIM R min")
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.out.r.men, types = names(i.D.list), metric = "ANOSIM R mean")
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.out.r.max, types = names(i.D.list), metric = "ANOSIM R max")
        
        rm(cluster_assignments)
        rm(cluster_indices)
        rm(results)
        rm(eval_results)
        
        # AUC zeta ----------------------------------------------------------------
        print(paste(i, "___ AUC ZETA START"))
        i.baseline <- vector(mode = "numeric", length = i.N)
        # baseline is the inter-type AUCζ that intra-type values can be compared to
        for (evaluation in 1:i.N){
                ev.types           <- i.D.list[[evaluation]]$types
                ev.baseline2       <- vector(mode = "numeric", length = 10)
                # 10 repetitions for each evaluation 
                for (evaluation2 in 1:10){
                        # compute the mean number of samples per type, rounded 
                        # to next integer. Rounding is done to ensure integers,
                        # which is necessary as prop_samples will end in an end-
                        # less while loop if N is not an integer. 
                        ev2.N  <- ev.types |> table() |> mean() |>  round()
                        #- draw ids to sample simulated data with. 
                        #- the number of samples from each type is proportional
                        #- to the total number is sites in each type.
                        ev2.id     <- prop_sample(
                                x = ev.types, 
                                N = ev2.N,
                                string = T
                        )
                        #- sample simulated data 
                        ev2.baseline.data <- i.D.list[[evaluation]]$Y[ev2.id, ]
                        #- compute AUCζ baseline for simulated data
                        ev.baseline2[evaluation2] <- 
                                Zeta.decline.ex(
                                        data.spec  = ev2.baseline.data,
                                        orders     = 1:10,
                                        plot       = FALSE,
                                        rescale    = TRUE)$zeta.val |> 
                                calculate_auc(x = 1:10)    
                        rm(list = ls()[grepl(pattern = "^ev2\\.", x = ls())])
                }
                i.baseline[evaluation] <- mean(ev.baseline2)
                rm(list = ls()[grepl(pattern = "^ev\\.", x = ls())])
        }
        # turn vector into data.frame with additional ID column
        i.baseline        <- data.frame(i.baseline, 1:i.N)
        names(i.baseline) <- c("baseline", "ID")

        # apply the custom function group_sites_by_cluster to create a list 
        # for each type that shows which sites are part of that type
        for (ii in 1:i.N) {
                if (ii == 1) i.cluster_indices <- list()
                i.cluster_indices[[ii]] <- group_sites_by_cluster(i.D.list[[ii]]$types)
        }; rm(ii)
        i.zeta <- list()
        for (ii in 1:i.N){
                ii.out <- list()
                ii.combination_set <- i.cluster_indices[[ii]]        
                ii.orders          <- sapply(ii.combination_set, function(x) min(length(x), 10))
                for (iii in 1:length(ii.combination_set)){
                        pa.id <- ii.combination_set[[iii]]
                        pa.zeta  <- Zeta.decline.ex(
                                data.spec = i.D.list[[evaluation]]$Y[pa.id,],
                                orders     = 1:ii.orders[iii],
                                plot       = FALSE,
                                rescale    = TRUE)
                        pa.zeta <- pa.zeta$zeta.val
                        pa.zeta <- calculate_auc(pa.zeta, x = 1:ii.orders[iii])
                        pa.zeta <- data.table(raw_zeta = pa.zeta, ID = ii)
                        ii.out[[length(ii.out) + 1]] <- pa.zeta
                        rm(list = ls()[grepl(x = ls(), pattern = "^pa\\.")])
                }
                i.zeta[[length(i.zeta) + 1]] <- rbindlist(ii.out)
                rm(list = ls()[grepl(x = ls(), pattern = "^ii\\.")])
        }
        i.zeta <- rbindlist(i.zeta)
        setDT(i.baseline)
        i.zeta2 <- i.baseline[i.zeta, on = "ID"]
        i.zeta2[, auczeta_relative := raw_zeta - baseline]
        i.zeta2[, min  :=  min(auczeta_relative), by = "ID"]
        i.zeta2[, mean := mean(auczeta_relative), by = "ID"]
        i.zeta2[, max  :=  max(auczeta_relative), by = "ID"]
        i.zeta2 <- unique(i.zeta2, by = "ID")
        
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.zeta2$min, types  = names(i.D.list), metric = "AUCζ min" )
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.zeta2$mean, types = names(i.D.list), metric = "AUCζ mean")
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.zeta2$max, types  = names(i.D.list), metric = "AUCζ max" )

# PERMANOVA ---------------------------------------------------------------
        print(paste(i, "___ PERMANOVA START"))
        i.r <- i.F <- 
                vector(mode = "numeric", length = i.N)
        for (ii in 1:i.N){
                ii.types <- i.D.list[[ii]]$types	
                ii.perma <- adonis2(formula = i.D.list[[ii]]$distance ~ ii.types, permutations = 1)
                i.r[ii]  <- ii.perma$R2[1]
                i.F[ii]  <- ii.perma$F[1]
                rm(ii.types, ii.perma)
        };rm(ii)
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.r, types = names(i.D.list), metric = "PERMANOVA R2")
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = i.F, types = names(i.D.list), metric = "PERMANOVA F")

# Indicator Value ---------------------------------------------------------
        # print(paste(i, "___ INDVAL START"))
        #  # Setup parallel backend
        # cl <- makeCluster(cores_to_use)
        # registerDoParallel(cl)
        # # Start parallel loop 
        # results <- foreach(evaluation = 1:i.N, .combine = "rbind", .packages = c('indicspecies', 'permute')) %dopar% {
        #         
        #         
        #         isa <- multipatt(
        #                 x = i.D.list[[evaluation]]$Y, 
        #                          i.D.list[[evaluation]]$types, 
        #                          func = "indval.g",
        #                          control = permute::how(nperm = 500))
        #         isa$sign$holm_p_value <- p.adjust(isa$sign$p.value, method = "holm")
        #         out1 <- nrow(isa$sign[which(isa$sign$holm_p_value<=0.05), ])/nrow(isa$sign)
        #         out2 <- mean(isa$sign$holm_p_value, na.rm = T)
        #         cbind(out1, out2)
        # }
        # #rm(isa, out1, out2, evaluation)
        # i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = results[,1], types = names(i.D.list), metric = "isa_number")
        # i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = results[,2], types = names(i.D.list), metric = "isa_avg_p")

# ISAMIC ------------------------------------------------------------------
        print(paste(i, "___ ISAMIC START"))
        j.isamic <- vector(mode = "numeric", length = i.N)
        for (evaluation in 1:i.N){
                isa <- isamic(
                        comm = i.D.list[[evaluation]]$Y, 
                        clustering = i.D.list[[evaluation]]$types)
                j.isamic[evaluation] <- mean(isa)
        }
        rm(isa, evaluation)
        i.out[[length(i.out)+1]] <- data.table(scheme = i.ss$scheme_id, value = j.isamic, types = names(i.D.list), metric = "isamic")
        rm(j.isamic)
        i.out <- rbindlist(i.out)
        out[[i]] <- i.out
        rm(list = ls()[grepl(x = ls(), pattern = "^i\\.")])
}
        
