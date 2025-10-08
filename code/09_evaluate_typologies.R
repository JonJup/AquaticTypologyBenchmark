################################################################################
# Script Name:        evaluate_typolohies.R
# Description:        evaluate coherence between simulated typologies and communities 
#                        
#
# Author:             Jonathan Jupke 
# Date Created:       2025-09-15
# Last Modified:      2025-09-22
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
#library(indicspecies)
library(philentropy)
library(labdsv)

# custom scripts ----------------------------------------------------------
source("code/functions/render_table.R")
source("code/functions/prop_sample.R")
source("code/functions/calculate_auc.R")
source("code/functions/group_sites_by_cluster.R")
source("code/functions/cv_d_squared.R")

# script wide parameters --------------------------------------------------
source("code/parameters/run_parameters.R")
cores_to_use <- detectCores() - 2

simulated_files <- list.files("data/simulated_data/", full.names = T)
ss_files        <- list.files("data/spatial_scale/", full.names = T)
vp.files        <- list.files("data/variation_partitioning/", full.names = T)

for (i in 1:length(simulated_files)) {
        if (i != 1) next()
        # update to console 
        paste(i, "___", Sys.time())
        
        # list to store all the results of this loop
        out <- list()
        
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

        i.file <- readRDS(simulated_files[i])
        i.N <- length(i.file$data)
        #i.file <- unlist(i.file, recursive = F)
        #- Identify and read the corresponding spatial scale file to add scheme ID to model
        i.ss.index <- sub("simulated_data", "spatial_scale", x = simulated_files[i])
        i.ss.index <- which(ss_files == i.ss.index)
        if (length(i.ss.index) == 0) break()
        i.ss <- readRDS(ss_files[i.ss.index])
        
        #- compute distance matrices
        if (i.ss$organismQuantityType == "presence"){
                i.distance.matrices <-
                        rapply(
                                i.file$data,
                                parallelDist,
                                method = "binary",
                                how = "list"
                        ) 
        } else if (i.ss$organismQuantityType == "individuals"){
                i.distance.matrices <-
                        rapply(
                                i.file$data,
                                parallelDist,
                                method = "bray",
                                how = "list"
                        ) 
        } else {
                stop("unknown organism Quantity type")
        }
        
      


        # classification strength ------------------------------------------------
        i.cs <- sapply(
                1:length(i.distance.matrices),
                        function(x) {
                                y <- meandist(
                                        dist     = i.distance.matrices[[x]],
                                        grouping = i.file$hard_cluster_assignment[[x]]
                                )
                        unlist(summary(y))["CS"]
                        }
                )
        i.cs <- render_table(i.cs, "classification strength")
        out[[length(out) + 1]] <- i.cs
        rm(i.cs)        
        # ANOSIM -----------------------------------------------------------------

        # Preallocate objects 
        i.out.r.min <- 
        i.out.r.men <- 
        i.out.r.max <-  numeric(length = i.N)
        
        i.combinations <- 
                lapply(
                        1:i.N, 
                        function(x) {
                                i.file$hard_cluster_assignment[[x]] |>
                                unique() |>
                                combn(m = 2) |>
                                t()
                        }
                )  
        # Pre-compute cluster indices
        cluster_indices <- list()
        # Loop through each unique cluster set
        for (ci in 1:i.N) {
                # Get the cluster assignment for this set
                cluster_assignments <- i.file$hard_cluster_assignment[[ci]]
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
                cluster_set <- evaluation
                combination_set <- i.combinations[[cluster_set]]
               
                pa_results <- 
                        t(
                                sapply(
                                        1:nrow(combination_set), 
                                        function(pa) {
                                                pa_id <- cluster_indices[[cluster_set]][[pa]]
                                                dist_matrix <- as.matrix(i.distance.matrices[[evaluation]])
                                                pa_dist <- as.dist(dist_matrix[pa_id, pa_id])

                                                pa_ano <- anosim(
                                                        x = pa_dist,
                                                        grouping = i.file$hard_cluster_assignment[[cluster_set]][pa_id],
                                                        permutations = 5
                                                )

                                                return(c(pa_ano$statistic, pa_ano$signif, evaluation))
                                        }
                                 )
                        )
                return(pa_results)
        }
        stopCluster(cl)
        for (evaluation in 1:i.N) {
                eval_results <- results[which(results[,3] == evaluation), 1:2]
                
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
                
                # i.out.r.min[evaluation] <- min(eval_results[1])
                # i.out.r.men[evaluation] <- mean(eval_results[1])
                # i.out.r.max[evaluation] <- max(eval_results[1])
                
                # i.out.p.min[evaluation] <- min(eval_results[, 2])
                # i.out.p.men[evaluation] <- mean(eval_results[, 2])
                # i.out.p.max[evaluation] <- max(eval_results[, 2])
        }
        # Move into final format 
        out[[length(out) + 1]] <- render_table(i.out.r.min, "ANSOIM R min")
        out[[length(out) + 1]] <- render_table(i.out.r.men, "ANSOIM R mean")
        out[[length(out) + 1]] <- render_table(i.out.r.max, "ANSOIM R max")
        # out[[length(out) + 1]] <- render_table(i.out.p.min, "ANSOIM p min")
        # out[[length(out) + 1]] <- render_table(i.out.p.men, "ANSOIM p mean")
        # out[[length(out) + 1]] <- render_table(i.out.p.max, "ANSOIM p max")

        rm(cluster_assignments)
        rm(cluster_indices)
        rm(results)
        rm(eval_results)


# AUC Zeta ----------------------------------------------------------------
        # The AUC baseline is the inter-type AUCζ that intra-type values can be compared to
        baseline <-numeric(i.N)
        # Pre-calculate constants outside loops where possible
        zeta_orders <- 1:10
        zeta_x <- 1:10
        n_reps <- 10
        for (evaluation in seq_len(i.N)) {
                # Get cluster assignments
                ev_types <- i.file$hard_cluster_assignment[[evaluation]]
                
                # Pre-calculate the mean sample size (avoid recalculating)
                type_counts <- table(ev_types)
                ev_N <- round(mean(type_counts))
                
                # Pre-allocate for inner loop results
                ev_baseline2 <- numeric(n_reps)
                
                # Get the data this evaluation
                eval_data <- i.file$data[[evaluation]]
                
                for (evaluation2 in seq_len(n_reps)) {
                        # Sample indices
                        ev2_id <- prop_sample(x = ev_types, N = ev_N)
                        
                        # Sample data and compute AUCζ
                        ev_baseline2[evaluation2] <- 
                                eval_data[ev2_id, ] |>
                                Zeta.decline.ex(
                                        orders = zeta_orders,
                                        plot = FALSE,
                                        rescale = TRUE
                                ) |>
                                _$zeta.val |>
                                calculate_auc(x = zeta_x)
                }
                baseline[evaluation] <- mean(ev_baseline2)
        }
        
        # for (evaluation in 1:i.N){
        #         #print(paste(evaluation, "/", i.N))
        #        
        #         ev.cluster_set     <- evaluation
        #         ev.types           <- i.file$hard_cluster_assignment[[ev.cluster_set]]
        #         ev.baseline2       <- vector(mode = "numeric", length = 10)
        #         # 10 repetitions for each evaluation 
        #         for (evaluation2 in 1:10){
        #                 # compute the mean number of samples per type, rounded 
        #                 # to next integer. Rounding is done to ensure integers,
        #                 # which is necessary as prop_samples will end in an end-
        #                 # less while loop if N is not an integer. 
        #                 ev2.N  <- ev.types |> table() |> mean() |> round()
        #                 #- draw ids to sample simulated data with. 
        #                 #- the number of samples from each type is proportional
        #                 #- to the total number is sites in each type.
        #                 ev2.id     <- prop_sample(
        #                         x = ev.types, 
        #                         N = ev2.N
        #                         )
        #                 #- sample simulated data 
        #                 ev2.baseline.data <- i.file$data[[evaluation]][ev2.id, ]
        #                 #- compute AUCζ baseline for simulated data
        #                 ev.baseline2[evaluation2] <- 
        #                         Zeta.decline.ex(
        #                                 data.spec  = ev2.baseline.data,
        #                                 orders     = 1:10,
        #                                 plot       = FALSE,
        #                                 rescale    = TRUE)$zeta.val |> 
        #                         calculate_auc(x = 1:10)    
        #                 rm(list = ls()[grepl(pattern = "^ev2\\.", x = ls())])
        #         }
        #         baseline[evaluation] <- mean(ev.baseline2)
        #         rm(list = ls()[grepl(pattern = "^ev\\.", x = ls())])
        # }
        # turn vector into data.frame with additional ID column
        baseline        <- data.frame(baseline, 1:i.N)
        names(baseline) <- c("baseline", "ID")
        # pre-allocate object 
        cluster_indices <- list()
        # apply the custom function group_sites_by_cluster to create a list 
        # for each type that shows which sites are part of that type
        for (evaluation in seq_along(i.file$hard_cluster_assignment)) {
                cluster_indices[[evaluation]] <- group_sites_by_cluster(i.file$hard_cluster_assignment[[evaluation]])
        }
        # Setup parallel backend
        cl <- makeCluster(cores_to_use)
        registerDoParallel(cl)
        results <- foreach(evaluation = 1:i.N, .combine = "rbind", .packages = c('zetadiv', 'dplyr')) %dopar% {

                cluster_set     <- evaluation
                combination_set <- cluster_indices[[cluster_set]]        
                orders = sapply(combination_set, function(x) min(length(x), 10))
                zetaout <- sapply(
                        #1:length(combination_set), 
                        2, 
                        function(pa) {
                                pa_id <- combination_set[[pa]]
                                if (sum(i.file$data[[evaluation]][pa_id,]) == 0){
                                        zeta = 0
                                        zeta
                                } else {
                                        zeta <- 
                                                Zeta.decline.ex(
                                                        data.spec = i.file$data[[evaluation]][pa_id,],
                                                        orders     = 1:orders[pa],
                                                        plot       = FALSE,
                                                        rescale    = TRUE)$zeta.val |> 
                                                calculate_auc(x = 1:orders[pa])
        
                                        zeta
                                }
                        }
                 )
                 matrix(c(zetaout, rep(evaluation, times = length(zetaout))), ncol = 2)  
        }
        stopCluster(cl)
        results2 <- as.data.frame(results)
        names(results2)[which(names(results2) == "V1")] <- "auczeta_raw"
        names(results2)[which(names(results2) == "V2")] <- "ID"
        setDT(results2)
        setDT(baseline)
        results2 <- baseline[results2, on = "ID"]
        results2[, auczeta_relative := auczeta_raw - baseline]
        results2[, min :=  min(auczeta_relative), by = "ID"]
        results2[, mean := mean(auczeta_relative), by = "ID"]
        results2[, max :=  max(auczeta_relative), by = "ID"]
        out[[length(out) + 1]] <- render_table(results2$min, "AucZeta min")        
        out[[length(out) + 1]] <- render_table(results2$mean, "AucZeta mean")        
        out[[length(out) + 1]] <- render_table(results2$max, "AucZeta max")        
        rm(baseline)
        #rm(baseline2)

        #——— Compute PERMANOVA ———————————————————————————————————————————————————
        i.r <- i.p  <- i.F  <- 
                vector(mode = "numeric", length = i.N)
        for (evaluation in 1:i.N){
                types <- i.file$hard_cluster_assignment[[evaluation]]	
                perma <- adonis2(formula = i.distance.matrices[[evaluation]] ~ types, permutations = 10)
                i.r[evaluation]  <- perma$R2[1]
                i.F[evaluation]  <- perma$F[1]
                #i.p[evaluation]  <- perma$`Pr(>F)`[1]
        }
        rm(perma)
        rm(evaluation)
        out[[length(out) + 1]] <- render_table(i.r, "PERMANOVA R2")
        #out[[length(out) + 1]] <- render_table(i.p, "PERMANOVA p")
        out[[length(out) + 1]] <- render_table(i.F, "PERMANOVA F")
        
        # #———— COMPUTE Silhouette Width ———————————————————————————————————————————
        # j.asw <- vector(mode = "numeric", length = neval)
        # for (asw in 1:neval){
        #         out <- silhouette(dist = j.distance.matrix[[asw]], 
        #                           x = j.cluster.assignments[[ceiling(asw/(neval/max.q))]])
        #         j.asw[asw] <- mean(out[, "sil_width"])
        # }
        # j.asw %<>% render_table("asw")

        #———— Indicator Value ————————————————————————————————————————————
        # # Setup parallel backend
        # cl <- makeCluster(cores_to_use)
        # registerDoParallel(cl)
        # # Prepare vectors for 
        # #j.isa1 <- j.isa2 <- vector(mode = "numeric", length = i.N)
        # # Start parallel loop 
        # results <- foreach(evaluation = 1:i.N, .combine = "rbind", .packages = c('indicspecies', 'permute')) %dopar% {
        # 
        #         
        #         isa <- multipatt(i.file$data[[evaluation]], 
        #                          i.file$hard_cluster_assignment[[evaluation]], 
        #                          func = "indval.g",
        #                          control = permute::how(nperm = 500))
        #         isa$sign$holm_p_value <- p.adjust(isa$sign$p.value, method = "holm")
        #         out1 <- nrow(isa$sign[which(isa$sign$holm_p_value<=0.05), ])/nrow(isa$sign)
        #         out2 <- mean(isa$sign$holm_p_value, na.rm = T)
        #         cbind(out1, out2)
        # }
        # #rm(isa, out1, out2, evaluation)
        # 
        # out[[length(out) + 1]] <- render_table(results[,1], "isa_number")
        # out[[length(out) + 1]] <- render_table(results[,2], "isa_avg_p")
        # rm(results)
        #——— COMPUTE ISAMIC ——————————————————————————————————————————————————————
        j.isamic <- vector(mode = "numeric", length = i.N)
        for (evaluation in 1:i.N){
                isa <- isamic(
                        comm = i.file$data[[evaluation]], 
                        clustering = i.file$hard_cluster_assignment[[ceiling(evaluation/within.q)]])
                j.isamic[evaluation] <- mean(isa)
        }
        rm(isa, evaluation)
        out[[length(out) + 1]] <- render_table(j.isamic, "isamic")
        rm(j.isamic)

        # fuzzy classification ---------------------------------------------------
        i.fuzzy.out <- numeric(i.N)
        i.fuzzy.out2 <- numeric(i.N)
        for (evaluation in 1:i.N){
               
                i.fuzzy.out[evaluation] <- 
                        cv_d_squared(
                                species_data = i.file$data[[evaluation]], 
                                fuzzy_memberships = as.matrix(unlist(i.file$fuzzy_cluster_assignment, recursive = F)[[evaluation]]$memb)
                                     )
                i.fuzzy.out2[evaluation] <- mantel(
                                xdis = vegdist(i.file$data[[evaluation]]), 
                                ydis = dist(as.matrix(unlist(i.file$fuzzy_cluster_assignment, recursive = F)[[evaluation]]$memb)),
                                permutations = 1
                        )[["statistic"]]
                        
                        # test <fuzzy_mantel_test(
                        # species_data = i.file$data[[evaluation]], 
                        # fuzzy_memberships = as.matrix(unlist(i.file$fuzzy_cluster_assignment, recursive = F)[[evaluation]]$memb)
        }

        out[[length(out) + 1]] <- render_table(i.fuzzy.out, "fuzzy_divergence")
        out[[length(out) + 1]] <- render_table(i.fuzzy.out2, "fuzzy_mantel")
        out2 <- rbindlist(out, fill =T )
        i.save.name <- sub(pattern = "simulated_data/", replacement = "evaluations/", x = simulated_files[i])
        saveRDS(out2, i.save.name)
        rm(out)
        rm(out2)
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])

}

