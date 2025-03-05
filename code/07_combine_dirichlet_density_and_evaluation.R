#- combine dirichlet density with typology evaluations

dirichlet   <- readRDS("data/dirichlet_weights.rds")
VP          <- readRDS("data/VP_results_polished.rds")
evaluations <- dir_ls("data/evaluations/") %>% lapply(readRDS)
evaluations <- unlist(evaluations, recursive = FALSE)
for (i in 1:length(evaluations)) evaluations[[i]]$run = i
evaluations <- rbindlist(evaluations)

joined <- merge(evaluations, dirichlet, by = "run")
joined <- merge(joined, VP, by = "run")

#- summarize results from different samples of the same run 
joined[,value := mean(value), by = c("run", "types", "metric", "simulation", "q", "k", "s")]
joined <- unique(joined, by = c("run", "types", "metric", "simulation", "q", "k", "s"))

#- add overall ranking of typology system into a five factor design
joined[, paste0("rank_", c("contraction_points", "contraction_centroids", "env_asw", "importance")) := 
               lapply(.SD, quintile_score), 
       .SDcols = c("contraction_points", "contraction_centroids", "env_asw", "importance")]

joined[, rank_sum := rowSums(.SD), .SDcols = patterns("^rank_")]
joined[, rank_sum2 := quintile_score(rank_sum)]



saveRDS(joined, "data/eval_w_dirichlet.rds")
