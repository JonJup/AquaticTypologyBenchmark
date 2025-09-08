#- combine Dirichlet density with typology evaluations and spatial scales

# setup -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())
library(data.table)
#source("code/functions/quintile_score.R")

# load data ---------------------------------------------------------------
#dirichlet   <- readRDS("data/dirichlet_weights.rds")
VP          <- readRDS("data/VP_results_polished.rds")
TX          <- rbindlist(lapply(list.files("data/taxonomic_resolution/", full.names = T), readRDS))
SS          <- rbindlist(lapply(list.files("data/spatial_scale/", full.names = T), readRDS))
evaluations <- list.files("data/evaluations/", full.names = TRUE) |> lapply(readRDS)
evaluations <- rbindlist(evaluations)
spatial.scales <- list.files("data/spatial_scale/", full.names = TRUE)
spatial.scales <- rbindlist(lapply(spatial.scales, readRDS))

# TODO remove! This is supposed to be a temporary work around 
evaluations[, model2 := ceiling(model/10)]
evaluations[, model3 := model + 10 - model2*10]
evaluations[, model := paste(model2, model3, sep = "_")]
evaluations[, c("model2", "model3") := NULL]
VP[, model2 := ceiling(run/10)]
VP[, model3 := run + 10 - model2*10]
VP[, model := paste(model2, model3, sep = "_")]
VP[, c("model2", "model3") := NULL]
### END OF TEMPORARY work around 

joined <- merge(evaluations, dirichlet, by = "model")
joined <- merge(joined, VP, by = "model")
#- summarize results from different samples of the same run 
#joined[,value := mean(value), by = c("run", "types", "metric", "simulation", "q", "k", "s")]
#joined <- unique(joined, by = c("run", "types", "metric", "simulation", "q", "k", "s"))

#- add overall ranking of typology system into a five factor design
joined[, paste0("rank_", c("contraction_points", "contraction_centroids", "env_asw", "importance")) := 
               lapply(.SD, quintile_score), 
       .SDcols = c("contraction_points", "contraction_centroids", "env_asw", "importance")]

joined[, rank_sum := rowSums(.SD), .SDcols = patterns("^rank_")]
joined[, rank_sum2 := quintile_score(rank_sum)]



# save to file ------------------------------------------------------------
saveRDS(joined, "data/evaluations_w_dirichlet.rds")
