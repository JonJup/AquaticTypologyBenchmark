#- combine dirichlet density with typology evaluations

dirichlet   <- readRDS("data/dirichlet_weights.rds")
evaluations <- readRDS("data/241212_first_results.rds")
VP          <- readRDS("data/VP_results_polished.rds")

joined <- merge(evaluations, dirichlet, by = "run")
joined <- merge(joined, VP, by = "run")

#- summarize results from diffent samples of the same run 
joined[,value := mean(value), by = c("run", "types", "statistic")]
joined <- unique(joined, by = c("run", "types", "statistic"))

saveRDS(joined, "data/eval_w_dirichlet.rds")
