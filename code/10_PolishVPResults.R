# ——————————————————————————————————————————————————— #
# ——— Compute dirichlet distribution of Paritions ——— # 
# ——————————————————————————————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 12.12.24

#library(groundhog)
# library(pacman)
# pkgs <- c(
#         "fs", 
#         "tidyr",
#         "data.table",
#         "DirichletReg",
#         "sirt"
# )
# lapply(pkgs, library, character.only = TRUE)
# #groundhog.library(pkgs, "2024-12-01")
# rm(pkgs)
library(data.table)
setwd(rstudioapi::getActiveProject())
# LOAD DATA -------------------------------------------------------------------------
f <- list.files("data/variation_partitioning/", full.names = T)
d <- lapply(f, readRDS)
TX  <- rbindlist(lapply(list.files("data/taxonomic_resolution/", full.names = T), readRDS))
SS  <- rbindlist(lapply(list.files("data/spatial_scale/", full.names = T), readRDS))
EVAL <- rbindlist(lapply(list.files("data/evaluations/", full.names = T), readRDS), fill = T)

# prepare data ------------------------------------------------------------
TXSS <- TX[SS, on = c("scheme_id", "taxon", "data.set", "year", "samples")]
d <- lapply(d, function(x) x$VP[, c("taxon", "driver", "scaled_values", "run")])
d <- lapply(d, function(df) {
        dt <- as.data.table(df)
        dcast(dt, taxon + run ~ driver, value.var = "scaled_values")
})
d <- rbindlist(d, use.names = TRUE)
d <- d[, lapply(.SD, mean), by = "run", .SDcols = c("env", "bio", "space", "stochastic")]
 
d <- cbind(d, TXSS) 

d2 <- d[EVAL, on = "scheme_id"]
saveRDS(d, "data/VP_results_polished.rds")
rm(d, f, SS, TX, TXSS, d2, EVAL)
