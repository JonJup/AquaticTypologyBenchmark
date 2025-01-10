# ——————————————————————————————————————————————————— #
# ——— Compute dirichlet distribution of Paritions ——— # 
# ——————————————————————————————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 12.12.24

# LOAD DATA -------------------------------------------------------------------------
f <- dir_ls("data/fitted_models/")
d <- lapply(f, readRDS)
d <- lapply(d, function(x)x[[4]])
d <- lapply(d, function(x) x[, c("taxon", "driver", "scaled_values", "run")])
d <- lapply(d, pivot_wider, id_cols = c("taxon", "run"), names_from = driver, values_from = scaled_values)
d <- lapply(d, setDT)
d <- rbindlist(d)
d <- d[, lapply(.SD, mean), by = "run", .SDcols = c("env", "bio", "space", "stochastic")]
 
 
# d     <- lapply(d, function(x) summarise(x, env = mean(env), bio = mean(bio), space = mean(space), stochastic = mean(stochastic)))
# 
# 
# env <- d[driver == "env", scaled_values] 
# bio <- d[driver == "bio", scaled_values]
# spa <- d[driver == "space", scaled_values]
# sto <- d[driver == "stochastic", scaled_values]
# 
# d2 <- data.frame(env, bio, spa, sto)
# dir_data <- DR_data(d2)
# fit <- DirichReg(dir_data ~ 1)
# summary(fit)

normalize <- function(x) {
        return((x - min(x)) / (max(x) - min(x)))
}

d3 <- copy(d)
d3[, c("run") := NULL]
d3 <- as.data.frame(d3)
d3 <- as.matrix(d3) 
fit <- dirichlet.mle(d3)
dens <- ddirichlet(x=d3, alpha = fit$alpha)
dens <- normalize(dens)
out <- data.table(run = d$run, 
                  normalized_dirichlet_density = dens)

# save to file ----------------------------------------------------------------------
saveRDS(out, "data/dirichlet_weights.rds")
saveRDS(d, "data/VP_results_polished.rds")
rm(d, d3, f, fit, dens, out, normalize)
