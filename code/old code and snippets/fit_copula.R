library(copula)
library(data.table)
library(fitdistrplus)

# load data -------------------------------------------------------------------------
d <- readRDS("data/prelim_vp_results.rds")
d

fit1 <- fitdist(d[driver == "bio"]$scaled_value+0.0000001, "beta")
fit2 <- fitdist(d[driver == "env"]$scaled_value+0.0000001, "beta")
fit3 <- fitdist(d[driver == "space"]$scaled_value+0.0000001, "beta")
fit4 <- fitdist(d[driver == "stochastic"]$scaled_value-0.000001, "beta")

plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)

#- check for alternative distributions or 1 or 0 inflated betas 

library(copula)

#- add goodnes of fit tests
x <- frankCopula(param = 1, dim = 4)

x2 <- mvdc(copula = x, 
             margins = rep("beta", 4), 
             paramMargins = list(
                  list(shape1 = unname(fit1$estimate[1]), shape2=unname(fit1$estimate[2])),   
                  list(shape1 = unname(fit2$estimate[1]), shape2=unname(fit2$estimate[2])),   
                  list(shape1 = unname(fit3$estimate[1]), shape2=unname(fit3$estimate[2])),   
                  list(shape1 = unname(fit4$estimate[1]), shape2=unname(fit4$estimate[2]))  
             ))

pMvdc(x = c(0.1,0.1,0.1,0.7), mvdc = x2)

library(tidyr)
d2 <- pivot_wider(d, id_cols = c("taxon", "run"), names_from = driver, values_from = scaled_values)
setDT(d2)
d2[, env := mean(env), by = "run"]
d2[, space :=      mean(space), by = "run"]
d2[, bio :=        mean(bio), by = "run"]
d2[, stochastic := mean(stochastic), by = "run"]
d2 <- unique(d2, by = "run")
d2[, total := env + space + bio + stochastic]
d2[, `:=` (env = env/total, space = space/total, bio = bio/total, stochastic = stochastic/bio)]

d2[, prop := pMvdc(c(bio, env, space, stochastic), mvdc = x2), by = "run"]
