### - Variation Partitioning - ###

setwd(rstudioapi::getActiveProject())
library(Hmsc)
library(data.table)
models <- list.files("data/fitted_hmsc_models/", full.names = T)

for (i in seq_along(models)){
        #if (i < 8) next()
        print(i)
        i.models    <- readRDS(models[i])
        i.preds     <- computePredictedValues(i.models)
        i.MF        <- evaluateModelFit(hM = i.models, 
                                        predY = i.preds)
        # how many spatial predictors are there?
        i.n.spatial <- sum(grepl(x = colnames(i.models$XData), pattern = "AEM|MEM"))
        
        if (i.n.spatial == 0){
                i.VP    <- computeVariancePartitioning(
                        i.models,
                        group = c(rep(
                                1, ncol(i.models$XData) - i.n.spatial
                        )),
                        groupnames = c("env")
                )  
                i.VP2   <- data.table(
                        taxon = rep(colnames(i.VP$vals), each = 2),
                        driver = rep(c("env",  "bio"), times = ncol(i.VP$vals)),
                        value  = c(i.VP$vals),
                        r2     = rep(i.MF$TjurR2, each = 2)
                )
                i.VP22 <- data.table(
                        taxon = colnames(i.VP$vals),
                        driver = "space",
                        value  = 0,
                        r2     = i.MF$TjurR2
                )
                i.VP2 <- rbindlist(list(i.VP2, i.VP22))
        } else {
                i.VP    <- computeVariancePartitioning(
                        i.models,
                        group = c(rep(
                                1, ncol(i.models$XData) - i.n.spatial
                        ), rep(2, i.n.spatial)),
                        groupnames = c("env", "space")
                )  
                i.VP2   <- data.table(
                        taxon = rep(colnames(i.VP$vals), each = 3),
                        driver = rep(c("env", "space", "bio"), times = ncol(i.VP$vals)),
                        value  = c(i.VP$vals),
                        r2     = rep(i.MF$TjurR2, each = 3)
                )
        }
        
        i.VP2[r2 < 0, r2 := 0]
        i.VP2[, scaled_values := value * r2]
        #- add stochasticity = 1-r2
        i.VP3 <- data.table(
                taxon = colnames(i.VP$vals),
                driver = "stochastic",
                value = 0,
                r2 = i.MF$TjurR2
        )
        i.VP3[r2 < 0, r2 := 0]
        i.VP3[, value := 1 - r2]
        i.VP3[, scaled_values := 1 - r2]
        i.VP4 <- rbindlist(list(i.VP2, i.VP3))
        i.VP4$run <- i
        #i.VP4$psrf_mean <- rep(i.means, each = 4)
        #i.VP4$psrf_max  <- rep(i.maxes, each = 4)
        
        #- new VP to determine relative predictor importance
        i.VP5    <- computeVariancePartitioning(i.models)
        
        i.VP6 <- rowSums(i.VP5$vals)
        #- drop morans eigenvectors
        i.VP6 <- i.VP6[!grepl(x = names(i.VP6), pattern = "MEM|AEM|Random")]
        i.VP6 <- i.VP6 / sum(i.VP6)
        
        i.save.name <- sub(x = models[i], pattern = "fitted_hmsc_models", replacement = "variation_partitioning")
        saveRDS(list(VP = i.VP4, importance = i.VP6), i.save.name)
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
}


# old code  ---------------------------------------------------------------

# 
# # -- for environmental vs space
# o.preds     <- computePredictedValues(o.m1.mcmc)
# o.MF        <- evaluateModelFit(hM = o.m1.mcmc, predY = o.preds)
# # how many spatial predictors are there?
# o.n.spatial <- sum(str_detect(names(o.env.samples), "AEM|MEM"))
# o.VP    <- computeVariancePartitioning(
#         o.m1.mcmc,
#         group = c(rep(
#                 1, length(o.predictors) - o.n.spatial
#         ), rep(2, o.n.spatial)),
#         groupnames = c("env", "space")
# )
# o.VP2   <- data.table(
#         taxon = rep(colnames(o.VP$vals), each = 3),
#         driver = rep(c("env", "space", "bio"), times = ncol(o.VP$vals)),
#         value  = c(o.VP$vals),
#         r2     = rep(o.MF$TjurR2, each = 3)
# )
# o.VP2[r2 < 0, r2 := 0]
# o.VP2[, scaled_values := value * r2]
# #- add stochasticity = 1-r2
# o.VP3 <- data.table(
#         taxon = colnames(o.VP$vals),
#         driver = "stochastic",
#         value = 0,
#         r2 = o.MF$TjurR2
# )
# o.VP3[r2 < 0, r2 := 0]
# o.VP3[, value := 1 - r2]
# o.VP3[, scaled_values := 1 - r2]
# o.VP4 <- rbindlist(list(o.VP2, o.VP3))
# o.VP4$run <- o
# #o.VP4$psrf_mean <- rep(o.means, each = 4)
# #o.VP4$psrf_max  <- rep(o.maxes, each = 4)
# 
# #- new VP to determine relative predictor importance
# o.VP5    <- computeVariancePartitioning(o.m1.mcmc)
# 
# o.VP6 <- rowSums(o.VP5$vals)
# #- drop morans eigenvectors
# o.VP6 <- o.VP6[str_detect(names(o.VP6), "MEM|AEM", negate =
#                                   T)]
# #- drop random sample
# o.VP6 <- o.VP6[str_detect(names(o.VP6), "Random", negate =
#                                   T)]
# o.VP6 <- o.VP6 / sum(o.VP6)
# #o.VP6 <- rank(o.VP6)
# 
# #o.groupnames = c("env", "space")
# 
# #- save fitted model
# o.out.list <-
#         list(
#                 biota       = o.data3,
#                 environment = o.env.samples,
#                 model       = o.m1.mcmc,
#                 VP          = o.VP4,
#                 coordinates = o.xy.mat,
#                 scale       = o.sf,
#                 importance  = o.VP6,
#                 scheme      = o.scheme
#         )
