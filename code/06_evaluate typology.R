
#' Explantion of variables and functions 
#' VARIABLES 



# setup loop ------------------------------------------------------------------------
model.files <- fs::dir_ls("data/fitted_models/")
n_types = 2:6
out <- list()
spatial_constrained <- TRUE
#- how many orders to evaluate for zeta diversity 
v.orders <- 10

model.files <- c("data/fitted_models/1.rds", "data/fitted_models/2.rds", "data/fitted_models/3.rds")

# Loop ------------------------------------------------------------------------------
for (i in 1:length(model.files)){
        
        print(paste("i =",i))
        
        i.model <- readRDS(model.files[i])
        i.nrow  <- nrow(i.model[[1]])
        i.impo  <- ceiling(length(i.model$importance)/3)
        i.impo2 <- sort(i.model$importance)
        i.impo3 <- list(
                i.impo2[length(i.impo2):1][1:i.impo],
                i.impo2[(1+i.impo):(2*i.impo)],
                i.impo2[1:i.impo]
        )
        #--- loop over different number of types 
        for (j in 1:length(n_types)){
                print(paste("j =",j))
                j.types <- n_types[j]
        
                #--- spatially constrain type membership
                if (spatial_constrained == TRUE){

                        #- create sf object to create spatial cluster
                        j.tv <- st_as_sf(data.frame(i.model$coordinates), coords = c("x", "y"))
                        j.tv <- st_distance(j.tv)
                        j.tv <- as.dist(j.tv)
                        j.tv <- cluster::agnes(j.tv, method = "complete")
                        j.tv <- cutree(j.tv, k = j.types)
                        
                #--- spatially unconstrained type membership
                } else { 
                        #- create type vector to count instances 
                        j.tv <- sort(rep(1:j.types, times = 50)[1:100])
                }
                
                #  — — — NEW PREDCICTOR VALUES  —  —  —  —
                #- extract predictor names 
                j.all.vars1 <- colnames(i.model[[3]]$X)[-c(1)]
                #- Remove spatial predictors 
                j.all.vars1 <- j.all.vars1[!str_detect(j.all.vars1, "MEM")]
                
                #- create a list to store results of different variable sets
                j.storage <- list()
                
                #- which variables follow the typology 
                for (q in 1:3){
                        print(paste("q =",q))
                        q.impo     <- i.impo3[[q]]
                        q.all.vars <- j.all.vars1[j.all.vars1 %in% names(q.impo)]
                        q.missing  <- j.all.vars1[which(!j.all.vars1 %in% q.all.vars)]
                        
                        
                        q.all.vars.list <- list()       
                        for (s in 1:3){
                                print(paste("s =",s))
                                #- simulate environmental variables using the simulate_env function defines in the script 00_setup.
                                #- this function takes a model object (x)
                                s.all.vars2 <- 
                                        lapply(q.all.vars, 
                                               function(lapp.x) simulate_env(
                                                       x = i.model, 
                                                       y = lapp.x, 
                                                       types = j.types,
                                                       strength = s)
                                        ) %>% 
                                        do.call(cbind, .) %>%
                                        data.frame
                                
                                
                                names(s.all.vars2)   <- q.all.vars
                                #s.all.vars2$strength <- 
                                q.all.vars.list[[s]] <- s.all.vars2
                                rm(list = ls()[grepl("^s\\.", x = ls())])
                        } # END LOOP OVER DIFFERENT STRENGTHs 
                        rm(s)
                        #- add missing variables. These values are also simulated but from one common distribution. 
                        q.all.vars2 <- 
                                lapply(q.missing, 
                                       function(lapp.x) simulate_env(
                                               x     = i.model, 
                                               y     = lapp.x, 
                                               types = 1,
                                               strength = 1,
                                               obs_per_type = 100)
                                ) %>% 
                                do.call(cbind, .)%>%
                                data.frame
                        names(q.all.vars2) <- q.missing
                        
                        q.all.vars3 <- lapply(q.all.vars.list, FUN = function(x) cbind(x, q.all.vars2))
                        
                        q.new_x <- lapply(
                                X = q.all.vars3,
                                FUN = function(x)
                                        as.matrix(cbind(1, x, select(
                                                i.model[[2]], contains("MEM")
                                        ))
                                        )
                                )
        
                        # j.new_x <- 
                        #         cbind(
                        #                 1, 
                        #                 j.all.vars2, 
                        #                 select(i.model[[2]], 
                        #                         contains("MEM")
                        #                 )
                        #         ) %>% 
                        #         as.matrix()
        
                        #  — — — predict new biota  —  —  —  —
                        q.out <- lapply(q.new_x, 
                                        FUN = function(pre) predict(object = i.model[[3]],
                                                                  X = pre))
                        
                        # j.out    <- predict(
                        #                         object = i.model[[3]],
                        #                         X      = j.new_x
                        #  
                        #)
                        #- this created many predictions - sample 10
                        #- a list of three (one for each strength) with a list of 10 inside. 
                        #- I.e. 10 data sets per strength.
                        q.out2 <- lapply(q.out, function(x) x[sample(1:length(x), 10)])
                        j.storage[[q]] <- q.out2
                        rm(list = ls()[grepl("^q\\.", x = ls())])
                }
                rm(q)

                #- compute distance matrix for each simulated community
                j.distance.matrix  <- 
                        rapply(
                                j.storage, 
                                parallelDist, 
                                method = "binary",
                                how = "list"
                        )

                #- compute Classification strength
                j.cs <- rapply(j.distance.matrix, my.cs)
                j.cs <- data.table(
                        value = j.cs, 
                        simulation = rep(1:10, times = 9),
                        q = rep(c(1,2,3), each = 30),
                        s = rep(rep(c(1,2,3), each = 10), times = 3),
                        metric = "cs"
                )

                #- compute ANOSIM 
                j.anosim <- rapply(j.distance.matrix, my.anosim)
                j.anosim <- data.table(
                        value = j.anosim, 
                        simulation = rep(rep(1:10, each = 2), times = 9),
                        q = rep(rep(c(1,2,3), each = 2), each = 30),
                        s = rep(rep(rep(c(1,2,3), each = 2), each = 10), times = 3),
                        metric = c("anosim.stat", "anosim.p")
                )
                
                
                # x.1 ——— compute area under the zeta diversity decline curve
                j.ut  <- unique(j.tv)
                j.ut2 <- table(j.tv)
                if (any(j.ut2 < 11)){
                        j.ut <- j.ut[-which(j.ut2<11)]
                }
                j.zeta <- list()
                #- LOOP OVER unique types. One AUCζ is computed for each type.  
                for (k in seq_along(j.ut)){
                        print(paste("k =",k))
                        
                        k.x <- 
                                rapply(j.storage,
                                       function(x) Zeta.decline.ex(data.spec = x[which(j.tv == j.ut[k]),],
                                                                  orders     = 1:v.orders,
                                                                  plot       = FALSE,
                                                                  rescale    = TRUE ),
                                       how = "list"
                                )
                        
                        k.x2 <- lapply(k.x, 
                                       function(x1) {
                                               lapply(x1, function(x2) {
                                                       lapply(x2, function(x3) {
                                                               x3[["zeta.val"]]
                                        })
                                })
                        })
                        
                        k.x2 <-
                                rapply(k.x2, function(x) calculate_auc(y = x, x = 1:v.orders))%>%
                                #do.call(rbind, .) %>%
                                data.frame %>%
                                mutate(type = k,
                                       simulation = rep(1:10, times = 9),
                                       q = rep(c(1,2,3), each = 30),
                                       s = rep(rep(c(1,2,3), each = 10), times = 3),
                                       metric = "auc_zeta"
                                       )
                                
                        k.x2 <- rename(k.x2, value = ".")
                        j.zeta[[k]] <- k.x2
                        rm(list = ls()[grepl("^k\\.", x = ls())])
                }
                rm(k)
                j.zeta %<>% rbindlist
                j.zeta[, value := mean(value), by = c("simulation", "q", "s")]
                j.zeta %<>% unique(by = c("simulation", "q", "s"))
                j.zeta[, type := NULL]
                
                
                #  ——— Compute Permanova 
                j.permanova <- rapply(j.distance.matrix, function(x) adonis2(formula = x ~j.tv), how = "list")
                j.r2 <- lapply(j.permanova, 
                               function(x1) {
                                       lapply(x1, function(x2) {
                                               lapply(x2, function(x3) {
                                                       x3[["R2"]][1]
                                               })
                                       })
                               }) %>% unlist
                j.p <- lapply(j.permanova, 
                               function(x1) {
                                       lapply(x1, function(x2) {
                                               lapply(x2, function(x3) {
                                                       x3[["Pr(>F)"]][1]
                                               })
                                       })
                               }) %>% unlist
                j.F <- lapply(j.permanova, 
                               function(x1) {
                                       lapply(x1, function(x2) {
                                               lapply(x2, function(x3) {
                                                       x3[["F"]][1]
                                               })
                                       })
                               }) %>% unlist
                
                j.r2 <- data.table(
                        value = j.r2, 
                        simulation = rep(1:10, times = 9),
                        q = rep(c(1,2,3), each = 30),
                        s = rep(rep(c(1,2,3), each = 10), times = 3),
                        metric = "PERMANOVA R2"
                )
                j.p <- data.table(
                        value = j.p, 
                        simulation = rep(1:10, times = 9),
                        q = rep(c(1,2,3), each = 30),
                        s = rep(rep(c(1,2,3), each = 10), times = 3),
                        metric = "PERMANOVA p"
                )
                j.F <- data.table(
                        value = j.F, 
                        simulation = rep(1:10, times = 9),
                        q = rep(c(1,2,3), each = 30),
                        s = rep(rep(c(1,2,3), each = 10), times = 3),
                        metric = "PERMANOVA F"
                )
                
                # silhouette width 
                j.asw <- rapply(j.distance.matrix, my.asw)
                j.asw <- data.table(
                        value = j.asw, 
                        simulation = rep(1:10, times = 9),
                        q = rep(c(1,2,3), each = 30),
                        s = rep(rep(c(1,2,3), each = 10), times = 3),
                        metric = "asw"
                )
                
                
                # Indicator value 
                 j.isa <- rapply(j.storage, my.isa)
                 j.isa <- data.table(
                         value = j.isa, 
                         simulation = rep(rep(1:10, each = 2), times = 9),
                         q = rep(rep(c(1,2,3), each = 2), each = 30),
                         s = rep(rep(rep(c(1,2,3), each = 2), each = 10), times = 3),
                         metric = c("isa.number", "isa.avg.p")
                 )
                
                # ISAMIC
                 j.isamic <- rapply(j.storage, my.isamic)
                 j.isamic <- data.table(
                         value = j.isamic, 
                         simulation = rep(1:10, times = 9),
                         q = rep(c(1,2,3), each = 30),
                         s = rep(rep(c(1,2,3), each = 10), times = 3),
                         metric = "isamic"
                 )
                
                #  — — — Compile output into a single table  — — —
                j.out.final <- 
                        rbindlist(
                                list(
                                        j.cs,
                                        j.zeta,
                                        j.anosim,
                                        j.r2,
                                        j.p,
                                        j.F,
                                        j.asw,
                                        j.isa,
                                        j.isamic
                                )
                        ) %>% 
                        .[, `:=` (run = i, types = j.types)]
                
                # j.out.final <- data.table(
                #                                 run = i, 
                #                                 types = j.types, 
                #                                 value = c(j.cs, 
                #                                           j.zeta$auczeta, 
                #                                           j.r2,
                #                                           j.F,
                #                                           j.p,
                #                                           j.anosimstat,
                #                                           j.anosimp),
                #                                 statistic = rep(c("CS", 
                #                                                   "AUCzeta",
                #                                                   "PERMANOVA R2",
                #                                                   "PERMANOVA F",
                #                                                   "PERMANOVA p",
                #                                                   "ANOSIM R",
                #                                                   "ANOSIM p"),each = 10)
                #                         )
                out[[length(out) + 1]] <- j.out.final
                rm(list = ls()[grepl("^j\\.", x = ls())])
        } # END OF LOOP j 
        rm(j)
        rm(list = ls()[grepl("^i\\.", x = ls())])
} # END OF LOOP i 
rm(i)
out.all <- rbindlist(out)
saveRDS(out.all, "data/241220_first_results.rds")
rm(out.all, out, model.files, n_types)

# library(ggplot2)
# out.all %>% 
#         ggplot(aes(x = model, y = value, group = model)) + 
#         geom_violin(aes(col = factor(model)), draw_quantiles = .5) + 
#         geom_jitter(width = 0.25, alpha = 0.2) +
#         facet_wrap(.~statistic, scales = "free")

# x.1 <- 
# x1 <- readRDS(model.files[1]) %>% .$VP %>% group_by(driver) %>% summarise(mean = mean(scaled_values)) %>% mutate(model = 1)
# x2 <- readRDS(model.files[2]) %>% .$VP %>% group_by(driver) %>% summarise(mean = mean(scaled_values)) %>% mutate(model = 2)
# x3 <- readRDS(model.files[3]) %>% .$VP %>% group_by(driver) %>% summarise(mean = mean(scaled_values)) %>% mutate(model = 3)
# x4 <- readRDS(model.files[4]) %>% .$VP %>% group_by(driver) %>% summarise(mean = mean(scaled_values)) %>% mutate(model = 4)

# xx <- bind_rows(x1, x2, x3, x4)
# xx2 <- tidyr::pivot_wider(xx, id_cols = model, names_from = driver, values_from = mean)


# out.all2 <- left_join(out.all, xx2, by = "model") 

# out.all2 %>% 
#         ggplot(aes(x = value, y = env)) + 
#         #geom_violin(aes(col = factor(model)), draw_quantiles = .5) + 
#         geom_jitter(width = 0.25, alpha = 0.2) +
#         facet_wrap(.~statistic, scales = "free")

# # load data ------------------------------------------------------------------------------------
# x <- readRDS("data/fitted_models/1.rds")
# 
# # create new values for predictor variables ----------------------------------------------------
# n_types <- 4
# n_row <- nrow(x[[1]])
# n_obs_per_type <- n_row/n_types
# 
# 
# 
# 
# # test <- simulate_env(y = "soil_pH")
# # test <- data.frame(type = rep(1:4, each =25), ph = test)
# # ggplot(test, aes(y = ph, x = type, group = type)) + geom_boxplot()
# 
# 
# #- sites should be in one type if they are close #- at least at landscape and ecoregion scale 
# #- add this later
# 
# 
# 
# #new_x <- x[[3]]$XData[1:3, -c(1,2)] %>% data.frame %>% as.matrix
# 
# 
# #- ten random samples for out
# 

