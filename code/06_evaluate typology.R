# Simulate new data
#' Explantion of variables and functions 
#' VARIABLES 
# setup ------------------------------------------------------------------
library(groundhog)
pkgs <- c(
        "cluster",
        "data.table",
        "dplyr",
        "fs",
        "Hmsc",
        "kernlab",
        "vegclust",
        "magrittr",
        "purrr",
        "stringr"

)
groundhog.library(pkgs,'2024-12-01')
rm(pkgs)


box::use(
        code/custom_functions[adjust_clusters, balance_clusters, ]
)


# setup loop ------------------------------------------------------------------------
model.files <- dir_ls("data/fitted_models/")
n_types = c(4,5,10)
# how many orders to evaluate for zeta diversity 
v.orders <- 10
# number of evaluations model
max.q <- 3
within.q <- 30
neval <- max.q * within.q

# Loop ------------------------------------------------------------------------------
# this could be parallelized
for (i in 1:length(model.files)){
        
        #- Update 
        print(paste("i =",i))
        
        #- Prepare list to store results of loop
        i.out <- vector(mode = "list", length = length(n_types))
        
        #- Load model results 
        i.model <- readRDS(model.files[i])
        
        #- Extract name of model 
        i.model.name <- model.files[i] %>% 
                str_remove("data/fitted_models/")

        #- Extract number of samples in model
        i.nrow  <- nrow(i.model[[2]])
        #- Which variables follow the typology?
        #- Here we make 5 subsets which contain: Four times random assortments and one time 
        #- all variables. 
        #- The strength of the typology is partly determined by how important the variables are. 
        #- This is measured by the fraction of VP scores captured by the included variables. 
        #- These vectors hold information on the simulations 
        i.n.variables <- 
        i.contraction.points <- 
        i.contraction.centroids <- 
        i.importance  <- 
        i.asw <- c()
        i.cluster.assignments <- list()
        i.fuzzy.assignments <- list()
        #--- loop over different number of types 
        for (j in 1:length(n_types)){
                
                #- Update 
                print(paste("i =",i, "j =",j))
                j.types <- n_types[j]
                #  — — — NEW PREDCICTOR VALUES  —  —  —  —
                #- extract predictor names 
                j.all.vars1 <- colnames(i.model[[3]]$X)[-c(1)]
                #- Remove spatial predictors (MEM = MORANS EIGENVECTOR MAPS) 
                j.all.vars1 <- j.all.vars1[!str_detect(j.all.vars1, "MEM")]
                j.all.vars1 <- j.all.vars1[!str_detect(j.all.vars1, "AEM")]
                #- create a list to store results of different variable sets
                j.predictions <- list()

                #- START LOOP OVER Q, varies the number of variables that 
                #- follow the typology
                for (q in 1:max.q){
                        print(paste("i =", i, "j =", j, "q =", q))
                        
                        # How many variables are included in this subset? 
                        q.nvariables <- sample(x = 2:(length(j.all.vars1) - 1),
                                               size = 1)
                        # Add number to vector for later reference
                        i.n.variables %<>% append(q.nvariables)
                        # Which variables are included 
                        q.all.vars   <- sample(j.all.vars1, 
                                size = q.nvariables, 
                                replace = F)
                        # Vector with names of all variables that do not follow 
                        # the typology
                        q.missing  <- j.all.vars1[which(!j.all.vars1 %in% q.all.vars)]
                        
                        # At this point we can already determine the typology 
                        # strength as judged by the selection of variables.
                        q.importance <- sum(i.model$importance[which(names(i.model$importance) %in% q.all.vars)])
                        i.importance %<>% append(q.importance) 
                        
                        q.all.vars.list <- list()

                        q.env <- copy(i.model$environment)
                        
                        #- cluster based on environment 
                        q.cluster.env <- copy(q.env)
                        q.cluster.env %<>%
                                select(all_of(q.all.vars)) %>% 
                                as.data.frame %>% 
                                as.matrix
                        q.type <- specc(x = q.cluster.env, centers = j.types)
                        q.clusters <- balance_clusters(
                                data = q.cluster.env,
                                clusters = q.type@.Data,
                                min_size = nrow(q.cluster.env) /
                                        j.types * 0.75
                        ) 

                        
                        #- store cluster assignment for later 
                        i.cluster.assignments[[q]] <- q.clusters
                        #j.fuzzy.assignments[[q]]   <- q.fc.b
                        q.env[, type := q.type]
                        #q.env[, type.fuzzy := q.fc.b]
                        # find centroid for each variable that follows the typology
                        q.centroids <- q.env[, lapply(.SD, mean), .SDcols = q.all.vars, by = "type"]
                        setorderv(q.centroids, "type")
                        q.centroids[, type := NULL]
                        q.centroids %<>% as.data.frame %>% as.matrix
                        
                        # prepare environmental variables for adjustments in script below 
                        q.observations <- as.matrix(as.data.frame(q.cluster.env))
                        
                        # TODO describe 
                        source("code/new_predictors_for_hard_classification.R")
                        
                        # fuzzy classification 
                        q.fuzzy <- lapply(1:within.q, 
                                function(x) {
                                        vegclust(q.newenv[[x]][, which(colnames(q.newenv[[x]]) %in% q.all.vars)],
                                        mobileCenters = j.types,
                                        method = "FCM",
                                        m = 1.5
                                        )
                                }        
                        )
                        i.fuzzy.assignments[[length(i.fuzzy.assignments) + 1]] <- q.fuzzy
                        # there used to be a balancing step here. 
                        # I removed it
                        # Sine the spectral clustering is already balanced, and the fuzzy classification uses 
                        # the results of the spectral classification, it is likely not necesarry.
                        q.posterior_samples <- poolMcmcChains(i.model[[3]]$postList)
                        q.selected_samples  <- q.posterior_samples[c(1, length(q.posterior_samples))]
                        #  — — — predict new biota  —  —  —  —
                        q.out <- lapply(q.newenv, 
                                        FUN = function(pre) {
                                                x <- predict(object = i.model[[3]], 
                                                        post = q.selected_samples, 
                                                        X = pre
                                                    )
                                                x <- x[[2]]
                                                return(x)
                                        }
                                        )
                        j.predictions [[q]] <- q.out
                        rm(list = ls()[grepl("^q\\.", x = ls())])
                }
                rm(q)
                i.out[[j]] <- j.predictions
                rm(list = ls()[grepl("^j\\.", x = ls())])
        }
        rm(j)
} 



# rm(i)
# out.all <- rbindlist(out)

# out.all %>% 
#         filter(metric == "cs") %>% 
#         ggplot(aes(y = value , x = importance)) + geom_point(aes(col = contraction)) + geom_smooth()
# 
# saveRDS(out.all, "data/250115_results.rds")
# rm(out.all, out, model.files, n_types)

