# d <- dir_ls("data/evaluations/")
# d <- lapply(d, readRDS)
# d <- lapply(d, rbindlist)
# d <- rbindlist(d)
# saveRDS(d, "rf_test_data.rds")


# setup -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())

library(data.table)
library(tidymodels)
library(rpart.plot)  # for visualizing a decision tree
library(vip)   
library(pdp)
library(directlabels)
# load data ---------------------------------------------------------------
d <- readRDS("data/VP_results_polished.rds")
e <- rbindlist(lapply(list.files("data/evaluations/", full.names = T), readRDS))

# prepare data ------------------------------------------------------------
de <- d[e, on = "scheme_id"]

# prepare model -----------------------------------------------------------
# Some parts of the model are agnostic toward which subset it is fit to and can be defined outside the loop. 
# Define a random forest model from {parsnip}. 
# The model will try to predict the different evaluation criteria
tree_tune <-
        rand_forest(
                mode   = "regression",
                trees  = tune(),
                mtry   = tune(),
                min_n  = tune()
        )
tree_tune <-
        set_engine(
                object = tree_tune,
                engine = "ranger",
                importance = "permutation"
        )
# Parameter grid for tuning         
predictor_vec <- c("env", "bio", "space", "stochastic", "samples", "species_rank", "genus_rank", "family_rank", "higher_rank", "min_distance", "median_distance", "max_distance", "n_types", "variables", "contraction_points", "contraction_centroids", "env_asw", "importance")


tree_grid <- grid_regular(
        trees(),
        finalize(mtry(),  x = de[, ..predictor_vec]),
        finalize(min_n(), x = de[, ..predictor_vec]),
        levels = 5)

#- create workflow 
tree_wf <- workflow() 
tree_wf <- add_model(x = tree_wf, spec = tree_tune)
tree_wf <- add_formula(x = tree_wf, value ~ variables + n_types + contraction_points + contraction_centroids + env_asw + importance + env + bio + space + stochastic + samples + species_rank + genus_rank + family_rank + higher_rank + min_distance + median_distance + max_distance)
vip.scores <- list()
# Loop Over taxa 
for (i in 1:4){
        
        i.d <- which(de$taxon == c("diatoms", "fish", "invertebrates", "macrophytes")[i])
        i.d <- de[i.d, ]
        
        # Loop over evaluation criteria 
        for (j in 1:uniqueN(d$metric)){
                
                j.d <- i.d[which(i.d$metric == unique(de$metric)[j]), ]
                
                set.seed(1)
                j.split <- initial_split(data = j.d)
                j.train <- training(j.split)
                j.folds <- vfold_cv(j.train)
                #- perform model tuning via grid search on tree_grid
                j.tune <- 
                        tune_grid(
                                object    = tree_wf,
                                resamples = j.folds,
                                grid      = tree_grid
                        )
                j.best      <- select_best(j.tune, metric = "rmse")
                j.wf        <- finalize_workflow(tree_wf, j.best)
                j.final.fit <- last_fit(j.wf, j.split)
                j.tree      <- extract_workflow(j.final.fit)
                j.model.fit <- extract_fit_parsnip(j.tree) 
                # Plots 
                ## Variance Importance Plot 
                j.vip.scores <- j.model.fit$fit$variable.importance
                vip.scores[[length(vip.scores) + 1]] <- data.frame(vip = j.vip.scores, 
                                           organism = c("diatoms", "fish", "invertebrates", "macrophytes")[i], 
                                           metric = unique(de$metric)[j], 
                                           scheme = unique(i.d$scheme_id),
                                           model.MSE = j.model.fit$fit$prediction.error,
                                           model.R2 = j.model.fit$fit$r.squared)
                j.vip.plot <-  vip(j.model.fit)
                ggsave(filename= paste0("figures/RF variance importance/",c("diatoms", "fish", "invertebrates", "macrophytes")[i],"_",unique(d$metric)[j],".png"))
                cut(j.d)
                ggplot(j.d, aes(x = value )) + geom_density() + geom_rug() + facet_wrap(.~contraction_centroids)
                
                
                
                covariate_vector <- c("variables",
                                      "n_types",
                                      "contraction_centroids",
                                      "env_asw",
                                      "importance",
                                      "env",
                                      "bio",
                                      "space",
                                      "stochastic",
                                      "contraction_points")
                covlist <- list()
                for (k in seq_along(covariate_vector)){
                        
                        k.pdp <- pdp::partial(j.model.fit, pred.var = covariate_vector[k], train = j.train)  
                        k.df <- k.pdp
                        k.df$pred.var <- covariate_vector[k]
                        names(k.df)[1] <- "value"
                        k.df[,1] <- scale(k.df$value)
                        covlist[[k]] <- k.df
                }
                pdp.plot <- rbindlist(covlist)
                pdp.plot <- direct.label(ggplot(pdp.plot, aes(x=value, y = yhat, col = pred.var)) + geom_smooth(se = FALSE) + xlim(-2, 3), "last.points")
                ggsave(plot = pdp.plot, 
                       filename = paste0("figures/partial dependence plots/",c("diatoms", "fish", "invertebrates", "macrophytes")[i],"_",unique(d$metric)[j],".png")
                       )    
        }
        
}
# # final_tree %>% 
# #         extract_fit_parsnip() %>% 
# #         vip()
# # ggsave(filename= "figures/RF variance importance/isamic_varimpo.png")
# final_tree %>% 
#         extract_fit_parsnip() 
#         
#        
# p1 <- pdp::partial(ef, pred.var = "variables", train = cs_train)  
# p2 <- pdp::partial(ef, pred.var = "n_types", train = cs_train)  
# p3 <- pdp::partial(ef, pred.var = "contraction_centroids", train = cs_train)  
# p4 <- pdp::partial(ef, pred.var = "env_asw", train = cs_train)  
# p5 <- pdp::partial(ef, pred.var = "importance", train = cs_train)  
# p6 <- pdp::partial(ef, pred.var = "env", train = cs_train)  
# p7 <- pdp::partial(ef, pred.var = "bio", train = cs_train)  
# p8 <- pdp::partial(ef, pred.var = "space", train = cs_train)  
# p9 <- pdp::partial(ef, pred.var = "stochastic", train = cs_train)  
# p10 <- pdp::partial(ef, pred.var = "contraction_points", train = cs_train)  
# 
# 
# 
# 
# ggplot(aes(x = contraction_points, y = yhat)) + geom_line(linewidth = 2) + ylim(0,1)
# ggsave(filename= "figures/partial dependence plots/isamic_pd_cp.png")
# 
# 
# final_tree %>% 
#         extract_fit_parsnip() %>%
#         pdp::partial(pred.var = "scale_min", train = cs_train) %>% 
#         ggplot(aes(x = scale_min, y = yhat)) + geom_line()
# ggsave(filename= "figures/partial dependence plots/isamic_pd_scale_min.png")
# final_tree %>% 
#         extract_fit_parsnip() %>%
#         pdp::partial(pred.var = "env_asw", train = cs_train) %>% 
#         ggplot(aes(x = env_asw, y = yhat)) + geom_line()
# ggsave(filename= "figures/partial dependence plots/isamic_pd_env_asw.png")
# 
# 
# cs <- filter(d, metric == "auc_zeta")
# 
# tree_tune <- 
#         rand_forest(trees = tune(), 
#                     mtry = tune(),
#                     min_n = tune()) %>% 
#         set_engine("ranger", importance = "permutation") %>% 
#         set_mode("regression")
# tree_grid <- grid_regular(
#         trees(),
#         finalize(mtry(), x = d[, -1]),
#         finalize(min_n(), x = d[, -1]),
#         levels = 5)
# 
# set.seed(234)
# cs_split <- initial_split(data = cs)
# cs_train <- training(cs_split)
# cs_folds <- vfold_cv(cs_train)
# 
# 
# #- create workflow 
# tree_wf <- workflow() %>%
#         add_model(tree_tune) %>%
#         add_formula(value ~ variables + contraction_points + contraction_centroids + env_asw + importance + scale_median + scale_min + scale_max + types + env + bio + space + stochastic)
# 
# tree_res <- 
#         tree_wf %>% 
#         tune_grid(
#                 resamples = cs_folds,
#                 grid      = tree_grid
#         )
# tree_res %>%
#         show_best(metric = "rmse")
# best_tree <- tree_res %>%
#         select_best(metric = "rmse")
# final_wf <- 
#         tree_wf %>% 
#         finalize_workflow(best_tree)
# final_fit <- 
#         final_wf %>%
#         last_fit(cs_split)
# final_tree <- extract_workflow(final_fit)
# # final_tree %>%
# #         extract_fit_engine() %>%
# #         rpart.plot(roundint = FALSE)
# final_tree %>% 
#         extract_fit_parsnip() %>% 
#         vip()
# ggsave(filename= "figures/RF variance importance/auc_varimpo.png")
# final_tree %>% 
#         extract_fit_parsnip() %>%
#         pdp::partial(pred.var = "env_asw", train = cs_train) %>% 
#         ggplot(aes(x = env_asw, y = yhat)) + geom_line(linewidth = 2)
# ggsave(filename= "figures/partial dependence plots/auc_pd_asw.png")
# final_tree %>% 
#         extract_fit_parsnip() %>%
#         pdp::partial(pred.var = "env", train = cs_train) %>% 
#         ggplot(aes(x = env, y = yhat)) + geom_line()
# ggsave(filename= "figures/partial dependence plots/auc_pd_env.png")
# final_tree %>% 
#         extract_fit_parsnip() %>%
#         pdp::partial(pred.var = "contraction_points", train = cs_train) %>% 
#         ggplot(aes(x = contraction_points, y = yhat)) + geom_line()
# ggsave(filename= "figures/partial dependence plots/cs_pd_cp.png")
