# d <- dir_ls("data/evaluations/")
# d <- lapply(d, readRDS)
# d <- lapply(d, rbindlist)
# d <- rbindlist(d)
# saveRDS(d, "rf_test_data.rds")

d <- readRDS("data/eval_w_dirichlet.rds")

library(tidymodels)
library(rpart.plot)  # for visualizing a decision tree
library(vip)   
library(pdp)

cs <- filter(d, metric == "isamic")

tree_tune <- 
                rand_forest(trees = tune(), 
                            mtry = tune(),
                            min_n = tune()) %>% 
                set_engine("ranger", importance = "permutation") %>% 
                set_mode("regression")
tree_grid <- grid_regular(
                          trees(),
                          finalize(mtry(), x = d[, -1]),
                          finalize(min_n(), x = d[, -1]),
                          levels = 5)

set.seed(234)
cs_split <- initial_split(data = cs)
cs_train <- training(cs_split)
cs_folds <- vfold_cv(cs_train)


#- create workflow 
tree_wf <- workflow() %>%
        add_model(tree_tune) %>%
        add_formula(value ~ variables + contraction_points + contraction_centroids + env_asw + importance + scale_median + scale_min + scale_max + types + env + bio + space + stochastic)

tree_res <- 
        tree_wf %>% 
        tune_grid(
                resamples = cs_folds,
                grid = tree_grid
        )
tree_res %>%
        show_best(metric = "rmse")
best_tree <- tree_res %>%
        select_best(metric = "rmse")
final_wf <- 
        tree_wf %>% 
        finalize_workflow(best_tree)
final_fit <- 
        final_wf %>%
        last_fit(cs_split)
final_tree <- extract_workflow(final_fit)
# final_tree %>%
#         extract_fit_engine() %>%
#         rpart.plot(roundint = FALSE)
final_tree %>% 
        extract_fit_parsnip() %>% 
        vip()
ggsave(filename= "figures/RF variance importance/isamic_varimpo.png")
final_tree %>% 
        extract_fit_parsnip() %>%
        pdp::partial(pred.var = "contraction_points", train = cs_train) %>% 
        ggplot(aes(x = contraction_points, y = yhat)) + geom_line(linewidth = 2)
ggsave(filename= "figures/partial dependence plots/isamic_pd_cp.png")
final_tree %>% 
        extract_fit_parsnip() %>%
        pdp::partial(pred.var = "scale_min", train = cs_train) %>% 
        ggplot(aes(x = scale_min, y = yhat)) + geom_line()
ggsave(filename= "figures/partial dependence plots/isamic_pd_scale_min.png")
final_tree %>% 
        extract_fit_parsnip() %>%
        pdp::partial(pred.var = "env_asw", train = cs_train) %>% 
        ggplot(aes(x = env_asw, y = yhat)) + geom_line()
ggsave(filename= "figures/partial dependence plots/isamic_pd_env_asw.png")


cs <- filter(d, metric == "auc_zeta")

tree_tune <- 
        rand_forest(trees = tune(), 
                    mtry = tune(),
                    min_n = tune()) %>% 
        set_engine("ranger", importance = "permutation") %>% 
        set_mode("regression")
tree_grid <- grid_regular(
        trees(),
        finalize(mtry(), x = d[, -1]),
        finalize(min_n(), x = d[, -1]),
        levels = 5)

set.seed(234)
cs_split <- initial_split(data = cs)
cs_train <- training(cs_split)
cs_folds <- vfold_cv(cs_train)


#- create workflow 
tree_wf <- workflow() %>%
        add_model(tree_tune) %>%
        add_formula(value ~ variables + contraction_points + contraction_centroids + env_asw + importance + scale_median + scale_min + scale_max + types + env + bio + space + stochastic)

tree_res <- 
        tree_wf %>% 
        tune_grid(
                resamples = cs_folds,
                grid = tree_grid
        )
tree_res %>%
        show_best(metric = "rmse")
best_tree <- tree_res %>%
        select_best(metric = "rmse")
final_wf <- 
        tree_wf %>% 
        finalize_workflow(best_tree)
final_fit <- 
        final_wf %>%
        last_fit(cs_split)
final_tree <- extract_workflow(final_fit)
# final_tree %>%
#         extract_fit_engine() %>%
#         rpart.plot(roundint = FALSE)
final_tree %>% 
        extract_fit_parsnip() %>% 
        vip()
ggsave(filename= "figures/RF variance importance/auc_varimpo.png")
final_tree %>% 
        extract_fit_parsnip() %>%
        pdp::partial(pred.var = "env_asw", train = cs_train) %>% 
        ggplot(aes(x = env_asw, y = yhat)) + geom_line(linewidth = 2)
ggsave(filename= "figures/partial dependence plots/auc_pd_asw.png")
final_tree %>% 
        extract_fit_parsnip() %>%
        pdp::partial(pred.var = "env", train = cs_train) %>% 
        ggplot(aes(x = env, y = yhat)) + geom_line()
ggsave(filename= "figures/partial dependence plots/auc_pd_env.png")
final_tree %>% 
        extract_fit_parsnip() %>%
        pdp::partial(pred.var = "contraction_points", train = cs_train) %>% 
        ggplot(aes(x = contraction_points, y = yhat)) + geom_line()
ggsave(filename= "figures/partial dependence plots/cs_pd_cp.png")
