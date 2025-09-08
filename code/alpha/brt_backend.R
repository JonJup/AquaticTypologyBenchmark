# Load required packages
library(groundhog)
pkgs <- c(
        "data.table",
        "tidymodels",
        "themis",
        "bonsai",
        "DALEXtra",
        "tidyr",
        "dplyr",
        "ggplot2",
        "xgboost"
)
groundhog.library(pkgs, "2024-12-01")
rm(pkgs)

# # Create example data (replace with your actual data)
# set.seed(123)
# n <- 1000

# data <- tibble(
#         simulation_id = 1:n,
#         scenario_id = rep(1:(n/10), each = 10),
#         spatial_scale = rnorm(n),
#         #taxonomic_group = factor(sample(c("plants", "invertebrates", "vertebrates"), n, replace = TRUE)),
#         metric_1 = rnorm(n),
#         metric_2 = rnorm(n),
#         metric_3 = rnorm(n),
#         metric_4 = rnorm(n),
#         quality_level = runif(n)
# )


data <- fs::dir_ls("data/evaluations/")
data <- lapply(data, readRDS) %>% rbindlist
typicality <- readRDS("data/dirichlet_weights.rds")
typicality <- typicality[1:4, ]
names(typicality)[1] <- "model"
# Prepare variable entries in "metric" to be suited for column names 
data[metric == "classification strength", metric := "classification_strength"]

# Pivot wider to have one column for each explanatory variable 
data <- data %>% 
        pivot_wider(
                id_cols = c("model", "n_types", "variables", "contraction_points", "contraction_centroids", "env_asw", "importance"), 
                names_from = metric, 
        values_from = value)

# Join typicality values to main data
data <- left_join(data, typicality, by = "model")

# This step will use the typicality value as importance weights
data <- mutate(data, normalized_dirichlet_density = importance_weights(normalized_dirichlet_density))
# Split data into training and testing sets
set.seed(234)
data_split <- initial_split(data, prop = 0.8)
train_data <- training(data_split)
test_data  <- testing(data_split)


recepie_cs <- recipe(
        x = data, 
        vars    = c("classification_strength", "contraction_points", "contraction_centroids", "env_asw", "importance", "n_types","normalized_dirichlet_density"), 
        roles   = c("outcome", "predictor", "predictor", "predictor", "predictor", "predictor", "case_weights")
)


# Create a recipe for preprocessing
model_recipe <- recipe(
                        x = train_data,
                        vars = c("quality_level", "typicality_value", "metric_1", "metric_2",  "metric_3",  "metric_4",  "spatial_scale"),
                        roles = c("outcome", "case_weight", "predictor", "predictor", "predictor", "predictor", "predictor")
                        ) 
prepped_recipe  <- prep(recepie_cs)

# Create a boosted tree model specification
brt_spec <- boost_tree(
        trees = tune(),
        tree_depth = tune(),
        learn_rate = tune(),
        min_n = tune(),
        loss_reduction = tune(),
        sample_size = tune()
) %>%
        set_engine("xgboost") %>%
        set_mode("regression")

# Create a workflow that combines the recipe and model
brt_workflow <- workflow() %>%
        add_recipe(recepie_cs) %>%
        add_model(brt_spec) %>% 
        add_case_weights(normalized_dirichlet_density)

# Setup Cross validation folds 
set.seed(345)
cv_folds <- vfold_cv(train_data, v = 5)

# Define the parameter grid
xgb_params <- parameters(
        trees(range = c(100, 1000)),
        tree_depth(range = c(2, 10)),
        learn_rate(range = c(-3, -1), trans = log10_trans()),
        min_n(range = c(2, 20)),
        loss_reduction(range = c(-4, 0), trans = log10_trans()),
        sample_prop(range = c(0.5, 1.0))
      )

# Create a grid of hyperparameter combinations
set.seed(456)
xgb_grid <- grid_space_filling(
  xgb_params,
  size = 20  # Number of parameter combinations to try
)
# Tune the model
set.seed(567)
tune_results <- tune_grid(
  brt_workflow,
  resamples = cv_folds,
  grid = xgb_grid,
  metrics = metric_set(rmse, rsq, mae),
  control = control_grid(save_pred = TRUE,
                         save_workflow = TRUE,
                         verbose = TRUE)
)

# Visualize tuning results
autoplot(tune_results)

# Show the best parameters
show_best(tune_results, metric = "rmse", n = 5)
# Select the best model based on RMSE
best_params <- select_best(tune_results, metric = "rmse")

# Finalize the workflow with the best parameters
final_workflow <- finalize_workflow(
        brt_workflow,
        best_params
      )

# Train the final model
set.seed(789)
final_fit <- final_workflow %>%
  fit(data = train_data)

# Extract the fitted model
xgb_model <- final_fit %>%
  extract_fit_engine()

# Evaluate the final model on the test set
final_results <- final_workflow %>%
  last_fit(data_split)

# Show test set metrics
collect_metrics(final_results)
vip_plot <- xgb_model %>%
        vip::vip(num_features = 10, geom = "point") +
        theme_minimal() +
        labs(title = "Variable Importance in the Boosted Model",
             subtitle = "Based on model's internal importance metric")
print(vip_plot)


# # Train the model with importance weights
# set.seed(345)
# brt_fit <- brt_workflow %>%
#         fit(data = train_data)

# # Extract the fitted model
# xgb_model <- brt_fit %>%
#         extract_fit_engine()





# # Variable importance based on the model's internal metrics
# #TODO move importance
# vip_plot <- xgb_model %>%
#         vip::vip(num_features = 10, geom = "point") +
#         theme_minimal() +
#         labs(title = "Variable Importance in the Boosted Model",
#              subtitle = "Based on model's internal importance metric")
# print(vip_plot)

# Create an explainer for Shapley values
# First, get the preprocessed training data
# train_processed <- brt_fit %>%
#         extract_recipe() %>%
#         prep(retain = TRUE) %>%
#         bake(new_data = train_data) %>%
#         select(-quality_level, -typicality_value) # Remove response variable and weights

# We will use the DALEX package to compute shapely values. 
# This package cannot deal with the importance weight column, as it has an unusual object class. 
# This would lead to an error if we do not omit this column here. 
# train_processed <- bake(prepped_recipe, new_data = train_data) %>% select(-typicality_value)
# test_processed  <- bake(prepped_recipe, new_data = test_data) %>%  select(-typicality_value)


# # Create explainer
# explainer <- DALEX::explain(
#         model = brt_fit,
#         data = train_processed,
#         y = as.numeric(train_data$quality_level) - 1, # Convert ordered factor to 0-based numeric
#         type = "regression",
#         label = "XGBoost"
# )




# shapely values ---------------------------------------------------------
# Currently of questionable additional value. 

# # # Calculate Shapley values for a set of observations
# # For demonstration, we'll use a subset of test observations
# set.seed(456)
# sample_rows <- sample(1:nrow(test_data), 20) # Sample 20 observations
# sample_data <- test_data[sample_rows, ]

# # Calculate Shapley values for each observation in the sample
# shapley_values <- list()
# for (i in 1:length(sample_rows)) {
#         # Process the observation using the same recipe
#         obs_processed <- bake(prepped_recipe, new_data = sample_data[i,])

#         # obs_processed <- brt_fit %>%
#         #         extract_recipe() %>%
#         #         prep() %>%
#         #         bake(new_data = sample_data[i, ])
        
#         # Calculate Shapley values
#         shap_i <- DALEX::predict_parts(
#                 explainer,
#                 new_observation = obs_processed,
#                 type = "shap",
#                 B = 50 # Number of samples for approximation
#         )
        
#         # Add typicality and scenario information
#         shap_i$typicality  <- sample_data$typicality_value[i]
#         shap_i$scenario_id <- sample_data$scenario_id[i]
        
#         shapley_values[[i]] <- shap_i
# }

# # Combine all Shapley values
# all_shapley <- bind_rows(shapley_values)

# # Calculate average absolute Shapley values (weighted by typicality)
# weighted_shapley <- all_shapley %>%
#         group_by(variable) %>%
#         mutate(typicality = as.numeric(typicality))%>%
#         summarize(
#                 mean_abs_contribution = weighted.mean(abs(contribution), typicality),
#                 sd_contribution = sd(contribution)
#         ) %>%
#         arrange(desc(mean_abs_contribution))

# # Plot weighted Shapley importance
# ggplot(weighted_shapley, aes(x = reorder(variable, mean_abs_contribution), 
#                              y = mean_abs_contribution)) +
#         geom_bar(stat = "identity", fill = "steelblue") +
#         geom_errorbar(aes(ymin = mean_abs_contribution - sd_contribution,
#                           ymax = mean_abs_contribution + sd_contribution), 
#                       width = 0.2) +
#         coord_flip() +
#         theme_minimal() +
#         labs(
#                 title = "Feature Importance Based on Shapley Values",
#                 subtitle = "Weighted by Scenario Typicality",
#                 x = "Feature",
#                 y = "Mean |Shapley Value|"
#         )

# # Calculate expected ranges for metrics by quality level (weighted by typicality)
# ranges_by_quality <- data %>%
#         group_by(quality_level) %>%
#         summarize(
#                 across(
#                         starts_with("metric"),
#                         list(
#                                 weighted_low = ~weighted.mean(quantile(., 0.05), w = typicality_value),
#                                 weighted_high = ~weighted.mean(quantile(., 0.95), w = typicality_value)
#                         )
#                 )
#         )

# print(ranges_by_quality)

# # Create partial dependence plots for the top metrics
# top_metrics <- weighted_shapley %>%
#         filter(grepl("metric", variable)) %>%
#         head(2) %>%
#         pull(variable)

# for (metric in top_metrics) {
#         # Create a pdp for the metric
#         pdp <- DALEX::model_profile(
#                 explainer,
#                 variables = metric,
#                 N = 100
#         )
        
#         # Plot the pdp
#         plot(pdp) +
#                 ggtitle(paste("Partial Dependence Plot for", metric)) +
#                 theme_minimal()
# }
