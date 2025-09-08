# Load required packages
library(mlr3)               # Core mlr3 package
library(mlr3learners)       # For ML algorithm implementations
#library(mlr3extralearners)  # For xgboost 
library(mlr3tuning)         # For hyperparameter tuning
library(paradox)            # For search spaces
library(iml)                # For Shapley values
library(data.table)         # For data manipulation
library(ggplot2)            # For visualization

# Create example data (replace with your actual data)
set.seed(123)
n <- 1000

data <- data.table(
  simulation_id = 1:n,
  scenario_id = rep(1:(n/10), each = 10),
  spatial_scale = factor(sample(c("local", "regional", "continental"), n, replace = TRUE)),
  taxonomic_group = factor(sample(c("plants", "invertebrates", "vertebrates"), n, replace = TRUE)),
  metric_1 = rnorm(n),
  metric_2 = rnorm(n),
  metric_3 = rnorm(n),
  metric_4 = rnorm(n),
  quality_level = factor(sample(c("poor", "mediocre", "good"), n, replace = TRUE),
                         levels = c("poor", "mediocre", "good"))
)

# Add typicality values at the scenario level
scenario_typicality <- data.table(
  scenario_id = 1:(n/10),
  typicality_value = rbeta(n/10, 2, 2)  # Beta distribution for values between 0-1
)

# Join typicality values back to main data
data <- merge(data, scenario_typicality, by = "scenario_id")

# Define the task
task <- as_task_classif(data, 
                        target = "quality_level",
                        positive = "good")

# Set observation weights based on typicality
task$set_weights(data$typicality_value)

# Create a learner with hyperparameters
learner <- lrn("classif.xgboost",
               predict_type = "prob",
               nrounds = 500,
               eta = 0.01,
               max_depth = 3,
               min_child_weight = 5,
               subsample = 0.8,
               colsample_bytree = 0.8)

# Split data for training and testing
splits <- partition(task, ratio = 0.8, stratify = TRUE)
task_train <- task$clone()$filter(splits$train)
task_test <- task$clone()$filter(splits$test)

# Train the model
learner$train(task_train)

# Make predictions on test set
predictions <- learner$predict(task_test)

# Calculate performance
measure <- msr("classif.acc")
performance <- predictions$score(measure)
print(paste("Test accuracy:", performance))

# Examine feature importance from XGBoost directly
importance <- learner$model$importance()
importance_dt <- data.table(
  Feature = names(importance),
  Importance = importance
)
setorder(importance_dt, -Importance)
print(importance_dt)

# Plot feature importance
ggplot(importance_dt[1:10], aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "XGBoost Feature Importance",
       x = "Feature",
       y = "Importance")

# Calculate Shapley values with iml
# Prepare test data without target variable
X_test <- task_test$data(cols = task_test$feature_names)

# Select a subset for Shapley (it's computationally intensive)
set.seed(456)
sample_rows <- sample(1:nrow(X_test), 20)
X_sample <- X_test[sample_rows]
weights_sample <- task_test$weights[sample_rows]

# Create a Predictor object
predictor <- Predictor$new(
  model = learner,
  data = X_test,
  y = as.numeric(task_test$truth()) - 1, # Convert factor levels to 0-based
  type = "prob"
)

# Calculate Shapley values
shapley <- Shapley$new(
  predictor = predictor,
  x.interest = X_sample,
  sample.size = 50  # Number of Monte Carlo samples
)

# Get Shapley results
shap_results <- shapley$results
shap_results$typicality <- weights_sample[match(shap_results$.id, 1:length(weights_sample))]

# Calculate weighted mean absolute Shapley values
shap_summary <- shap_results[, .(
  mean_abs_value = weighted.mean(abs(phi), w = typicality),
  sd_value = sd(phi)
), by = "feature"]
setorder(shap_summary, -mean_abs_value)

# Plot Shapley importance
ggplot(shap_summary, aes(x = reorder(feature, mean_abs_value), y = mean_abs_value)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_abs_value - sd_value, 
                    ymax = mean_abs_value + sd_value),
                width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Feature Importance Based on Shapley Values",
    subtitle = "Weighted by Scenario Typicality",
    x = "Feature",
    y = "Mean |Shapley Value|"
  )

# Calculate expected ranges for metrics by quality level
ranges_by_quality <- data[, .(
  weighted_low = weighted.mean(quantile(get(metric), 0.05), w = typicality_value),
  weighted_high = weighted.mean(quantile(get(metric), 0.95), w = typicality_value)
), by = .(quality_level, metric = c("metric_1", "metric_2", "metric_3", "metric_4"))]

print(ranges_by_quality)

# Create partial dependence plots for top metrics
effect_plots <- list()
top_metrics <- shap_summary[grepl("metric", feature)][1:2, feature]

for (metric in top_metrics) {
  # Create a FeatureEffect object
  effect <- FeatureEffect$new(
    predictor = predictor,
    feature = metric,
    method = "pdp"  # Partial Dependence Plot
  )
  
  # Plot the effect
  effect_plots[[metric]] <- plot(effect) +
    ggtitle(paste("Effect of", metric, "on Quality Level")) +
    theme_minimal()
  
  print(effect_plots[[metric]])
}