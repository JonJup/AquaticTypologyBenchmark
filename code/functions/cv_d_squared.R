cv_d_squared <- function(species_data, fuzzy_memberships, k=5, family="binomial") {
        set.seed(123)
        n <- nrow(species_data)
        folds <- sample(1:k, n, replace=TRUE)
        
        # Store test set performance
        fold_d2 <- numeric(k)
        
        for(cv in 1:k) {
          # Split data
          train_idx <- folds != cv
          test_idx  <- folds == cv
          
          # Create training and test datasets
          train_species <- species_data[train_idx,]
          test_species  <- species_data[test_idx,]
          
          # Prepare fuzzy membership data frames with proper column names
          fuzzy_df <- as.data.frame(fuzzy_memberships)
          # Ensure consistent naming. This is important for the 
          # model formula later. 
          colnames(fuzzy_df) <- paste0("M", 1:ncol(fuzzy_df))  
          
          # split membership data into training and test data
          train_fuzzy <- fuzzy_df[train_idx,]
          test_fuzzy  <- fuzzy_df[test_idx,]
          
          # Convert to mvabund objects
          train_mv <- mvabund(train_species)
          test_mv  <- mvabund(test_species)
          
          # Fit null model
          null_model <- manyglm(train_mv ~ 1, family=family)
          
          # Calculate null deviance for test data
          test_null_dev <- sum((test_species - colMeans(train_species))^2)
          
          # Get the column names as a string for formula construction
          predictor_cols <- paste(colnames(train_fuzzy), collapse="+")
          
          # Create and evaluate formula in the model environment
          formula_str  <- paste("train_mv ~", predictor_cols)
          full_formula <- formula(formula_str)
          
          # Fit the model
          full_model <- try(manyglm(full_formula, data=train_fuzzy, family=family), silent=TRUE)
          
          if(inherits(full_model, "try-error")) {
            warning(paste("Error in fold", cv, "- skipping"))
            fold_d2[cv] <- NA
            next
          }
          
          # Calculate predictions using a safer approach
          # Instead of using predict.manyglm directly, we'll manually calculate
          # predictions using the model coefficients
          fitted_values <- try(predict.manyglm(object = full_model, 
                                               newdata = test_fuzzy, 
                                               type = "response"))
          
          if(inherits(fitted_values, "try-error")) {
                  # Get coefficients
                  beta <- coef(full_model)
                  
                  # Prepare test design matrix (ensuring column order matches coefficients)
                  X_test <- model.matrix(~ ., data=test_fuzzy)
                  # Ensure column order matches coefficient order
                  X_test <- X_test[, rownames(beta), drop=FALSE]
                  # Calculate linear predictor for each species
                  fitted_values <- matrix(0, nrow=nrow(X_test), ncol=ncol(beta))
                  for(j in 1:ncol(beta)) {
                          fitted_values[,j] <- X_test %*% beta[,j, drop=FALSE]
                  }
                  # Convert from linear predictor to response scale
                  if(family == "negative.binomial" || family == "poisson") {
                          fitted_values <- exp(fitted_values)
                  } else if(family == "binomial") {
                          fitted_values <- 1/(1 + exp(-fitted_values))
                  }
          }
          # Calculate squared residuals
          test_resid_dev <- sum((test_species - fitted_values)^2)
          
          # Calculate D² for this fold
          fold_d2[cv] <- 1 - test_resid_dev/test_null_dev
        }
        
        # Return average D² across folds (excluding NAs)
        return(mean(fold_d2, na.rm=TRUE))
}
