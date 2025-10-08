#' Apply constraints to ensure variables stay within valid bounds
#' 
#' @param data Matrix or vector of data to constrain
#' @param variable_names Character vector of variable names
#' 
#' @return Constrained data
apply_constraints <- function(data, variable_names = NULL) {
        # Convert to matrix if needed for uniform handling
        if (is.vector(data)) {
                data <- matrix(data, nrow = 1)
                data <- data.frame(data)
                colnames(data) <- variable_names
        } else if ("data.table" %in% class(data)){
                data <- as.data.frame(data)
        }
        if (is.null(variable_names)) variable_names = colnames(data)
        
        # Identify column indices for different constraint types
        area_cols <- grep("area", variable_names)
        if ("upstream_catchment_area" %in% variable_names){
                upstream_col <- grep("upstream", variable_names)
                area_cols <- area_cols[- which(area_cols == upstream_col)]
        }
        area_special_cols <- grep("area_(calcareous|siliceous|sediment)", variable_names)
        saturated_col <- grep("saturated_soil_water_content", variable_names)
        bioclim06_col <- grep("bioclim06", variable_names)
        positive_cols <- setdiff(1:length(variable_names), bioclim06_col)
        
        # Apply constraints row by row
        for (aci in 1:nrow(data)) {
                # 1. Constrain area variables and saturated_soil_water_content to [0, 1]
                if (length(area_cols) > 0) {
                        data[aci, area_cols] <- pmax(0, pmin(1, data[aci, area_cols]))
                }
                if (length(saturated_col) > 0) {
                        data[aci, saturated_col] <- pmax(0, pmin(1, data[aci, saturated_col]))
                }
                
                # 2. Ensure all variables except bioclim06 stay positive
                if (length(positive_cols) > 0) {
                        data[aci, positive_cols] <- pmax(0, data[aci, positive_cols])
                }
                
                # 3. Normalize area_calcareous, area_siliceous, and area_sediment to sum to 1
                if (length(area_special_cols) == 3) {
                        current_sum <- sum(data[aci, area_special_cols])
                        if (current_sum == 0){
                                data[aci, area_special_cols] <- 1
                                current_sum <- sum(data[aci, area_special_cols])
                        }
                        if (current_sum != 1) {
                                # Normalize to sum to 1
                                data[aci, area_special_cols] <- data[aci, area_special_cols] / current_sum
                        } 
                }
        }
        
        return(data)
}

#' Adjust cluster centroids and observation distances with constraints
#' 
#' This function performs two adjustments while maintaining variable constraints:
#' 1. Moves cluster centroids further apart or closer together
#' 2. Adjusts how close observations are to their assigned centroids
#' 
#' Constraints applied:
#' - Variables with "area_" in name: bounded between 0 and 1
#' - saturated_soil_water_content: bounded between 0 and 1
#' - All variables except bioclim06: must stay positive
#' - area_calcareous + area_siliceous + area_sediment must sum to 1
#' 
#' @param centroids Matrix where each row is a centroid and columns are variables
#' @param observations Matrix where each row is an observation and columns are variables
#' @param cluster_assignments Vector indicating which cluster each observation belongs to
#' @param separation_factor Numeric controlling how centroids move relative to each other:
#'        positive = move centroids further apart, negative = move centroids closer together
#' @param compactness_factor Numeric controlling how observations move relative to their centroids:
#'        positive = move observations closer to centroids, negative = move observations away from centroids
#' @param apply_constraints_flag Logical indicating whether to apply constraints (default TRUE)
#' 
#' @return List containing the adjusted centroids and observations
#' 
#' @export
adjust_clusters <- function(centroids, observations, cluster_assignments,
                            separation_factor = 0, compactness_factor = 0,
                            apply_constraints_flag = TRUE) {
        
        # Get variable names for constraint checking
        variable_names <- colnames(observations)
        if (is.null(variable_names)) {
                variable_names <- colnames(centroids)
        }
        if (is.null(variable_names)) {
                warning("No column names found. Constraints cannot be properly applied.")
                apply_constraints_flag <- FALSE
        }
        
        # Step 1: Calculate the grand centroid (center of all centroids)
        grand_centroid <- colMeans(centroids)
        
        # Step 2: Move all centroids relative to grand centroid
        scaling_factor <- 1 + separation_factor
        new_centroids <- grand_centroid + scaling_factor * (centroids - grand_centroid)
        
        # Apply constraints to new centroids if flag is set
        if (apply_constraints_flag) {
                new_centroids <- apply_constraints(new_centroids)
        }
        
        # Step 3: Move observations with their centroids and adjust their distances to centroids
        new_observations <- matrix(0, nrow = nrow(observations), ncol = ncol(observations))
        for (löö in 1:nrow(observations)) {
                cluster_id <- cluster_assignments[löö]
                
                # Get original and new centroid for this observation
                orig_centroid <- centroids[cluster_id, ]
                new_centroid  <- new_centroids[cluster_id, ]
                
                # Calculate the vector from original centroid to observation
                centroid_to_obs_vector <- observations[löö, ] - orig_centroid
                
                # Apply compactness adjustment: scale the distance from centroid to observation
                # compactness_factor > 0: observations move closer to centroids
                # compactness_factor < 0: observations move further from centroids
                adjusted_vector <- centroid_to_obs_vector * (1 - compactness_factor)
                
                # Place observation relative to new centroid position
                new_observations[löö, ] <- new_centroid + adjusted_vector
        }
        
        # Apply constraints to new observations if flag is set
        if (apply_constraints_flag) {
                new_observations <- apply_constraints(new_observations, variable_names)
        }
        
        # Preserve column names
        colnames(new_observations) <- colnames(observations)
        colnames(new_centroids) <- colnames(centroids)
        
        # Return adjusted observations (following original function signature)
        return(list(observations = new_observations, centroids = new_centroids))
}

#' Validate that constraints are satisfied
#' 
#' This function checks if the data satisfies all constraints and returns
#' a report of any violations.
#' 
#' @param data Matrix of data to validate
#' @param variable_names Character vector of variable names
#' 
#' @return List with validation results
#' 
#' @export
validate_constraints <- function(data, variable_names = NULL) {
        if (is.vector(data)) {
                data <- matrix(data, nrow = 1)
        }
        if (is.null(variable_names)) variable_names = colnames(data)
        violations <- list()
        
        # Check area variables and saturated_soil_water_content bounds
        area_cols <- grep("area", variable_names)
        if ("upstream_catchment_area" %in% variable_names){
                upstream_col <- grep("upstream", variable_names)
                area_cols <- area_cols[- which(area_cols == upstream_col)]
        }
        saturated_col <- grep("saturated_soil_water_content", variable_names)
        
        if (length(area_cols) > 0) {
                out_of_bounds <- which(data[, area_cols] < 0 | data[, area_cols] > 1, arr.ind = TRUE)
                if (nrow(out_of_bounds) > 0) {
                        violations$area_bounds <- out_of_bounds
                }
        }
        
        if (length(saturated_col) > 0) {
                out_of_bounds <- which(data[, saturated_col] < 0 | data[, saturated_col] > 1)
                if (length(out_of_bounds) > 0) {
                        violations$saturated_bounds <- out_of_bounds
                }
        }
        
        # Check positivity constraint
        bioclim06_col <- grep("bioclim06", variable_names)
        positive_cols <- setdiff(1:length(variable_names), bioclim06_col)
        
        if (length(positive_cols) > 0) {
                negative_vals <- which(data[, positive_cols] < 0, arr.ind = TRUE)
                if (nrow(negative_vals) > 0) {
                        violations$negative_values <- negative_vals
                }
        }
        
        # Check sum constraint for special area variables
        area_special_cols <- grep("area_(calcareous|siliceous|sediment)", variable_names)
        if (length(area_special_cols) == 3) {
                sums <- rowSums(data[, area_special_cols])
                not_one <- which(abs(sums - 1) > 1e-10)
                if (length(not_one) > 0) {
                        violations$sum_not_one <- list(rows = not_one, sums = sums[not_one])
                }
        }
        
        is_valid <- length(violations) == 0
        
        return(list(
                is_valid = is_valid,
                violations = violations
        ))
}
#' Calculate effective adjustment factors after constraints
#' 
#' This function calculates what the actual adjustment factors were after
#' constraints have been applied, comparing the constrained positions to
#' the original positions.
#' 
#' @param original_centroids Original centroid positions
#' @param constrained_centroids Constrained centroid positions
#' @param original_observations Original observation positions
#' @param constrained_observations Constrained observation positions
#' @param cluster_assignments Vector indicating which cluster each observation belongs to
#' 
#' @return List containing effective_centroid_factor and effective_compactness_factor
#' 
#' @export
calculate_effective_factors <- function(original_centroids, constrained_centroids,
                                        original_observations, constrained_observations,
                                        cluster_assignments) {
        
        # Calculate effective centroid adjustment factor
        grand_centroid_orig <- colMeans(original_centroids)
        
        # Calculate average scaling for centroids
        centroid_distances_orig <- sqrt(rowSums((original_centroids - grand_centroid_orig)^2))
        centroid_distances_new <- sqrt(rowSums((constrained_centroids - grand_centroid_orig)^2))
        
        # Avoid division by zero
        valid_indices <- centroid_distances_orig > 1e-10
        if (sum(valid_indices) > 0) {
                scaling_ratios <- centroid_distances_new[valid_indices] / centroid_distances_orig[valid_indices]
                effective_centroid_factor <- mean(scaling_ratios) - 1
        } else {
                effective_centroid_factor <- 0
        }
        
        # Calculate effective compactness factor
        compactness_ratios <- numeric()
        
        for (i in 1:nrow(original_observations)) {
                cluster_id <- cluster_assignments[i]
                
                # Original distance to centroid
                orig_dist <- sqrt(sum((original_observations[i, ] - original_centroids[cluster_id, ])^2))
                
                # New distance to constrained centroid
                new_dist <- sqrt(sum((constrained_observations[i, ] - constrained_centroids[cluster_id, ])^2))
                
                # Only consider points that were not at the centroid
                if (orig_dist > 1e-10) {
                        compactness_ratios <- c(compactness_ratios, new_dist / orig_dist)
                }
        }
        
        if (length(compactness_ratios) > 0) {
                # The compactness factor represents how much closer (positive) or further (negative) points are
                # A ratio of 0.8 means points are 20% closer, corresponding to compactness_factor of 0.2
                effective_compactness_factor <- 1 - mean(compactness_ratios)
        } else {
                effective_compactness_factor <- 0
        }
        
        return(list(
                effective_centroid_factor = effective_centroid_factor,
                effective_compactness_factor = effective_compactness_factor
        ))
}
#' Iteratively adjust factors to achieve target effective factors
#' 
#' This function iteratively adjusts the input factors to achieve desired
#' effective factors after constraints are applied. This compensates for
#' the reduction in movement caused by constraints.
#' 
#' @param centroids Matrix where each row is a centroid and columns are variables
#' @param observations Matrix where each row is an observation and columns are variables
#' @param cluster_assignments Vector indicating which cluster each observation belongs to
#' @param target_separation_factor Target effective centroid adjustment factor
#' @param target_compactness_factor Target effective compactness factor
#' @param max_iterations Maximum number of iterations to try (default 10)
#' @param tolerance Tolerance for considering factors close enough (default 0.05)
#' @param verbose Logical indicating whether to print progress (default FALSE)
#' 
#' @return List containing the adjusted observations, centroids, and achieved factors
#' 
#' @export
adjust_clusters_iterative <- function(centroids, observations, cluster_assignments,
                                      target_separation_factor = 0, 
                                      target_compactness_factor = 0,
                                      max_iterations = 10,
                                      tolerance = 0.05,
                                      verbose = FALSE) {
        
        # Initialize with the target factors
        current_separation_input <- target_separation_factor
        current_compactness_input <- target_compactness_factor
        
        # Learning rates for adjustment
        centroid_learning_rate <- 1.5
        compactness_learning_rate <- 1.5
        
        for (iter in 1:max_iterations) {
                # Apply adjustment with current input factors
                result <- adjust_clusters(
                        centroids = centroids,
                        observations = observations,
                        cluster_assignments = cluster_assignments,
                        separation_factor = current_separation_input,
                        compactness_factor = current_compactness_input,
                        apply_constraints_flag = TRUE
                )
                result2 <- calculate_effective_factors(original_centroids   = centroids, 
                                                       constrained_centroids = result$centroids,
                                                       original_observations = observations, 
                                                       constrained_observations = result$observations,
                                                       cluster_assignments = cluster_assignments)
                # Check if we're close enough to targets
                centroid_error <- abs(result2$effective_centroid_factor - target_separation_factor)
                compactness_error <- abs(result2$effective_compactness_factor - target_compactness_factor)
                if (verbose) {
                        cat(sprintf("Iteration %d:\n", iter))
                        cat(sprintf("  Input factors: centroid=%.3f, compactness=%.3f\n", 
                                    current_separation_input, current_compactness_input))
                        cat(sprintf("  Effective factors: centroid=%.3f, compactness=%.3f\n",
                                    result$effective_centroid_factor, result$effective_compactness_factor))
                        cat(sprintf("  Errors: centroid=%.3f, compactness=%.3f\n", 
                                    centroid_error, compactness_error))
                }
                
                # Check convergence
                if (centroid_error < tolerance && compactness_error < tolerance) {
                        if (verbose) {
                                cat("Converged!\n")
                        }
                        result$iterations_used <- iter
                        result$converged <- TRUE
                        return(result)
                }
                
                # Adjust input factors based on the error
                # If effective < target, we need to increase the input
                if (abs(target_separation_factor) > 1e-10) {
                        centroid_ratio <- target_separation_factor / 
                                (result2$effective_centroid_factor + 1e-10)
                        current_separation_input <- current_separation_input + 
                                centroid_learning_rate * (centroid_ratio - 1) * current_separation_input
                }
                
                if (abs(target_compactness_factor) > 1e-10) {
                        compactness_ratio <- target_compactness_factor / 
                                (result2$effective_compactness_factor + 1e-10)
                        current_compactness_input <- current_compactness_input + 
                                compactness_learning_rate * (compactness_ratio - 1) * current_compactness_input
                }
                
                # Decay learning rate
                centroid_learning_rate <- centroid_learning_rate * 0.9
                compactness_learning_rate <- compactness_learning_rate * 0.9
        }
        
        # If we didn't converge, return the last result
        if (verbose) {
                cat("Maximum iterations reached without convergence.\n")
        }
        result$iterations_used <- max_iterations
        result$converged <- FALSE
        return(result)
}

