#' Adjust cluster centroids and observation distances
#' 
#' This function performs two adjustments:
#' 1. Moves cluster centroids further apart or closer together
#' 2. Adjusts how close observations are to their assigned centroids
#' 
#' @param centroids Matrix where each row is a centroid and columns are variables
#' @param observations Matrix where each row is an observation and columns are variables
#' @param cluster_assignments Vector indicating which cluster each observation belongs to
#' @param centroid_adjustment_factor Numeric controlling how centroids move relative to each other:
#'        positive = move centroids further apart, negative = move centroids closer together
#' @param compactness_factor Numeric controlling how observations move relative to their centroids:
#'        positive = move observations closer to centroids, negative = move observations away from centroids
#' 
#' @return List containing the adjusted centroids and observations
#' 
adjust_clusters <- function(centroids, observations, cluster_assignments, 
                            centroid_adjustment_factor = 0, compactness_factor = 0) {
        
        # Step 1: Calculate the grand centroid (center of all centroids)
        grand_centroid <- colMeans(centroids)
        
        # Step 2: Move all centroids relative to grand centroid
        scaling_factor <- 1 + centroid_adjustment_factor
        new_centroids <- grand_centroid + scaling_factor * (centroids - grand_centroid)
        
        # Step 3: Move observations with their centroids and adjust their distances to centroids
        new_observations <- matrix(0, nrow = nrow(observations), ncol = ncol(observations))
        for (i in 1:nrow(observations)) {
                cluster_id <- cluster_assignments[i]
                
                # Get original and new centroid for this observation
                orig_centroid <- centroids[cluster_id, ]
                new_centroid <- new_centroids[cluster_id, ]
                
                # Calculate the vector from original centroid to observation
                centroid_to_obs_vector <- observations[i, ] - orig_centroid
                
                # Apply compactness adjustment: scale the distance from centroid to observation
                # compactness_factor > 0: observations move closer to centroids
                # compactness_factor < 0: observations move further from centroids
                adjusted_vector <- centroid_to_obs_vector * (1 - compactness_factor)
                
                # Place observation relative to new centroid position
                new_observations[i, ] <- new_centroid + adjusted_vector
        }
        colnames(new_observations) <- colnames(observations)
        #return(list(centroids = new_centroids, observations = new_observations))
        return(new_observations)
}
