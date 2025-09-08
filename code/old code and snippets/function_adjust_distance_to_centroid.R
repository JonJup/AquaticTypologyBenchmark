#' Adjust observation distances to centroids
#' 
#' This function moves observations closer to or further from their assigned centroids
#' by a consistent fraction of their current distance.
#' 
#' @param centroids Matrix where each row is a centroid and columns are variables
#' @param observations Matrix where each row is an observation and columns are variables
#' @param cluster_assignments Vector indicating which cluster each observation belongs to
#' @param proximity_factor Numeric value controlling how observations move relative to centroids:
#'        positive = move observations closer to centroids
#'        negative = move observations away from centroids
#'        0 = no change
#' 
#' @return Matrix of adjusted observations
#' 
adjust_observation_proximity <- function(centroids, observations, cluster_assignments, proximity_factor) {
        # Ensure inputs are matrices
        centroids <- as.matrix(centroids)
        observations <- as.matrix(observations)
        
        # Initialize matrix for adjusted observations
        adjusted_observations <- matrix(0, nrow = nrow(observations), ncol = ncol(observations))
        
        # For each observation
        for (i in 1:nrow(observations)) {
                # Get assigned cluster
                cluster_id <- cluster_assignments[i]
                
                # Get centroid for this observation
                assigned_centroid <- centroids[cluster_id, ]
                
                # Calculate vector from centroid to observation
                centroid_to_obs <- observations[i, ] - assigned_centroid
                
                # Adjust the distance: 
                # - proximity_factor > 0: move closer to centroid
                # - proximity_factor < 0: move away from centroid
                adjusted_vector <- centroid_to_obs * (1 - proximity_factor)
                
                # Calculate new position
                adjusted_observations[i, ] <- assigned_centroid + adjusted_vector
        }
        
        return(adjusted_observations)
}
