#' Adjust distances between cluster centroids
#' 
#' This function moves cluster centroids either further apart or closer together based on an
#' adjustment factor, while maintaining the relative distances between observations and their
#' assigned centroids.
#' 
#' @param centroids Matrix where each row is a centroid and columns are variables
#' @param observations Matrix where each row is an observation and columns are variables
#' @param cluster_assignments Vector indicating which cluster each observation belongs to
#' @param adjustment_factor Numeric value controlling movement strength and direction:
#'        positive = move centroids further apart, negative = move centroids closer together
#' 
#' @return List containing the adjusted centroids and observations
#' 
adjust_centroid_locations    <- function(centroids, observations, cluster_assignments, adjustment_factor) {
        
        if(!"matrix" %in% class(observations)){
                observations <- as.data.frame(observations)
                observations <- as.matrix(observations)
        }
        if(!"matrix" %in% class(centroids)){
                centroids <- as.data.frame(centroids)
                centroids <- as.matrix(centroids)
        }
        # Calculate the grand centroid (center of all centroids)
        centroids <- centroids[order(centroids[,1]),]
        centroids <- centroids[, -1]
        grand_centroid <- colMeans(centroids)
        
        # Move all centroids relative to grand centroid
        # adjustment_factor > 0: centroids move further apart
        # adjustment_factor < 0: centroids move closer together
        scaling_factor <- 1 + adjustment_factor
        new_centroids  <- grand_centroid + scaling_factor * (centroids - grand_centroid)
        
        # Move observations with their centroids
        new_observations <- matrix(0, nrow = nrow(observations), ncol = ncol(observations))
        for (iobs in 1:nrow(observations)) {
                cluster_id <- cluster_assignments[iobs]
                # Calculate the displacement of the centroid
                centroid_displacement <- new_centroids[cluster_id, ] - centroids[cluster_id,]
                # Apply the same displacement to the observation
                new_observations[iobs, ] <- observations[iobs, ] + centroid_displacement
        }
        new_centroids <- cbind(1:j.types, new_centroids)
        new_observations <- cbind(cluster_assignments, new_observations)
        colnames(new_centroids)[1] <- "type"
        colnames(new_observations) <- c("type", colnames(observations))
        
        # new_centroids <- setDT(as.data.frame(new_centroids))
        # new_observations <- setDT(as.data.frame(new_observations))
        # 
        # out <- new_centroids[new_observations, on = "type"]
        # return(out)
        
        return(list(centroids = new_centroids, observations = new_observations))
}
