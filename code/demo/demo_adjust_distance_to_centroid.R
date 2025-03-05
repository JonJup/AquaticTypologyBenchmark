# Generate example data
set.seed(123)
n_per_cluster <- 50
n_clusters <- 4
n_dims <- 2  # Using 2 dimensions for easier visualization

# Create cluster centroids
centroids <- matrix(c(
        1, 1,    # Centroid 1
        -1, 1,   # Centroid 2
        -1, -1,  # Centroid 3
        1, -1    # Centroid 4
), nrow = n_clusters, byrow = TRUE)

# Create observations around each centroid
observations <- matrix(0, nrow = n_clusters * n_per_cluster, ncol = n_dims)
cluster_assignments <- rep(1:n_clusters, each = n_per_cluster)

# Generate observations with varying distances around centroids
for (i in 1:nrow(observations)) {
        cluster_id <- cluster_assignments[i]
        observations[i, ] <- centroids[cluster_id, ] + rnorm(n_dims, 0, 0.3)
}

# Standardize all variables
observations <- scale(observations)
centroids <- scale(centroids)

# Create visualizations of different proximity adjustments
par(mfrow = c(2, 2))

# Original clusters
plot(observations, col = cluster_assignments, pch = 16, 
     main = "Original Clusters", xlab = "Variable 1", ylab = "Variable 2")
points(centroids, col = 1:4, pch = 8, cex = 2)

# Make clusters more compact (move observations closer to centroids)
closer_observations <- adjust_observation_proximity(centroids, observations, 
                                                    cluster_assignments, 
                                                    proximity_factor = 0.5)
plot(closer_observations, col = cluster_assignments, pch = 16, 
     main = "More Compact Clusters (proximity = 0.5)", 
     xlab = "Variable 1", ylab = "Variable 2")
points(centroids, col = 1:4, pch = 8, cex = 2)

# Make clusters very compact
very_compact_observations <- adjust_observation_proximity(centroids, observations, 
                                                          cluster_assignments, 
                                                          proximity_factor = 0.8)
plot(very_compact_observations, col = cluster_assignments, pch = 16, 
     main = "Very Compact Clusters (proximity = 0.8)", 
     xlab = "Variable 1", ylab = "Variable 2")
points(centroids, col = 1:4, pch = 8, cex = 2)

# Make clusters more dispersed (move observations away from centroids)
dispersed_observations <- adjust_observation_proximity(centroids, observations, 
                                                       cluster_assignments, 
                                                       proximity_factor = -0.5)
plot(dispersed_observations, col = cluster_assignments, pch = 16, 
     main = "More Dispersed Clusters (proximity = -0.5)", 
     xlab = "Variable 1", ylab = "Variable 2")
points(centroids, col = 1:4, pch = 8, cex = 2)
