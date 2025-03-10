
# describe ---------------------------------------------------------------
#' @export
calculate_auc <- function(x, y) {
        # Sort x and y together (in case they're not already ordered)
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        
        # Calculate width of each trapezoid
        dx <- diff(x)
        
        # Calculate mean height of each trapezoid
        mean_height <- (y[-1] + y[-length(y)]) / 2
        
        # Calculate area as sum of trapezoid areas
        area <- sum(dx * mean_height)
        
        return(area)
}


# describe ---------------------------------------------------------------
#' @export
find_max_consecutive_sum <- function(x, window_size = 3) {
        # Input validation
        if (!is.numeric(x)) stop("Input must be numeric")
        if (length(x) < window_size) stop("Input length must be >= window_size")
        if (window_size < 1) stop("Window size must be positive")
        
        # Create all possible windows
        n_windows <- length(x) - window_size + 1
        window_sums <- numeric(n_windows)
        window_indices <- list()
        
        # Calculate sum for each window
        for (i in 1:n_windows) {
                current_window <- x[i:(i + window_size - 1)]
                window_sums[i] <- sum(current_window)
                window_indices[[i]] <- i:(i + window_size - 1)
        }
        
        # Find the maximum sum and its location
        max_index <- which.max(window_sums)
        
        # Return results as a list
        return(list(
                max_sum = window_sums[max_index],
                values = x[window_indices[[max_index]]],
                positions = window_indices[[max_index]],
                all_sums = window_sums
        ))
}

# asd --------------------------------------------------------------------
#' @export
render_table <- function(x, variable){
        x2  <- data.table(
                        value                 = x, 
                        simulation            = rep(1:max.ss, times = neval/max.ss),
                        q                     = rep(1:max.q, each  = neval/max.q),
                        cluster_adjustment    = rep(1:within.q, times = max.q),
                        metric                = variable, 
                        variables             = rep(j.n.variables, each = neval/max.q),
                        contraction_points    = j.contraction.points,
                        contraction_centroids = j.contraction.centroids,
                        env_asw               = j.env.asw,
                        importance            = rep(j.importance, each = neval/max.q)
                )
        return(x2)
}



# Assign score based on qunatiles ---------------------------------------------------
# This funciton assigns a score between 1 and 5 for each metric based on the quintile it is in.
#' @export
quintile_score <- function(x) {
        y <- findInterval(x, quantile(x, probs = seq(0, 1, 0.2))) 
        y[which(y == 6)] <- 5
        return(y)
}


#' @export
prop_sample <- function(x, N) {
        # Calculate original proportions
        props <- table(x) / length(x)
        
        # Calculate how many of each number we need
        # Using round to ensure we get whole numbers
        n_each <- round(N * props)
        
        # Make sure we get exactly N samples (might be off by 1-2 due to rounding)
        while(sum(n_each) != N) {
                if(sum(n_each) < N) {
                        # Add 1 to the value with the biggest difference from target
                        target_props <- props * N
                        diff <- target_props - n_each
                        n_each[which.max(diff)] <- n_each[which.max(diff)] + 1
                } else {
                        # Subtract 1 from the value with the smallest difference from target
                        target_props <- props * N
                        diff <- n_each - target_props
                        n_each[which.max(diff)] <- n_each[which.max(diff)] - 1
                }
        }
        # Sample from each group
        result <- numeric(0)
        for(num in names(n_each)) {
                positions <- which(x == as.numeric(num))
                result <- c(result, sample(positions, n_each[num]))
        }
        
        # Return sorted positions
        sort(result)
}


# HMSC ------------------------------------------------------------------------------

#  ——— determine spatial scale 
#' @export
determine_spatial_scale <- function(x){
        o.sf <- st_as_sf(unique(x, by = "original_site_name"),
                         coords = c("x.coord", "y.coord"),
                         crs = 3035)
        o.sf <- drop_units(st_distance(o.sf))
        o.sf <- o.sf[lower.tri(o.sf)]
        o.sf <- list(
                min = min(o.sf),
                max = max(o.sf),
                mean = mean(o.sf),
                median = median(o.sf)
        )
        return(o.sf)
}
#  ——— PSRF Check 
#' @export
gelman_check <- function(posterior){
        o.ge    <- gelman.diag(x = posterior$Beta)
        #- species names
        o.species_names  <- o.ge$psrf %>%
                rownames %>%
                str_extract(pattern = ",\\ .*\\(S") %>%
                str_remove(",") %>%
                str_remove("\\(S") %>%
                str_trim %>%
                unique
        #- Create empty vectors to store results
        o.means <- c()
        #- Calculate statistics for each species
        for (j in seq_along(o.species_names)) {
                # Get rows corresponding to current species
                j.species_rows <- grep(o.species_names[j], rownames(o.ge$psrf))
                # Calculate mean and max of point estimates for these rows
                o.means %<>% append(c(o.ge$psrf[j.species_rows, "Point est."]))
                rm(list = ls()[grepl("^j\\.", x = ls())])
        }
        rm(j)
        return(o.means)
}




# AEM -------------------------------------------------------------------------------
#' @export
load_river_networks <- function(data, base_path = "E:/Arbeit/data/river_network/hydrography90m/network_in_catchments/") {
        catchments <-
                unique(data$ID) %>% 
                str_remove("[0-9].*") %>% 
                unique
        files <- list.files(base_path, full.names = TRUE) %>%
                .[grep(paste(catchments, collapse = "|"), ., ignore.case = TRUE)]
        
        rivers <- map(files, st_read_parquet) %>%
                bind_rows() %>%
                filter(strahler > 2)
        
        return(rivers)
}
#' @export
create_directional_incidence <- function(cost_matrix, sites) {
        # Create adjacency matrix based on flow direction
        n_sites <- nrow(cost_matrix)
        adj_matrix <- matrix(0, n_sites, n_sites)
        
        # Determine flow direction based on elevation
        for(i in 1:n_sites) {
                for(j in 1:n_sites) {
                        if(is.finite(cost_matrix[i,j])) {
                                # If site i is higher than site j, flow goes from i to j
                                if(sites$elevation[i] > sites$elevation[j]) {
                                        adj_matrix[i,j] <- 1
                                }
                        }
                }
        }
        
        # Create directed graph
        g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
        edges <- E(g)
        edge_list <- as_edgelist(g)
        
        # Initialize matrices
        n_links <- length(edges)
        incidence <- matrix(0, nrow = n_sites, ncol = n_links)
        weights <- numeric(n_links)
        
        # Fill matrices considering flow direction
        for(i in seq_len(n_links)) {
                v1 <- edge_list[i,1]  # upstream site
                v2 <- edge_list[i,2]  # downstream site
                
                # Get all vertices reachable downstream from v1
                reachable <- subcomponent(g, v1, mode = "out")
                
                # Set 1s for all downstream vertices in this link's column
                incidence[reachable, i] <- 1
                
                # Weight can include both distance and elevation difference
                dist_weight <- cost_matrix[v1, v2]
                elev_diff <- sites$elevation[v1] - sites$elevation[v2]
                
                # Calculate slope (rise over run)
                slope <- elev_diff / dist_weight
                # Adjust weight inversely with slope
                # Using exponential decay function: weight = distance * exp(-k * slope)
                # k is a scaling factor that determines how quickly weights decrease with slope
                k <- 2  
                weights[i] <- dist_weight * exp(-k * slope)
        }
        
        return(list(incidence = incidence, weights = weights))
}
#' @export
prepare_directional_network <- function(rivers, sites) {
        
        site_buffer <- sites %>%
                st_bbox() %>%
                st_as_sfc() %>%
                st_as_sf() %>%
                st_buffer(100000) %>%
                st_transform(4326)
        # Filter and prepare network
        network <- rivers %>%
                st_filter(site_buffer) %>%
                st_cast("LINESTRING") %>%
                as_sfnetwork(directed = TRUE) %>% 
                activate("edges") %>%
                mutate(
                        weight = edge_length(),
                        # Add flow direction based on elevation difference
                        flow_direction = source_elev - outlet_elev 
                )
        
        return(network)
}

# asd --------------------------------------------------------------------
#' @export
balance_clusters <- function(data, clusters, min_size) {
        # Count members in each cluster
        cluster_counts <- table(clusters)
        small_clusters <- which(cluster_counts < min_size)
        
        # Identify clusters that need more members
        for (small_cluster in small_clusters) {
                needed <- min_size - cluster_counts[small_cluster]
                
                # Find cluster centers
                centers <- stats::aggregate(data, by = list(clusters), FUN = mean)
                
                # For each small cluster, find closest points from large clusters
                large_clusters <- which(cluster_counts > min_size + needed)
                if (length(large_clusters) == 0) next
                
                # Points in large clusters
                large_indices <- which(clusters %in% large_clusters)
                
                # Calculate distances to the small cluster center
                small_center <- centers[small_cluster, -1]
                distances <- apply(data[large_indices, ], 1, function(x) 
                        sqrt(sum((x - small_center)^2)))
                
                # Find the closest points and reassign them
                closest <- large_indices[order(distances)][1:needed]
                clusters[closest] <- small_cluster
        }
        
        return(clusters)
}

# asd --------------------------------------------------------------------
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
#' @export
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
        # return(list(centroids = new_centroids, observations = new_observations))
        return(new_observations)
}


