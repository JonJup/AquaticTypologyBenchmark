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