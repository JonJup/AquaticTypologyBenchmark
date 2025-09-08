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