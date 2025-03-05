
#TODO remove flow direction - I think it doesn't do anything. 

# load the data 
bio <- read_parquet("data/fish_w_env.parquet")
schemes <- readRDS("data/schemes_w_scale.rds")
schemes_small <- schemes[scale_max == min(scale_max)]

bio_small <- bio[data.set == schemes_small$data.set & year == schemes_small$year]
bio_small <- bio_small[month(date) %in% schemes_small$focal_months[[1]]]

#- in which catchments are these sampling locations? 
rivers <- load_river_networks(bio_small)
sites  <- bio_small %>%
        unique(by = "comb_site_id") %>%
        st_as_sf(coords = c("x.coord", "y.coord"), crs = 3035) %>%
        st_transform(4326) 
network <- prepare_directional_network(rivers)
# Calculate network costs
cost_matrix <-
        st_network_cost(network,
                        from = sites,
                        to = sites,
                        weights = "weight") %>%
        drop_units()
matrices <- create_directional_incidence(cost_matrix, sites)

aem_result <- aem(binary.mat = matrices$incidence,
                  weight = matrices$weights)

