#- ID spatial scales 
bio <- read_parquet("data/fish_w_env.parquet")
schemes <- readRDS("data/schemes_fish.rds")
#- loop over schemes
for (o in 1:nrow(schemes)) {
        print(paste(o, "of", nrow(schemes)))
        o.scheme <- schemes[o, ]
        o.data  <- bio[data.set == o.scheme$data.set & year == o.scheme$year]
        o.data  <- o.data[month(date) %in% as.numeric(o.scheme$focal_months[[1]])]
        o.sample.ids <- sample(
                x = unique(o.data$comb_sample_id),
                size = o.scheme$samples,
                replace = F
        )
        o.data <- o.data[comb_sample_id %in% o.sample.ids]
        #- establish spatial scale
        o.sf <- st_as_sf(unique(o.data, by = "original_site_name"),
                         coords = c("x.coord", "y.coord"),
                         crs = 3035)
        o.sf <- drop_units(st_distance(o.sf))
        o.sf <- o.sf[lower.tri(o.sf)]
        # o.sf <- list(
        #         min = min(o.sf),
        #         max = max(o.sf),
        #         mean = mean(o.sf),
        #         median = median(o.sf)
        # )
        schemes[o, `:=` (
                scale_min = min(o.sf), 
                scale_mean = mean(o.sf), 
                scale_median = median(o.sf),
                scale_max = max(o.sf)
                )]
        
}
saveRDS(schemes, "data/schemes_w_scale.rds")
