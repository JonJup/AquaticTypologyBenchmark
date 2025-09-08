# load data -------------------------------------------------------------------------

files      <- dir_ls("data/biota", regexp = "00_")
bio.names  <- sapply(files, function (x) str_remove(x, "data/biota/00_combined_")) %>%  sapply(function(x) str_remove(x, ".rds"))
vector_data <- dir_ls("D:/Arbeit/data/river_network/eu_hydro_dem/parquet/")

for (i.biota in bio.names){
        print(paste("STARING", i.biota))
        biota       <- readRDS(paste0("data/biota/00_combined_", i.biota,".rds"))
        
        # prepare sites  ----------------------------------------------------------
        sites <- unique(biota, by = c("comb_site_id")) %>% 
                filter(!is.na(x.coord)) %>%  
                st_as_sf(coords = c("x.coord", "y.coord"), crs = 3035) %>% 
                st_transform(4326)
        out.ls <- vector(mode = "list", length = length(vector_data))
        for (i in 1:length(vector_data)){
                print(i)
                i.vec  <- st_read_parquet(vector_data[i])
                if (i == 13) i.vec <- st_make_valid(i.vec)
                i.join <- st_join(sites, i.vec)
                
                out.ls[[i]] <- 
                        i.join %>% 
                        filter(!is.na(ID)) %>%
                        select(c("comb_site_id", "ID")) %>%
                        st_drop_geometry()
                
                
        }#;beepr::beep()
        
        out.ls <- bind_rows(out.ls)
        sites  <- left_join(sites, out.ls, by = c("comb_site_id"))
        sites  <- select(sites, comb_site_id, ID) %>% st_drop_geometry %>% setDT
        biota   <- sites[biota, on = "comb_site_id"]
        
        saveRDS(biota, paste0("data/biota/01_", i.biota,"_w_catchment_id.rds"))
        rm(out.ls, biota, sites)
}
