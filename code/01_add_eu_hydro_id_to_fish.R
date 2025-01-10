
# setup -----------------------------------------------------------------------------
library(arrow)
library(data.table)
library(dplyr)
library(sf)
library(sfarrow)

# load data -------------------------------------------------------------------------
fish <- read_parquet("E:/Arbeit/data/biota/fish data/complete_data.parquet")
vector_data <- fs::dir_ls("E:/Arbeit/data/river_network/eu_hydro_dem/parquet/")


# prepare fish sites ----------------------------------------------------------------
sites <- unique(fish, by = c("comb_site_id"))
sites <- st_as_sf(sites, coords = c("x.coord", "y.coord"), crs = 3035)
sites <- st_transform(sites, 4326)

out.ls <- vector(mode = "list", length = length(vector_data))
for (i in 1:length(vector_data)){
        i.vec  <- st_read_parquet(vector_data[i])
        if (i == 14) i.vec <- st_make_valid(i.vec)
        i.join <- st_join(sites, i.vec)
        
        out.ls[[i]] <- 
                i.join %>% 
                filter(!is.na(ID)) %>%
                select(c("comb_site_id", "ID")) %>%
                st_drop_geometry()
         
        
};beepr::beep()

out.ls2 <- bind_rows(out.ls)
sites2  <- left_join(sites, out.ls2, by = c("comb_site_id"))
sites3 <- select(sites2, comb_site_id, ID) %>% st_drop_geometry %>% setDT
fish2 <- sites3[fish, on = "comb_site_id"]

write_parquet(fish2, "data/fish_w_catchment_id.parquet")
