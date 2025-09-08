# setup -----------------------------------------------------------------------------
library(arrow)
library(data.table)
library(dplyr)
library(sf)
library(sfarrow)

# load data -------------------------------------------------------------------------
data     <- readRDS("E:/Arbeit/data/biota/invertebrate data/00_combine_invertebrates/combined_invertebrate_data.rds")
vector_data <- fs::dir_ls("E:/Arbeit/data/river_network/eu_hydro_dem/parquet/")


# prepare sites ----------------------------------------------------------------
sites <- unique(data, by = c("site_id"))
sites <- sites[!is.na(x.coord)]
sites <- st_as_sf(sites, coords = c("x.coord", "y.coord"), crs = 3035)
sites <- st_transform(sites, 4326)

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
         
        
}

out.ls2 <- bind_rows(out.ls)
sites2  <- left_join(sites, out.ls2, by = c("comb_site_id"))
sites3 <- select(sites2, comb_site_id, ID) %>% st_drop_geometry %>% setDT
data2 <- sites3[data, on = "comb_site_id"]

saveRDS(data2, "data/invertebrates_w_catchment_id.rds")
