# ————————————————————————————————————————— #
# ——— Add environmental data to biota   ——— # 
# ————————————————————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 26.11.2024

# load and prepare data -------------------------------------------------------------------------
files      <- dir_ls("data/biota", regexp = "01_")
bio.names  <- sapply(files, function (x) str_remove(x, "data/biota/01_")) %>%  sapply(function(x) str_remove(x, "_w_catchment_id.rds"))
catchments <- dir_ls("data/catchments") 
bio.list <- lapply(files, readRDS)
n_bio_datasets <- length(bio.list)
out_list <- vector("list", n_bio_datasets)
names(out_list) <- paste0("out", 1:n_bio_datasets)

# loop over catchments  ---------------------------------------------------
for (i in 1:length(catchments)) {
        i.vec  <- read_parquet(catchments[i])
        i.data <- lapply(bio.list, function(x) x[ID %in% i.vec$ID])
        out <- lapply(i.data, function(x) merge(x, i.vec, by = "ID", all.x = T))
        # Store results for each biological dataset
        for (j in 1:n_bio_datasets) {
                out_list[[j]][[i]] <- out[[j]]
        }
        rm(list = ls()[grepl("^i\\.", x = ls())])
        rm(i)
}
# Combine results
final_results <- lapply(out_list, function(x) rbindlist(x, fill = T))

final_results %<>% lapply(function(x) x[abundance != 0])
final_results %<>% lapply(function(x) x[, PA := 1])

# save to file ----------------------------------------------------------------------
lapply(1:n_bio_datasets,
       function(x)
               saveRDS(
                        final_results[[x]], 
                        paste0("data/biota/02_",bio.names[x],"_w_environment.rds")
               )
       )

