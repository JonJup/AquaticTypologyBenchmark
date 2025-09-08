# ————————————————————————————————————————— #
# ——— Add environmental data to diatoms ——— # 
# ————————————————————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 26.11.2024

# load data -------------------------------------------------------------------------
data       <- readRDS("data/diatoms_w_catchment_id.rds")
catchments <- fs::dir_ls("data/catchments") 
out        <- list()
for (i in 1:length(catchments)) {
        i.vec  <- read_parquet(catchments[i])
        i.data <- data[ID %in% i.vec$ID]
        out[[i]]  <- merge(i.data, i.vec, by = "ID", all.x = T)
        rm(list = ls()[grepl("^i\\.", x = ls())])
        rm(i)
}

out2 <- rbindlist(out, fill = T)
#- write 
out2 <- out2[abundance != 0]
out2[, PA := 1]

# save to file ----------------------------------------------------------------------

saveRDS(out2, "data/diatoms_w_env.rds")
rm(out, out2, data, catchments); gc()
