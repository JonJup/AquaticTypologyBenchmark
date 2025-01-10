# —————————————————— #
# ——— Add environmental data to fish ——— # 
# —————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 26.11.2024


# load data -------------------------------------------------------------------------
data       <- read_parquet("data/fish_w_catchment_id.parquet")
catchments <- dir_ls("data/catchments") 
out        <- list()
for (i in 1:length(catchments)) {
        i.vec  <- read_parquet(catchments[i])
        if (str_detect(catchments[i], "Danube")){
                i.vec %<>% rename(upstream_catchment_area = upstream_catchment_area.x)
                i.vec %<>% select(!upstream_catchment_area.y)
        }
        i.data <- data[ID %in% i.vec$ID]
        out[[i]]  <- merge(i.data, i.vec, by = "ID", all.x = T)
        rm(list = ls()[grepl("^i\\.", x = ls())])
        rm(i)
}
out2 <- rbindlist(out, fill = T)

#sapply(out, function(x) any(str_detect(names(x), "upstream_catchment_area.x")))

#- write 
out2 <- out2[abundance != 0]
out2[, PA := 1]

# save to file ----------------------------------------------------------------------

write_parquet(out2, "data/fish_w_env.parquet")
rm(out, out2, data, catchments); gc()
