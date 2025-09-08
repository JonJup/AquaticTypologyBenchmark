load_river_networks <- function(data, base_path = "E:/Arbeit/data/river_network/hydrography90m/network_in_catchments/") {
        catchments <-
                unique(data$ID) %>% 
                str_remove("[0-9].*") %>% 
                unique
        files <- list.files(base_path, full.names = TRUE) %>%
                .[grep(paste(catchments, collapse = "|"), ., ignore.case = TRUE)]
        
        rivers <- map(files, st_read_parquet) %>%
                bind_rows() %>%
                filter(strahler > 2)
        
        return(rivers)
}