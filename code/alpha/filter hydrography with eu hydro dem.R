#- catchments and h18v02
library(sf)
library(sfarrow)
library(fs)
library(stringr)
library(dplyr)

catchments <- dir_ls("E://Arbeit/Data/river_network/eu_hydro_dem/parquet/")
rivers     <- dir_ls("E://Arbeit/Data/river_network/hydrography90m/parquet/")
#@8 - 25

#- loop over rivers 
for (i in 8:length(rivers)){
        i.river <- st_read_parquet(rivers[i])
        i.name <- regmatches(rivers[i], regexpr("h\\d+v\\d+", rivers[i]))
        #- loop over catchments 
        for (j in 1:length(catchments)){
                if (i == 8 & j < 25) next()
                print(paste(i, j))
                j.catchment <- st_read_parquet(catchments[j])
                j.catchment <- st_make_valid(j.catchment)
                j.name <- str_remove(catchments[j], "\\.parquet")
                j.name <- str_remove(j.name, ".*/parquet/")
                i.x <- st_filter(i.river, j.catchment)
                if (nrow(i.x) == 0) next()
                st_write_parquet(i.x, paste0("E:Arbeit/Data/river_network/hydrography90m/network_in_catchments/", i.name, "_", j.name, ".parquet"))
                rm(j.catchment)
                rm(i.x)
                rm(j.name)
                gc()
        }
        rm(i.river)
        gc()
}

#- check for multiples. 
made.files <- dir_ls("E://Arbeit/data/river_network/hydrography90m/network_in_catchments/")


#- loop over all catchmens. For each catchment identify all rivernetworks that belong to
#- this catchment. If there are more than one, join them, save the joined object to file, 
#- remove the old files. If there is only one, rename to only the name of the catchment 
#- without reference to the hydropgraphy90m layer. 

for (i in 4:length(catchments)){
        print(paste(i, "/", length(catchments)))
        i.cat <- catchments[i]
        i.name <- 
                str_remove(i.cat, "\\.parquet") %>%  
                str_remove(".*/parquet/")
        i.river.files <- which(str_detect(made.files, i.name))
        i.river.files <- made.files[i.river.files]
        i.out <- vector(mode = "list", length = length(i.river.files))
        i.out <- lapply(i.river.files, st_read_parquet)
        i.out <- do.call(rbind, i.out)
        i.out.name <- paste0(i.name, "_hydrography90m_segments.parquet")
        st_write_parquet(i.out, paste0("E:/Arbeit/Data/river_network/hydrography90m/network_in_catchments/",i.out.name))
        #- remove old files 
        #file_delete(i.river.files)
        #- clean loop
        rm(list = ls()[grepl(pattern = "^i\\.", x = ls())])
        gc()
}
