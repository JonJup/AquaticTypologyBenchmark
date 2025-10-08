library(sf)
library(arrow)
library(data.table)
library(isotree)

setwd(rstudioapi::getActiveProject())
# load data ---------------------------------------------------------------
env <- list.files("C:/Users/jonjupke/OneDrive - University of Helsinki/PULSE/WP1/AquaticTypologyBenchmark/data/catchments/", full.names = T)
enz <- readRDS("E://Arbeit/data/typology_systems/environmental_zones/processedData/eu_hydro_dem_w_enz.rds")
schemes <- lapply(list.files("data/biota/", pattern = "03", full.names = T), readRDS)
schemes <- rbindlist(schemes, fill = T)
unique.enz <- unique(schemes$catchments)

out <- list()
for (i in seq_along(env)){
        print(i)
        i.cat <- read_parquet(env[i])
        i.join <- enz[i.cat, on = "ID"]
        i.join$EnZ_name <- factor(i.join$EnZ_name)
        i.join <- i.join[!is.na(EnZ_name)]
        i.join <- na.omit(i.join)
        out[[i]] <- i.join
}
out <- rbindlist(out, use.names = T)
out.sp <- split(out, by = "EnZ_name")
out2 <- list()
isotree.list <- list()
isotree.predictions <- list()
cov.list <- list()
envmean.list <- list()
mahathres.list <- list()
for (i in seq_along(unique.enz)){
        print(i)
        print(Sys.time())
        index <- unlist(unique.enz[i])
        i.out.sp <- which(names(out.sp) %in% index)
        i.out.sp <- out.sp[i.out.sp]
        i.out.sp <- rbindlist(i.out.sp)
        i.out.sp[, c("ID", "EnZ_name") := NULL]
        i.out.sp <- as.data.frame(i.out.sp)
        i.out.sp <- as.matrix(i.out.sp)
        i.out.sp <- i.out.sp[, - which(colnames(i.out.sp) %in% c("spi_max", "VBM_max"))]
        
        out2[[i]] <- i.out.sp
        names(out2)[i] <- toString(index)
        isotree.list[[i]] <- isolation.forest(
                out2[[i]],
                ntrees = 100,
                sample_size = 256,  # Subsample for efficiency
                ndim = ncol(out2[[i]]),
                prob_pick_pooled_gain = 0  # Use all dimensions
        )
        isotree.predictions[[i]] <- predict(isotree.list[[i]], out2[[i]])
        # Compute covariance for mahalanobis distance 
        cov.list[[i]] <- cov(i.out.sp)
        envmean.list[[i]] <- colMeans(i.out.sp)  
        mahathres.list[[i]] <- qchisq(0.99, df = ncol(i.out.sp))
}


names(isotree.list) <- 
names(isotree.predictions) <-
names(cov.list) <- 
names(envmean.list) <-
names(mahathres.list) <- 
names(out2)


iso_maha <- list(isotree  = isotree.list,
                 predictions = isotree.predictions,
                 covariance = cov.list,
                 means = envmean.list,
                 thresholds = mahathres.list)

# save to file 
saveRDS(iso_maha, "data/isomaha.rds")
# saveRDS(isotree.list, "data/isolation_forests_enz.rds")
# saveRDS(isotree.predictions, "data/isolation_forests_enz_predictions.rds")


