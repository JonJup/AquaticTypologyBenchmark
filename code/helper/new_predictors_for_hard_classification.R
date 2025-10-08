source("code/functions/adjust_clusters.R")

# What is the original level of 'normal'?
# load unscaled environment 
q.env <- list.files("data/unscaled_environments/", full.names = T)
q.ue <- grep(x = q.env, pattern = model.name)
q.ue  <- readRDS(q.env[q.ue])


# TODO remove me 
q.ue[, c("VBM_max", "spi_max") := NULL]
if ("comb_date_id" %in% names(q.ue)){
        q.ue[, comb_date_id := NULL]
}
if (
        !"glacial_area" %in% names(q.ue) | 
        !"mean_snow_equivalent" %in% names(q.ue) 
    )
        {
        
        if (!"glacial_area" %in% names(q.ue)) q.ue[, glacial_area := 0]
        if (!"mean_snow_equivalent" %in% names(q.ue)) q.ue[, mean_snow_equivalent := 0]
        setcolorder(q.ue, names(isomaha$means[[1]]))

}
# load environmental zones 
q.enz <- unique(unlist(i.schemes$catchments))
q.enz <- toString(q.enz)
# evaluate original environment with isolation forest 
q.idx <- which(names(isomaha$isotree) == q.enz)
q.iso <- isomaha$isotree[[q.idx]]
q.pre <- isomaha$predictions[[q.idx]]
q.scores <- predict(q.iso, q.ue)
q.scores.mean  <- mean(q.scores)
q.scores.sd  <- sd(q.scores)

# TODO JUST FOR TESTS! 
within.q.run = 300
within.q = 30

param_combinations <- data.frame(
  separation   = runif(n = within.q.run, min = 0.0, max = 0.25), 
  compactness  = runif(n = within.q.run, min = -.5, max = .5)
)

# validate param combinations against real environment  -------------------
q.names.ue <- names(q.ue)
q.ue[, type := q.clusters]


#- find the centroid for each variable that follows the typology
q.ue.centroids <- q.ue[, lapply(.SD, mean), .SDcols = q.all.vars, by = "type"]
setorderv(q.ue.centroids, "type")
q.ue.centroids[, type := NULL]
q.ue.centroids <- as.matrix(as.data.frame(q.ue.centroids))
q.ue.cluster.env <- as.data.frame(q.ue)
q.ue.cluster.env <- q.ue.cluster.env[, q.all.vars]
q.ue.cluster.env <- as.matrix(q.ue.cluster.env)
q.ue.missing2 <- q.ue[, .SD, .SDcols = !q.all.vars]

q.ue.missing2[, `:=` (type = NULL)]

q.newenv <- list()
for (adjclu in 1:within.q.run) {
        result <-
                adjust_clusters(
                        centroids                 = q.ue.centroids,
                        observations              = q.ue.cluster.env,
                        cluster_assignments       = q.clusters,
                        separation_factor  = param_combinations$separation[adjclu],
                        compactness_factor = param_combinations$compactness[adjclu]
                )
        cef <- calculate_effective_factors(
                original_centroids   = q.ue.centroids,
                constrained_centroids = result$centroids,
                original_observations = q.ue.cluster.env,
                constrained_observations = result$observations,
                cluster_assignments = q.clusters
        )
        param_combinations$separation[adjclu] <- cef$effective_centroid_factor
        param_combinations$compactness[adjclu] <- cef$effective_compactness_factor
        q.newenv[[adjclu]] <- result$observations
}
rm(cef)
rm(adjclu)
rm(result)

# Isolation forest --------------------------------------------------------
q.rating.out <- list()
for (jj in 1:within.q.run){
                
        # add back variables that are not affected by the dilations  
        jj.data            <- cbind(q.newenv[[jj]], q.ue.missing2)
        # remove Moran's Eigenvector Maps 
        jj.columns_to_keep <- setdiff(names(jj.data), 
                                              grep("MEM", 
                                                   names(jj.data), 
                                                   value = TRUE)
                                              )
        jj.data            <- jj.data[, ..jj.columns_to_keep]
        if ("type" %in% q.names.ue) q.names.ue <- q.names.ue[-which(q.names.ue == "type")]
        jj.data            <- jj.data[, ..q.names.ue]
        jj.scores             <- predict(q.iso, jj.data)
        jj.ratings <- data.table(
                scores = jj.scores, 
                scores_centered = jj.scores - q.scores.mean,
                separation = param_combinations$separation[jj],
                compactness = param_combinations$compactness[jj]
                   )
        
        jj.ratings[, scores_z := scores_centered/q.scores.sd]
        jj.ratings[, id := jj]
        q.rating.out[[jj]] <- jj.ratings
        rm(list = ls()[grepl(pattern = "^jj\\.", x = ls())])
        if (jj == within.q.run) q.rating.out <- rbindlist(q.rating.out)
}

# find those combinations that consistently stayed within 1 SD 
q.rating.out2 <- q.rating.out[scores_z <= 1]
q.rating.out2[, id_count := .N, by = "id"]
q.rating.out2[, avg_z := mean(scores_z), by = "id"]
q.rating.out2 <- q.rating.out2[id_count >= floor(i.schemes$samples - 0.25 * i.schemes$samples)]
q.rating.out2 <- unique(q.rating.out2, by = "id")
# print(paste(nrow(q.rating.out2), "passed isolation forest"))
setorderv(q.rating.out2, "avg_z")

# subset simulated environments to good ones 
q.newenv2 <- q.newenv[q.rating.out2$id]
if (length(q.newenv2) == 0) {
        next()
}

# Mahalanobis -------------------------------------------------------------

for (jj in 1:length(q.newenv2)){
        
        # add back variables that are not dilated 
        jj.data            <- cbind(q.newenv2[[jj]], q.ue.missing2)
        # remove Moran's Eigenvector Maps 
        jj.columns_to_keep <- setdiff(names(jj.data), 
                                      grep("MEM", 
                                           names(jj.data), 
                                           value = TRUE)
        )
        jj.data            <- jj.data[, ..jj.columns_to_keep]
        if ("type" %in% q.names.ue) q.names.ue <- q.names.ue[-which(q.names.ue == "type")]
        jj.data            <- jj.data[, ..q.names.ue]
        cor_mat <- cov2cor(isomaha$covariance[[q.idx]])
        if (any(cor_mat > 0.98)) {
                
                covmat <- isomaha$covariance[[q.idx]] + diag(1e-6, nrow(isomaha$covariance[[q.idx]]))
                m_dist <- mahalanobis(jj.data, center = isomaha$means[[q.idx]], cov = covmat)
        } else{
                m_dist <- mahalanobis(jj.data, center = isomaha$means[[q.idx]], cov = isomaha$covariance[[q.idx]])
        }
        rm(cor_mat, covmat)
        q.rating.out2$id[jj] <- jj
        # how many sites are above the predetermined threshold? 
        q.rating.out2$failsum[jj] <- sum(m_dist > isomaha$thresholds[[q.idx]])
        rm(m_dist)
        # jj.ratings <- data.table(
        #         scores = jj.scores, 
        #         scores_centered = jj.scores - q.scores.mean,
        #         separation = param_combinations$separation[jj],
        #         compactness = param_combinations$compactness[jj]
        # )
        # 
        # jj.ratings[, scores_z := scores_centered/q.scores.sd]
        # jj.ratings[, id := jj]
        # q.rating.out[[jj]] <- jj.ratings
        rm(list = ls()[grepl(pattern = "^jj\\.", x = ls())])
        # if (jj == within.q.run) q.rating.out <- rbindlist(q.rating.out)
}
rm(jj)
# In which samples are more than 75% of samples within 99% of the ChiSQ distribution.
q.rating.out2[, passed_mahalanobis := failsum < (0.1 * i.schemes$samples)]
q.rating.out2 <- q.rating.out2[passed_mahalanobis == 1]
# print(paste(nrow(q.rating.out2), "passed mahalanobis"))
if (nrow(q.rating.out2) == 0) next()
# subset to 30 if more are still available
if (nrow(q.rating.out2) > 30){
        q.rating.out2 <- q.rating.out2[1:30, ]
}
q.newenv3 <- q.newenv2[q.rating.out2$id]

if ("matrix" %in% class(q.newenv3)){
        q.newenv3 <- list(q.newenv3)
}

#- compute the silhouette widths to judge quality of the clusters
asw <- sapply(q.newenv3, function(x) {
        y = stats::dist(x)
        y = silhouette(q.clusters, dist = y)
        y = y[,3]
        y = mean(y)
        return(y)
})
i.asw <- append(i.asw, asw)
#- add other variables back to the mix
q.missing2 <- q.env[, .SD, .SDcols = !q.all.vars]
q.missing2[, `:=` (type = NULL)]
q.missing2 <- as.matrix(data.frame(q.missing2))
q.newenv <- lapply(q.newenv3, function(x) cbind(x, q.missing2))
s.order <- names(i.model$XData)
q.newenv <- lapply(q.newenv, function(x) {
        x <- data.frame(x)
        setDT(x)
        x <- x[, ..s.order]
        x <- as.matrix(data.frame(x))
        x <- cbind(rep(1, nrow(x)), x)
        colnames(x)[1] <- "(Intercept)"
        return(x)
})

i.contraction.points    <- append(i.contraction.points, q.rating.out2$compactness)
i.contraction.centroids <- append(i.contraction.centroids, q.rating.out2$separation)
rm(s.order)
rm(asw)
rm(param_combinations)
rm(q.missing2)
rm(q.newenv2)
rm(q.newenv3)