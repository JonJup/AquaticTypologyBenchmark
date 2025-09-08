# setup -------------------------------------------------------------------
source("code/functions/adjust_clusters.R")


param_combinations <- data.frame(
  separation   = runif(n = within.q, min = -1, max = 1), 
  compactness  = runif(n = within.q, min = -1, max = 1)
)
q.newenv <- mapply(
        function(separation, compactness) {
                adjust_clusters(
                        centroids                  = q.centroids,
                        observations               = q.observations,
                        cluster_assignments        = q.clusters,
                        centroid_adjustment_factor = centroid_adjustment,
                        compactness_factor         = compactness
                )
        },
        separation = param_combinations$separation,
        compactness = param_combinations$compactness,
        SIMPLIFY = FALSE
)

#- compute the silhouette widths to judge quality of the clusters 
asw <- sapply(q.newenv, function(x) {
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
q.newenv <- lapply(q.newenv, function(x) cbind(x, q.missing2))
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

i.contraction.points <- append(i.contraction.points, param_combinations$compactness)
i.contraction.centroids <- append(i.contraction.centroids, param_combinations$separation)
rm(s.order)
rm(asw)
rm(param_combinations)

