param_combinations <- data.frame(
  centroid_adjustment = runif(n = within.q, min = -1, max = 1), 
  compactness         = runif(n = within.q, min = -1, max = 1)
)
q.newenv <-
  pmap(
    param_combinations,
    function(centroid_adjustment, compactness) {
      adjust_clusters(
        centroids                  = q.centroids,
        observations               = q.observations,
        cluster_assignments        = q.clusters,
        centroid_adjustment_factor = centroid_adjustment,
        compactness_factor         = compactness
      )
    }
  )


# silhouette widths 
asw <- sapply(q.newenv, function(x) {
        y = stats::dist(x)
        y = silhouette(q.clusters, dist = y)
        y = y[,3]
        y = mean(y)
        return(y)
})
i.asw %<>% append(asw)
# add other variables back to the mix 
q.missing2 <- q.env[, .SD, .SDcols = !q.all.vars]
q.missing2[, `:=` (type = NULL)]
q.missing2 %<>% data.frame %>% as.matrix
q.newenv %<>% lapply(function(x) cbind(x, q.missing2))
s.order <- names(i.model$environment)
q.newenv %<>% lapply(function(x) {
        x <- data.frame(x)
        x <- setDT(x)
        x <- x[, ..s.order]
        x %<>% data.frame
        x %<>% as.matrix
        x <- cbind(rep(1, nrow(x)), x)
        colnames(x)[1] <- "(Intercept)"
        return(x)
})
i.contraction.points %<>% append(param_combinations$compactness)
i.contraction.centroids %<>% append(param_combinations$centroid_adjustment)
