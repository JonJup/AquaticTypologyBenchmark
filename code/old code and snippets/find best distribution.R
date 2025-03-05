#  — — — find best distribution
find_best_distribution <- function(x){
        #- selection of all evaluated distributions 
        distributions <- c("norm", "weibull", "gamma", 
                           "lnorm", "beta", "pois", 
                           "exp", "logis", "unif", 
                           "chisq", "t", "f", 
                           "nbinom", "geom", "hyper", 
                           "binomial")
        # Perform fits and store results
        fits <- lapply(distributions, function(dist) {
                tryCatch(
                        fitdist(x, dist),
                        error = function(e) {
                                return(NULL)
                        }
                )
        })
        results <- data.frame(
                Distribution = distributions,
                AIC = sapply(fits, function(x) ifelse(!is.null(x), x$aic, NA)),
                para = I(sapply(fits, function(x) coef(x)))
        )
        return(results[which.min(results$AIC), ])
}
