# setup -----------------------------------------------------------------------------
library(groundhog)
pkgs <- c(
        "adespatial", 
        "arrow", 
        "cluster",
        "data.table",
        "DirichletReg",
        "fitdistrplus",
        "fs",
        "Hmsc", 
        "indicspecies",
        "labdsv",
        "magrittr",
        "parallelDist",
        "sirt",
        "sf",
        "sfarrow",
        "tidyverse",
        "units",
        "vegan",
        "zetadiv"
)
groundhog.library(pkgs,'2024-12-01')
rm(pkgs)

# DEFINE CUSTOM FUNCTIONS --------------------------------------------------------------

#' TABLE OF CONTENTS: 
#' 1. find_best_distribution()
#' 2. 

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

#  — — — simulate environmental variables
#' @param x model output from 04_HMSC.R
#' @param y which variables to use? Names of columns
#' @param types integer. How many types should the artificial typology system have?
#' @param obs_per_type integer. How many observations are included in each type?
simulate_env <- function(x,y,types,strength,obs_per_type = j.tv){

        if (length(obs_per_type) == 1){
                n_obs_per_type <- obs_per_type
        } else {
                n_obs_per_type <- table(obs_per_type)
        }
        
        breaks <- 1:types/(types+1)
        y2     <- x[[3]]$X[, which(colnames(x[[3]]$X) == y)]
        
        #- find most appropriate distribution. 
        #- This uses the find_best_distribution_function defined above. 
        dist <- find_best_distribution(y2)
        if (dist$Distribution %in% c("norm", "weibull")){
                
                means   <- sapply(breaks, 
                                  FUN = function(k) quantile(y2, k)
                )
                denominator <-
                        case_when(
                                strength == 1 ~ 1, 
                                strength == 2 ~ j.types, 
                                strength == 3 ~ j.types * 10
                        )
                new_val <- map2(
                        means,
                        n_obs_per_type,
                        ~rnorm(
                                n = .y,
                                mean = .x,
                                sd = sd(y2)/denominator
                        )
                )
                
        } else if (dist$Distribution == "exp") {
                
                base_values <- rexp(1000, rate = dist$para[[1]])
                multiplicator <- case_when(
                        strength == 1 ~ 2,
                        strength == 2 ~ 1.2,
                        strength == 3 ~ 1,
                        )
                
                means   <- sapply(
                        breaks,
                        FUN = function(k)
                                quantile(base_values, k)
                )
                new_val <- map2(means,
                                n_obs_per_type,
                                ~ sample(
                                        x = base_values[which(abs(base_values - .x) < abs(.x - sd(base_values) * multiplicator))],
                                        size = .y,
                                        replace = T
                                        
                                ))
        } else if (dist$Distribution == "lnorm"){
                means <- sapply(breaks, 
                                  FUN = function(k) quantile(y2, k)
                )
                denominator <-
                        case_when(
                                strength == 1 ~ 1, 
                                strength == 2 ~ j.types, 
                                strength == 3 ~ j.types * 10
                        )
                new_val <- map2(
                        means,
                        n_obs_per_type,
                        ~rlnorm(
                                n = .y,
                                mean = log(.x),
                                sd = dist$para[[1]][2]/denominator
                        )
                )
                
                new_val <- unlist(new_val)
                new_val <- unname(new_val)
                check <- data.frame(value = c(new_val, y2), type = rep(c("sim", "meas"), each = 100))
                check$class <- c(j.tv, rep(0, 100))
                check$rank <- rank(check$value)
                # ggplot(check, aes(x = rank, y = value, col = type)) + 
                #         geom_jitter(height = 1000)
                ggplot(check, aes(x = rank, y = value, col = factor(class))) + 
                        geom_jitter(height = 1000)
                
        } else {
                stop(paste(dist$Distribution, "distribution not yet implemented"))
        }

        new_val <- unlist(new_val)
        new_val <- unname(new_val)
        return(new_val)
}
return(simulate_env)

calculate_auc <- function(x, y) {
        # Sort x and y together (in case they're not already ordered)
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        
        # Calculate width of each trapezoid
        dx <- diff(x)
        
        # Calculate mean height of each trapezoid
        mean_height <- (y[-1] + y[-length(y)]) / 2
        
        # Calculate area as sum of trapezoid areas
        area <- sum(dx * mean_height)
        
        return(area)
}
return(calculate_auc)
my.cs <- function(x){
        out <- meandist(dist = x, grouping = j.tv)
        out <- unlist(summary(out))["CS"]
        return(out)
}
return(my.cs)
#- wrapper to compute average silhouette width
my.asw <- function(x){
        out <- silhouette(dist = x, x = j.tv)
        out <- mean(out[, "sil_width"])
        return(out)
}
return(my.asw)
#- wrapper for metrics based on indicator species analysis
my.isa <- function(x){
        out <- multipatt(x, 
                         j.tv, 
                         func = "indval",
                         control = how(nperm = 500))
        out1 <- nrow(out$sign[which(out$sign$p.value<=0.05), ])/nrow(out$sign)
        out2 <- mean(out$sign$p.value, na.rm = T)
        return(c(out1, out2))
}
return(my.isa)
#- wrapper for ISAMIC
my.isamic <- function(x){
        out <- isamic(
                comm = x, 
                clustering = j.tv)
        out <- mean(out)
        return(out)
}
return(my.isamic)
#- wrapper for ANOSIM
my.anosim <- function(x){
        out <- anosim(x, grouping = j.tv, permutations = 500, parallel = 6)
        out1 <- out$statistic
        out2 <- out$signif
        return(c(out1, out2))
}
return(my.anosim)
