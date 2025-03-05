#' This function was used to simulate new values from an environmental variable. 
#' I have replaced it with using the actually measured variables. 


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