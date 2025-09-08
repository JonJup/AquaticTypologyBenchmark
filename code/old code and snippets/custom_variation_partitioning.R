#- custom variance partitioning 

#- works but I wont implement this. 
#- not worth the hassle. 

# Assuming hmsc_model is your fitted model
# First ensure you've run convertToCodaObject to prepare posteriors
hmsc_model_post <- convertToCodaObject(o.m1.mcmc)

# Set dimensions
n_samples <- niter(hmsc_model_post$Beta) # Number of posterior samples
n_sites <- nrow(o.m1.mcmc$Y)            # Number of sites
n_species <- ncol(o.m1.mcmc$Y)          # Number of species
n_covariates <- ncol(o.m1.mcmc$X)       # Number of covariates

# Initialize arrays for storing variance components
var_partition_posterior <- array(0, dim=c(n_samples, n_species, n_covariates))
var_linear_predictor_posterior <- array(0, dim=c(n_samples, n_species))

# Loop through each posterior sample
for(s in 1:n_samples) {
        # Extract Beta (regression coefficients) for this sample
        # Here we assume Beta is the first chain for simplicity
        beta_s <- hmsc_model_post$Beta[[1]][s,]
        
        # Reshape Beta to match dimensions [covariates, species]
        beta_matrix <- matrix(beta_s, nrow=n_covariates, ncol=n_species)
        
        # For each species
        for(sp in 1:n_species) {
                # Calculate linear predictor components
                linear_terms <- matrix(0, nrow=n_sites, ncol=n_covariates)
                
                for(j in 1:n_covariates) {
                        linear_terms[,j] <- o.m1.mcmc$X[,j] * beta_matrix[j,sp]
                }
        
                # Calculate total linear predictor for this species
                linear_predictor <- rowSums(linear_terms)
                
                # Store variance of the linear predictor
                var_linear_predictor_posterior[s,sp] <- var(linear_predictor)
                
                # Calculate and store variance of each component
                for(j in 1:n_covariates) {
                        var_partition_posterior[s,sp,j] <- var(linear_terms[,j])
                }
        }
}

# Normalize to get variance partitions
var_partition_normalized <- var_partition_posterior / 
        array(var_linear_predictor_posterior, dim=c(n_samples, n_species, n_covariates))
#- one row per species, one column per variable 
mean_var_partition <- apply(var_partition_normalized, c(2,3), mean)
rownames(mean_var_partition) <- colnames(o.m1.mcmc$Y)
colnames(mean_var_partition) <- colnames(o.m1.mcmc$X)
# 95% credible intervals
lower_CI <- apply(var_partition_normalized, c(2,3), quantile, 0.025)
upper_CI <- apply(var_partition_normalized, c(2,3), quantile, 0.975)
