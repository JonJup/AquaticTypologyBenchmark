gelman_check <- function(posterior){
        o.ge    <- gelman.diag(x = posterior$Beta, multivariate = FALSE)
        #- species names
        temp <- rownames(o.ge$psrf)
        temp <- regmatches(temp, regexpr(",\\s.*\\(S", temp))
        temp <- gsub("^,\\s*", "", temp)
        temp <- gsub("\\(S$", "", temp)
        o.species_names <- unique(temp)
        #- Create empty vectors to store results
        o.means <- c()
        #- Calculate statistics for each species
        for (j in seq_along(o.species_names)) {
                # Get rows corresponding to current species
                j.species_rows <- grep(o.species_names[j], rownames(o.ge$psrf))
                # Calculate mean and max of point estimates for these rows
                o.means <- append(o.means, c(o.ge$psrf[j.species_rows, "Point est."]))
                rm(list = ls()[grepl("^j\\.", x = ls())])
        }
        rm(j)
        return(o.means)
}