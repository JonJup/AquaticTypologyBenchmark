# n_types = c(3,5,10)
# max.q <- 3
# within.q <- 30
# neval <- max.q * within.q
# n_type_scenarios <- max.q * length(n_types)
# n_total <- length(n_types) * neval 


# how many types do the typologies have?
n_types = c(3,5,10)
# how many random selections of variables 
max.q <- 5
# how many simulated environments in each setting? These differ in how strong the classifications are.  
within.q <- 18
# so how many evaluations 
neval <- max.q * within.q
# what are the scenarios? 
n_type_scenarios <- max.q * length(n_types)
# another ????
n_total <- length(n_types) * neval 