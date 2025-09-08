### FIT HMSC Models 


# set working directory ---------------------------------------------------
setwd(rstudioapi::getActiveProject())

# load libraries ----------------------------------------------------------
library(Hmsc)
library(coda)

models <- list.files("data/unfitted_hmsc_models/", full.names = T)
source("code/functions/gelman_check.R")

## SET HMSC MCMC PARAMETERS
nChains = 3
samples = 2000
thin = 2
transient = round(samples*thin)


# 21 test error 
for (i in seq_along(models)){
        if (i == 1) next()
      print(i)        
      o.mod <- readRDS(models[i])  
      while.condition <- TRUE
      i.counter = 1
      while (while.condition) {
              print(paste("fitting model", i, i.counter))
              # ———— fit model —————
             # Some models returned Failed updaters and their counts in chain x  ( 8000  attempts): GammaEta
             # Gamma is the influence of traits on niches (i.e. beta coefficients). 
             # Setting the updater to false is recommended in https://www.helsinki.fi/assets/drupal/2022-08/hmscdevelopment.pdf
               o.m1.mcmc = sampleMcmc(
                      o.mod,
                      thin      = thin * i.counter,
                      samples   = samples * i.counter,
                      transient = transient * i.counter,
                      nChains   = nChains,
                      nParallel = nChains,
                      updater = list(GammaEta = FALSE)
              )
              # ————— evaluate model fit ————
              #- check model fit with the potential scale reduction factor
              o.mpost <- convertToCodaObject(o.m1.mcmc)
              # thinned to speed up 
              # Thin every 10th sample
              o.thinned_mcmc <- window(o.mpost, thin = 10)
              o.means        <- gelman_check(o.thinned_mcmc)
              o.means.rate <- sum(o.means >= 1.1) / length(o.means) 
              if (o.means.rate >= 0.1) {
                      while.condition <- TRUE
                      i.counter <- i.counter + 1
              } else {
                      # print("accepted")
                      while.condition <- FALSE
              }
              if (i.counter == 6) {
                      while.condition == F
                      o.model.fit = FALSE
              }
              o.save.name <- sub(x = models[i], pattern = "unfitted", replacement = "fitted")
              saveRDS(o.m1.mcmc, o.save.name)
              rm(list = ls()[grepl(pattern = "^o\\.", x = ls())])
      }           
}
