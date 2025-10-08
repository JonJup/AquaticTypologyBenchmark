### FIT HMSC Models 

# set working directory ---------------------------------------------------
setwd(rstudioapi::getActiveProject())

# load libraries ----------------------------------------------------------
library(Hmsc)
library(coda)
library(data.table)

models <- list.files("data/unfitted_hmsc_models/", full.names = T)
source("code/functions/gelman_check.R")

## SET HMSC MCMC PARAMETERS
nChains = 3
samples = 2000
thin = 2
transient = round(samples*thin)


# 21 test error 
#for (i in seq_along(models)){
for (i in c(1,16)){
        
      print(i)        
      i.mod <- readRDS(models[i])  
      while.condition <- TRUE
      i.counter = 1
      while (while.condition) {
              print(paste("fitting model", i, "; try:", i.counter))
              # ———— fit model —————
             # Some models returned Failed updaters and their counts in chain x  ( 8000  attempts): GammaEta
             # Gamma is the influence of traits on niches (i.e. beta coefficients). 
             # Setting the updater to false is recommended in https://www.helsinki.fi/assets/drupal/2022-08/hmscdevelopment.pdf
             fit = sampleMcmc(
                      i.mod,
                      thin      = thin * i.counter,
                      samples   = samples * i.counter,
                      transient = transient * i.counter,
                      nChains   = nChains,
                      nParallel = nChains,
                      updater   = list(GammaEta = FALSE)
              )
              
             # ————— evaluate model fit ————
              #- check model fit with the potential scale reduction factor
              o.mpost <- convertToCodaObject(fit)
              # thinned to speed up 
              # Thin every 20th sample
              o.thinned_mcmc <- window(o.mpost, thin = 20)
              o.means        <- gelman_check(o.thinned_mcmc)
              rate <- sum(o.means$taxon$psrf>= 1.1) / nrow(o.means$taxon)
              
              if (rate >= 0.1) {
                      i.counter <- i.counter + 1
                      # identify the worst taxon 
                      setorderv(o.means$taxon, "taxon_mean")
                      worst.taxa <- tail(o.means$taxon, 3)
                      worst.taxa <- worst.taxa[taxon_mean >= 1.1, taxon]
                      
                      # # identify the worst variable 
                      setorderv(o.means$variable, "var_mean")
                      
                      worst.vars <- o.means$variable[var_mean >= 1.1]
                      worst.vars <- worst.vars[variable != "(Intercept)"]
                      worst.vars <- tail(worst.vars$variable, 3)
                      
                      if (length(worst.vars) > 0) {
                              newX <- i.mod$XData[, -which(colnames(i.mod$XData) %in% worst.vars)] 
                              newFormula <- paste("-", worst.vars, collapse = " ")
                              newFormula <- as.formula(paste("~ .", newFormula))
                              newFormula <- update(i.mod$XFormula, newFormula)
                      } 
                      if (length(worst.taxa) > 0) {
                              newY <- i.mod$Y[, -which(colnames(i.mod$Y) %in% worst.taxa)]
                      }
                      studyDesign <- data.frame(sample = as.factor(1:nrow(i.mod$XData)))
                      rL          <- HmscRandomLevel(units = o.studyDesign$sample)
                      
                      i.mod <- Hmsc(
                              Y           = newY,
                              XData       = newX,
                              XFormula    = newFormula,
                              studyDesign = studyDesign,
                              ranLevels   = list("sample" = rL),
                              distr       = ifelse(any(grepl("poisson",i.mod$call)), "lognormal poisson", "probit")
                              
                      )   
                      
                      
                      # i.mod$XScalePar <- i.mod$XScalePar[, -which(colnames(i.mod$XData) %in% worst.vars$variable)]
                      # i.mod$XData <- i.mod$XData[, -which(colnames(i.mod$XData) %in% worst.vars$variable)] 
                      # 
                      # i.mod$X <- i.mod$X[, -which(colnames(i.mod$X) %in% worst.vars$variable)] 
                      # i.mod$XScaled <- i.mod$XScaled[, -which(colnames(i.mod$XScaled) %in% worst.vars$variable)] 
                      # i.mod$covNames <- i.mod$covNames[- which(i.mod$covNames %in% worst.vars$variable)]
                      # i.mod$nc <- length(i.mod$covNames)
                      # i.mod$mGamma <- i.mod$mGamma[1:i.mod$nc]
                      # i.mod$UGamma <- i.mod$UGamma[1:i.mod$nc, 1:i.mod$nc] 
                      # i.mod$V0 <- i.mod$V0[1:i.mod$nc, 1:i.mod$nc]  
                      # 
                      # 
                      # i.mod$Y     <- i.mod$Y[, -which(colnames(i.mod$Y) %in% worst.taxa$taxon)]
                      # i.mod$YScaled     <- i.mod$YScaled[, -which(colnames(i.mod$YScaled) %in% worst.taxa$taxon)]
                      # i.mod$YScalePar <- i.mod$YScalePar[,1:ncol(i.mod$YScaled)]
                      # i.mod$spNames <- i.mod$spNames[- which(i.mod$spNames %in% worst.taxa$taxon)]
                      # i.mod$ny <- i.mod$ns <- length(i.mod$spNames)
                      # i.mod$aSigma <- i.mod$aSigma[1:i.mod$ny]
                      # i.mod$bSigma <- i.mod$bnSigma[1:i.mod$ny]
                      # 
                      # 
                      # removal_string <- paste("-", worst.vars$variable, collapse = " ")
                      # update_formula <- as.formula(paste(". ~ .", removal_string))
                      # i.mod$XFormula <- update(i.mod$XFormula, update_formula)
                      
              } else {
                      # print("accepted")
                      while.condition <- FALSE
                      out <- list(model = fit, converged = "yes")
              }
      
              #if (i.counter == 4) {
              #        while.condition == F
              #        out <- list(model = fit, converged = "no")
              #}
      }
              o.save.name <- sub(x = models[i], pattern = "unfitted", replacement = "fitted")
              saveRDS(out, o.save.name)
                
}
