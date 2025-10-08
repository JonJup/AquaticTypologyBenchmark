################################################################################
# Script Name:        build_hmsc_hpc_models.R
# Description:        Variation of the build_hmsc_models.R script that builds models specifically for the HMSC-HPC version.
#                     See github: https://github.com/hmsc-r/hmsc-hpc and paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011914                      
#
#                     This script always represents the start of a model fitting iteration. 
#                     This iteration goes through multiple stages:
#                       1. creating the unfitted models and initialized models. 
#                       2. fitting the initialized models with Hmsc-HPC on server resulting in fitted models 
#                       3. evaluating fitted models. 
#                               |--- If they pass, they are saved in converged_models. 
#                               |--- If they fail, they are their entry in mcmc_parameters.rds is updated 
#                       4. Start again here. Refit all models that failed in 3. 

# Author:             Jonathan Jupke
# Date Created:       2025-09-16
# Last Modified:      2025-09-16
#
# R Version:          R 4.5.1
# Required Packages:  data.table, lubridate, sf, Hmsc, adespatial
#
# Notes:              This script is run locally. It is not computationally intensive. 
################################################################################

## DEFINED VARAIBLES ----- 

## LOOPS ----
# The script contains nested loops. These use to following indices to loop over the following objects:
# b - taxonomic groups 
# o - sampling schemes 

# setwd -------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())

# load packages -----------------------------------------------------------
library(data.table)
library(lubridate)
library(sf)
library(Hmsc)
library(adespatial)
library(jsonify)

# load custom functions ---------------------------------------------------
source("code/functions/determine_spatial_scale.R")


# data input ---------------------------------------------------
files        <- list.files("data/biota", pattern = "02_", full.names = TRUE)
schemes      <- list.files("data/biota", pattern = "03_", full.names = TRUE)
bio.list     <- lapply(files, readRDS)
schemes.list <- lapply(schemes, readRDS)
mcmcPar      <- readRDS("data/mcmc_parameters.rds")
if ("changes_for_refit.rds" %in% list.files("data")) varTodrop <- readRDS("data/changes_for_refit.rds")

# empty all folders that are part of the loop described above. 
# unlink("data/fitted_hmsc_models/*")
# unlink("data/unfitted_hmsc_models/*")
# unlink("data/initialized_hmsc_models/*")


# prepare fitting  --------------------------------------------------------
bio.names    <- 
        sapply(
                files, 
                function (x) sub(x = x, pattern = "data/biota/02_", replacement = "")
        ) |> 
        sapply(
                function(x) sub(x = x, pattern = "_w_environment.rds", replacement = "")
        )


# Prepare models  -------------------------------------------------------------------
# The first loop over b, loops over the taxonomic groups diatoms, fish, invertebrates, and macrophytes. 
# In that order.
# TODO describe what happens in it
for (b in 1:length(bio.names)){
        if (b < 4) next()
        # select list of sampling schemes for focal taxonomic group from schemes.list 
        b.scheme <- schemes.list[[b]]
        # select sample data for focal taxonomic group
        b.bio    <- bio.list[[b]]
        
        # The + in the names of some diatom complexes causes errors later on.
        # The + is removed from the names. 
        if (b == 1) b.bio[, working.taxon := gsub(x = working.taxon, pattern = "\\+", replacement = "")]
        
        # This loop is nested in b and loops over o, i.e., sampling schemes. 
        # These are the rows of b.scheme
        for (o in 1:nrow(b.scheme)) {
                
                
                # select sampling scheme 
                o.scheme <- b.scheme[o, ]
                # IF this model was already successfully fit skip it 
                o.par    <- mcmcPar[scheme_id == o.scheme$scheme_id] 
                if (o.par$converged) next()
                
                if(o > 5) next
                # TODO DEBUG CODE DELETE FOR GITHUB 
                #if (b == 1 & o == 1) next()
                
                # create number for filename 
                o.scheme.number <- ifelse(o > 9, as.character(o), paste0(0,o))
                o.scheme.number <- ifelse(o > 99, as.character(o.scheme.number), paste0(0,o.scheme.number))
                o.scheme.number <- ifelse(o > 999, as.character(o.scheme.number), paste0(0,o.scheme.number))
                
                # Set random number generation seed to ensure reproducability of random subsamples
                set.seed(o)
                
                # print update to console 
                print(paste("scheme", o, "of", nrow(b.scheme), "-", bio.names[b]))
                
                # select sampling scheme 
                o.scheme <- b.scheme[o, ]
                o.data   <- b.bio[data.set == o.scheme$data.set & 
                                        eventYear == o.scheme$eventYear & 
                                        month(eventDate) %in% as.numeric(o.scheme$focal_months[[1]])]
                # select MCMC parameters 
                o.para <- mcmcPar[scheme_id == o.scheme$scheme_id]
                
                # sample data according to scheme
                if (o.scheme$sample_type != "full"){
                        
                        o.sample.ids <- sample(
                                x       = unique(o.data$eventID),
                                size    = o.scheme$samples,
                                replace = F
                        )
                        o.data <- o.data[eventID %in% o.sample.ids]
                } 
                
                #- Remove rare taxa (less that 5 sites in total)
                #- Such a low number of occurrences would render the model inaccurate
                o.table <- sort(table(o.data$working.taxon))
                o.data  <- o.data[working.taxon %in% names(which(o.table >= 0.05 * uniqueN(o.data$eventID)))]
                
                # additionally remove taxa so that there are not more taxa than sampling sites 
                if (uniqueN(o.data$working.taxon) > o.scheme$samples){
                        o.table <- sort(table(o.data$working.taxon))
                        delta <- uniqueN(o.data$working.taxon) - o.scheme$samples
                        o.table <- names(o.table[1:delta])
                        o.data <- o.data[!working.taxon %in% o.table]
                }
                
                
                #- establish spatial scale
                o.sf <- determine_spatial_scale(o.data)
                
                #- save spatial scale to file
                saveRDS(object = o.sf,
                        file = paste0("data/spatial_scale/", bio.names[b], "_", o.scheme.number,".rds"))
                #- establish taxonomic resolution 
                o.tx <- data.table(
                        scheme_id       = paste0(bio.names[b], "_", o.scheme.number),
                        taxon           = o.scheme$taxon,
                        data.set        = o.scheme$data.set,
                        year            = o.scheme$eventYear,
                        samples         = o.scheme$samples,
                        organismQuantityType = o.scheme$organismQuantityType,
                        species_rank    = o.data[!is.na(species), .N]/nrow(o.data),
                        genus_rank      = o.data[is.na(species) & !is.na(genus), .N]/nrow(o.data),
                        family_rank     = o.data[is.na(species) & is.na(genus) & !is.na(family), .N]/nrow(o.data),
                        higher_rank     = o.data[is.na(species) & is.na(genus) & is.na(family), .N]/nrow(o.data)
                )
                saveRDS(object = o.tx,
                        file = paste0("data/taxonomic_resolution/", bio.names[b], "_", o.scheme.number,".rds"))
                #- extract types 
                o.types <- copy(o.data)
                o.types <- unique(o.types, by = "eventID")
                o.types <- o.types[, c("FEOW", "GLORiC", "BGR", "BRT", "EnZ", "HER", "IFE")]
                saveRDS(object = o.types,
                        file = paste0("data/scheme_types/", bio.names[b], "_", o.scheme.number,".rds"))
                # prepare data for HMSC ----------------------------------------------------------------------
                # HMSC needs a Y response matrix.
                # That matrix has one sample per row and one column per taxon
                # We need to reshape the data for this.
                # At this step we can check whether the
                # if (unique(o.data$abundance_type) == "relative abundance"){
                #         ?
                # }
                if (o.scheme$organismQuantityType == "presence"){
                        o.data2 <- dcast(
                                o.data,
                                formula = eventID ~ working.taxon,
                                value.var = "PA",
                                fill = 0
                        )
                } else if (o.scheme$organismQuantityType == "individuals"){
                        o.data2 <- dcast(
                                o.data,
                                formula = eventID ~ working.taxon,
                                value.var = "organismQuantity",
                                fill = 0
                        )
                }
                # remove site id to get a purely numeric matrix 
                o.data3 <- o.data2[,-1]
                o.data3 <- as.matrix(o.data3)
                #- prepare environmental variables
                o.env.samples <- unique(o.data, by = "eventID")

                o.env.samples <- o.env.samples[, c(
                        "x.coord", "y.coord","slope", "lake_area",
                        "soil_oc", "soil_pH", "inundation_depth", "flooded_area",
                        "spi_mean", "area_calcareous", "area_siliceous",
                        "area_sediment", "glacial_area", "bioclim05", "bioclim06",
                        "bioclim13", "bioclim14", "Rfactor_max", "Rfactor_avg",
                        "Rfactor_min", "roughness", "elevation", "mean_snow_equivalent",
                        "mean_discharge", "max_discharge", "min_discharge", "groundwater_table",
                        "upstream_catchment_area", "saturated_soil_water_content", "VBM_mean"
                )]
                #- If glacial area is zero everywhere, this leads to issues later on.
                #- Therefore, this variable is removed in such cases.
                o.cols.to.keep <- o.env.samples[, names(which(lapply(.SD, uniqueN) > 1))]
                o.env.samples  <- o.env.samples[, ..o.cols.to.keep]
                rm(o.cols.to.keep)
                
                # Saturated soil water content has some missing values 
                # These are imputed here with a OLS Regression
                if (any(is.na(o.env.samples$saturated_soil_water_content))) {
                        o.row <- which(is.na(o.env.samples$saturated_soil_water_content))
                        o.pre <- o.env.samples[-o.row, ]
                        o.mod <- lm(saturated_soil_water_content ~ ., data = o.pre)
                        o.pre <- o.pre[o.row]
                        o.pre <- predict(object = o.mod, newdata =  o.pre)
                        o.env.samples[is.na(saturated_soil_water_content), saturated_soil_water_content := o.pre]
                        rm(o.row, o.pre,o.mod)
                }
                #- scale predictor variables 
                o.env.samples.xy <- o.env.samples[, c("x.coord", "y.coord")]
                o.env.samples[, c("x.coord", "y.coord") := NULL]
                saveRDS(o.env.samples, paste0("data/unscaled_environments/", bio.names[b], "_", o.scheme.number,".rds"))
                
                o.env.samples <- scale(o.env.samples)
                attributes(o.env.samples)$`scaled:scale` <- NULL
                attributes(o.env.samples)$`scaled:center` <- NULL
                o.env.samples <- cbind(o.env.samples.xy, o.env.samples)
                
                # Morans Eigenvector Maps 
                o.xy.mat <- matrix(
                        c(o.env.samples$x.coord, o.env.samples$y.coord),
                        ncol = 2,
                        byrow = F
                )
                colnames(o.xy.mat) <- c("x", "y")
                o.x       <- dbmem(
                        xyORdist    = o.xy.mat,
                        MEM.autocor = "non-null",
                        store.listw = T
                )
                #- Following Bini et al (2009) [10.1111/j.1600-0587.2009.05717.x]. In the paper this is called SEVM-3. 
                #- Loops over species 
                # Initialize a counter for MEM significance before your loop
                o.mem.significance.counts <- numeric(ncol(o.x))
                names(o.mem.significance.counts) <- colnames(o.x)
                for (focal.spe in 1:ncol(o.data3)){
                        
                        # If a species is always present we cannot model it here. 
                        # The corresponding iteration is skipped 
                        if (o.scheme$organismQuantityType == "presence" & all(o.data3[, focal.spe] == 1)) next()
                        
                        # remove coordinates
                        fs.env.samples <- copy(o.env.samples)
                        fs.env.samples[, c("x.coord", "y.coord") := NULL]
                        
                        if (sum(o.data3[, focal.spe] != 0)/nrow(o.data3) < 0.05) next()
                        fs.env.samples <- cbind(o.data3[, focal.spe], fs.env.samples)
                        fs.predictors  <- names(fs.env.samples)[-1]
                        fs.formulas    <- paste("V1 ~", paste(fs.predictors, collapse = "+"))
                        fs.formulas    <- as.formula(fs.formulas)
                        if (o.scheme$organismQuantityType == "presence"){
                                fs.env.samples$V1 <- factor(fs.env.samples$V1)
                                fs.model       <- glm(fs.formulas, data = fs.env.samples, family = "binomial")
                        } else if (o.scheme$organismQuantityType == "individuals"){
                                fs.model <- glm(fs.formulas, data = fs.env.samples, family = "poisson")
                        }
                        
                        fs.residuals   <- residuals(fs.model)
                        fs.p.values    <- apply(o.x, 2, function(x) cor.test(x,fs.residuals)$p.value)
                        fs.p.values    <- p.adjust(fs.p.values, method = "holm")
                        # Update the counter for significant MEMs
                        if (any(fs.p.values < 0.005)) {
                                fs.sp <- which(fs.p.values < 0.005)
                                o.mem.significance.counts[fs.sp] <- o.mem.significance.counts[fs.sp] + 1
                         }
                        rm(list = ls()[grepl(pattern = "^fs\\.", x = ls())])
                }
                rm(focal.spe)
                o.top5.mems <- sort(o.mem.significance.counts, decreasing = TRUE)[1:5]
                if (any(o.top5.mems == 0)) o.top5.mems <- o.top5.mems[-which(o.top5.mems == 0)]
                if (length(o.top5.mems) != 0){
                        o.top5.mem.names   <- names(o.top5.mems)
                        o.top5.mem.indices <- which(names(o.mem.significance.counts) %in% o.top5.mem.names) 
                        o.MEM            <- o.x[, o.top5.mem.indices] |> 
                                                 data.frame() |>
                                                setDT()
                        names(o.MEM) <- o.top5.mem.names
                        o.env.samples <- cbind(o.env.samples, o.MEM)
                        o.n.env       <- ncol(o.env.samples) - 2 - length(o.top5.mem.names)
                } else {
                  o.n.env <- ncol(o.env.samples)
                }
                
                # Prepare HMSC model 
                #- establish a site level random factor
                o.studyDesign <- data.frame(sample = as.factor(1:nrow(o.data3)))
                o.rL          <- HmscRandomLevel(units = o.studyDesign$sample)
                
                #- create model formula
                o.env.samples[, c("x.coord", "y.coord") := NULL]
                o.predictors  <- names(o.env.samples)
                o.formulas    <- paste("~", paste(o.predictors, collapse = "+"))
                o.formulas    <- as.formula(o.formulas)
                
                #- define HMSC model
                if (o.scheme$organismQuantityType == "presence") {
                        o.mod1 <- Hmsc(
                                Y           = o.data3,
                                XData       = o.env.samples,
                                XFormula    = o.formulas,
                                studyDesign = o.studyDesign,
                                ranLevels   = list("sample" = o.rL),
                                distr       = "probit"
                                
                        )
                } else if (o.scheme$organismQuantityType == "individuals"){
                        o.mod1 <- Hmsc(
                                Y           = o.data3,
                                XData       = o.env.samples,
                                XFormula    = o.formulas,
                                studyDesign = o.studyDesign,
                                ranLevels   = list("sample" = o.rL),
                                distr       = "lognormal poisson"
                                
                        )   
                }
                o.save.name <-paste0("data/unfitted_hmsc_models/", bio.names[b], "_", o.scheme.number,".rds")
                saveRDS(o.mod1, o.save.name)
                # This step is the main difference to using vanilla Hmsc.
                # We initialize the model without fitting it. 
                # This is forced by the engine = "HPC" argument.
                o.init_obj = sampleMcmc(
                        o.mod1,
                        samples = o.para$nSamples,
                        thin = o.para$thin,
                        transient = o.para$transient,
                        nChains = o.para$nChains,
                        verbose = 100,
                        engine = "HPC"
                )
               
                o.init_obj <- to_json(o.init_obj)
                o.save.name <-paste0("data/initialized_hmsc_models/", bio.names[b], "_", o.scheme.number,".rds")
                saveRDS(o.init_obj, o.save.name)
                rm(list = ls()[grepl("^o\\.", x = ls())])
        }
}
