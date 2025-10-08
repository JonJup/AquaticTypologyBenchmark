# Asymmetric Eigenvector Maps -------------------------------------------------------
                # Load river network for the EU HYDRO DEM catchments in which samples are located.
                # o.rivers <- load_river_networks(o.data)
                # o.sites  <- o.data |>
                #         unique(by = "sample_id") |>
                #         st_as_sf(coords = c("x.coord", "y.coord"), crs = 3035) |>
                #         st_transform(4326) 
                # o.network <- prepare_directional_network(o.rivers, o.sites)
                # # Calculate network costs
                # o.cost_matrix <-
                #         st_network_cost(x       = o.network,
                #                         from    = o.sites,
                #                         to      = o.sites,
                #                         weights = "weight") |>
                #         drop_units()
                # o.cost_matrix <- create_directional_incidence(o.cost_matrix, o.sites)
                # o.aem_result <- aem(binary.mat = o.cost_matrix$incidence,
                #                     weight     = o.cost_matrix$weights)
                # #- select AEM vectors 
                # o.signi.x <- c()
                # for (focal.spe in 1:ncol(o.data3)){
                #         
                #         fs.env.samples <- o.env.samples[, -c(1, 2)]
                #         fs.env.samples <- cbind(o.data3[, focal.spe], fs.env.samples)
                #         fs.env.samples$V1 %<>% factor
                #         fs.predictors  <- names(fs.env.samples)[-1]
                #         fs.formulas    <- paste("V1 ~", paste(fs.predictors, collapse = "+"))
                #         fs.formulas    <- as.formula(fs.formulas)
                #         fs.model       <- glm(fs.formulas, data = fs.env.samples, family = "binomial")
                #         fs.residuals   <- residuals(fs.model)
                #         fs.p.values    <- apply(o.aem_result$vectors, 2, function(x) cor.test(x,fs.residuals)$p.value)
                #         fs.p.values    <- p.adjust(fs.p.values, method = "holm")
                #         if (any(fs.p.values < 0.05)){
                #                 fs.sp <- which(fs.p.values < 0.05)
                #                 o.signi.x <- append(o.signi.x, fs.sp)
                #                 o.signi.x <- unique(o.signi.x)
                #         }
                #         rm(list = ls()[grepl(pattern = "^fs\\.", x = ls())])
                # }
                # o.AEM     <- 
                #         o.aem_result$vectors[, o.signi.x] |> 
                #         data.frame |> 
                #         setDT
                # names(o.AEM) <- paste0("AEM", 1:ncol(o.AEM))
                # o.env.samples <- cbind(o.env.samples, o.AEM)