# Evaluate simulated communities 

# Evaluate hard classification ------------------------------------------------------
                #- compute distance matrix for each simulated community
                j.distance.matrix  <- 
                        rapply(
                                j.storage, 
                                parallelDist, 
                                method = "binary",
                                how = "list"
                        )
                
                j.distance.matrix <- unlist(j.distance.matrix, recursive = F)
                j.distance.matrix <- unlist(j.distance.matrix, recursive = F)
                j.storage2 <- unlist(j.storage, recursive = F)
                j.storage2 <- unlist(j.storage2, recursive = F)
                #- compute Classification strength
                for (cs in 1:neval){
                        if (cs == 1) j.cs <- vector(mode = "numeric", length = neval)
                        out <- meandist(dist = j.distance.matrix[[cs]], grouping = j.cluster.assignments[[ceiling(cs/(neval/max.q))]])
                        j.cs[cs] <- unlist(summary(out))["CS"]
                }
                j.cs %<>% render_table(variable = "cs")
                #  ——— compute pairwise ANOSIM
                #- how many combinations are there
                j.perms <- 
                        permutations(
                                v = unique(
                                        j.cluster.assignments[[1]]
                                ),
                                k = 2
                                ) %>% 
                        data.frame %>%
                        filter(X1 < X2)
                
                
                for (evaluation in 1:neval) {
                        print(evaluation)
                        if (evaluation == 1){
                                j.out.p.min <- 
                                j.out.p.men <- 
                                j.out.p.max <- 
                                j.out.r.min <-
                                j.out.r.men <-
                                j.out.r.max <-
                                c()
                        } 
                        for (pa in 1:nrow(j.perms)) {
                                if (pa == 1){
                                        pa.out.p <- 
                                                pa.out.r <- 
                                                c()
                                }
                                        
                                pa.id <- which(
                                        j.cluster.assignments[[
                                                ceiling(evaluation / 25)
                                                ]] %in% 
                                                (j.perms[pa, ]
                                                 )
                                )
                                
                                pa.dist <-
                                        j.distance.matrix[[evaluation]] %>%
                                        as.matrix %>%
                                        .[pa.id, pa.id] %>% 
                                        as.dist
                                pa.ano <- anosim(
                                        x        = pa.dist,
                                        grouping = j.cluster.assignments[[
                                                ceiling(evaluation /25)]][pa.id],
                                        permutations = 100,
                                        parallel = 10
                                )
                                pa.out.r[pa] <- pa.ano$statistic 
                                pa.out.p[pa] <- pa.ano$signif
                                
                        }
                        j.out.p.min[evaluation] <- min(pa.out.p) 
                        j.out.p.men[evaluation] <- mean(pa.out.p) 
                        j.out.p.max[evaluation] <- max(pa.out.p) 
                        j.out.r.min[evaluation] <- min(pa.out.r)
                        j.out.r.men[evaluation] <- mean(pa.out.r)
                        j.out.r.max[evaluation] <- max(pa.out.r)
                }
                 
                j.out.p.min %<>% render_table(variable = "pairwise anosim minimum p value")
                j.out.p.men %<>% render_table(variable = "pairwise anosim mean p value")
                j.out.p.max %<>% render_table(variable = "pairwise anosim maximum p value")
                j.out.r.min %<>% render_table(variable = "pairwise anosim minimum R value")
                j.out.r.men %<>% render_table(variable = "pairwise anosim mean R value")
                j.out.r.max %<>% render_table(variable = "pairwise anosim maximum R value")
                
                #  ——— compute ANOSIM 
                # j.anosim  <- vector(mode = "numeric", length = neval)
                # j.anosim2 <- vector(mode = "numeric", length = neval)
                # for (anosim in 1:neval){
                #         out <- anosim(x            = j.distance.matrix[[anosim]],
                #                       grouping     = j.cluster.assignments[[ceiling(anosim/25)]],
                #                       permutations = 500,
                #                       parallel     = 6)
                #         j.anosim[anosim] <- out$statistic 
                #         j.anosim[anosim] <- out$signif 
                # }
                # j.anosim  %<>% render_table(variable = "anosim R")
                # j.anosim2 %<>% render_table(variable = "anosim p")
                
                # ——— compute area under the zeta diversity decline curve ————————————————
                j.ut  <- 1:j.types
                # j.ut2 <- table(j.tv)
                # if (any(j.ut2 < 11)){
                #         j.ut <- j.ut[-which(j.ut2<11)]
                # }
                j.zeta <- vector(mode = "list", length = j.types)
                #- LOOP over unique types. 
                #- establish inter-type AUCζ as a baseline 
                baseline <- vector(mode = "numeric", length = neval)
                for (zeta in 1:neval){
                        
                        z.types <- j.cluster.assignments[[ceiling(zeta/(neval/max.q))]]
                        # 10 repetitions for each envaluation 
                        for (zeta2 in 1:10){
                                if (zeta2 == 1) baseline2 <- vector(mode = "numeric", length = 10)
                                
                                z2id <- prop_sample(z.types, mean(table(z.types)))
                                z2data <- j.storage2[[zeta]][z2id, ]
                                baseline2[zeta2] <- Zeta.decline.ex(
                                        data.spec  = z2data,
                                        orders     = 1:v.orders,
                                        plot       = FALSE,
                                        rescale    = TRUE)$zeta.val %>% 
                                        calculate_auc(x = 1:v.orders)
                                
                                if (zeta2 == 10) baseline[zeta]<-mean(baseline2)
                                
                        }
                }
                
                #- One AUCζ value is computed for each type.  
                for (k in seq_along(j.ut)){
                        
                        #- prepare vector to store AUCζ for this iteration 
                        k.out <- vector(mode = "numeric", length = neval)
                        #-
                        for (zeta in 1:neval){
                                k.out[[zeta]] <- 
                                        Zeta.decline.ex(
                                                data.spec = j.storage2[[zeta]][which(j.cluster.assignments[[ceiling(zeta/(neval/max.q))]] == j.ut[k]),],
                                                orders     = 1:v.orders,
                                                plot       = FALSE,
                                                rescale    = TRUE)$zeta.val %>% 
                                        calculate_auc(x = 1:v.orders)
                        }
                        k.out2 <- k.out/baseline
                        j.zeta[[k]] <- render_table(k.out2, "auc_zeta")
                        rm(list = ls()[grepl("^k\\.", x = ls())])
                }
                rm(k)
                j.zeta %<>% rbindlist
                j.zeta[, value := mean(value), by = c("simulation", "q", "s", "k")]
                j.zeta %<>% unique(by = c("simulation", "q", "s", "k"))
                j.zeta[, type := NULL]
                
                #——— Compute PERMANOVA ———————————————————————————————————————————————————
                j.r <- j.p  <- j.F  <- 
                        vector(mode = "numeric", length = neval)
                for (prmnv in 1:neval){
                        out <- adonis2(formula = j.distance.matrix[[prmnv]] ~ j.cluster.assignments[[ceiling(prmnv/(neval/max.q))]])
                        j.r[prmnv] <- out$R2[1]
                        j.F[prmnv]  <- out$F[1]
                        j.p[prmnv]  <- out$`Pr(>F)`[1]
                }
                rm(out)
                rm(prmnv)
                j.r %<>% render_table("PERMANOVA R2")
                j.p %<>% render_table("PERMANOVA p")
                j.F %<>% render_table("PERMANOVA F")
                
                #———— COMPUTE Silhouette Width ———————————————————————————————————————————
                j.asw <- vector(mode = "numeric", length = neval)
                for (asw in 1:neval){
                        out <- silhouette(dist = j.distance.matrix[[asw]], 
                                          x = j.cluster.assignments[[ceiling(asw/(neval/max.q))]])
                        j.asw[asw] <- mean(out[, "sil_width"])
                }
                j.asw %<>% render_table("asw")

                #———— COMPUTE Indicator Value ————————————————————————————————————————————
                j.isa1 <- j.isa2 <- vector(mode = "numeric", length = neval)
                for (lp in 1:neval){
                        
                        out <- multipatt(j.storage2[[lp]], 
                                         j.cluster.assignments[[ceiling(lp/(neval/max.q))]], 
                                         func = "indval.g",
                                         control = how(nperm = 500))
                        out$sign$holm_p_value <- p.adjust(out$sign$p.value, method = "holm")
                        out1 <- nrow(out$sign[which(out$sign$holm_p_value<=0.05), ])/nrow(out$sign)
                        out2 <- mean(out$sign$holm_p_value, na.rm = T)
                        j.isa1[lp] <- out1
                        j.isa2[lp] <- out2
                }
                rm(out, out1, out2, lp)
                j.isa1 %<>% render_table("isa_number")
                j.isa2 %<>% render_table("isa_avg_p")
                
                #——— COMPUTE ISAMIC ——————————————————————————————————————————————————————
                j.isamic <- vector(mode = "numeric", length = neval)
                for (lp in 1:neval){
                        out <- isamic(
                                comm = j.storage2[[lp]], 
                                clustering = j.cluster.assignments[[ceiling(lp/(neval/max.q))]])
                        j.isamic[lp] <- mean(out)
                }
                rm(out, lp)
                j.isamic %<>% render_table("isamic")
                #  — — — Combine output into a single table  — — —
                j.out.final <- 
                        rbindlist(
                                list(
                                        j.cs,
                                        j.zeta,
                                        j.anosim,
                                        j.anosim2,
                                        j.r,
                                        j.p,
                                        j.F,
                                        j.asw,
                                        j.isa1,
                                        j.isa2,
                                        j.isamic
                                )
                        ) %>% 
                        .[, `:=` (run = i, types = j.types)]                
                i.out[[j]] <- j.out.final
                rm(list = ls()[grepl("^j\\.", x = ls())])
        } # END OF LOOP j over number of types 
        rm(j)
        
        #- save results of i loop to file 
        #i.model.number <- ifelse(i < 10, paste0("0",i), i)
        saveRDS(i.out, paste0("data/evaluations/eval_", i.model.name))
        #- Clean environment from i loop variables 
        rm(list = ls()[grepl("^i\\.", x = ls())])
        
} # END OF LOOP i 