# Render table 
render_table <- function(x, variable){
        x2  <- data.table(
                        scheme_id             = i.ss$scheme_id,
                        value                 = x, 
                        #model                 = model_names[i],
                        n_types               = i.file$number_of_clusters,
                        metric                = variable, 
                
                        variables             = rep(i.file$number_of_variables, each = neval/max.q),
                        contraction_points    = i.file$contraction_points,
                        contraction_centroids = i.file$contraction_centroids,
                        env_asw               = i.file$asw,
                        importance            = rep(i.file$variable_importance, each = neval/max.q)
                )
        return(x2)
}

## alternative version; copied from setup script on the 14.08.25
# render_table <- function(x, variable){
#         x2  <- data.table(
#                 value                 = x, 
#                 simulation            = rep(1:max.ss, times = neval/max.ss),
#                 q                     = rep(1:max.q, each  = neval/max.q),
#                 cluster_adjustment    = rep(1:within.q, times = max.q),
#                 metric                = variable, 
#                 variables             = rep(j.n.variables, each = neval/max.q),
#                 contraction_points    = j.contraction.points,
#                 contraction_centroids = j.contraction.centroids,
#                 env_asw               = j.env.asw,
#                 importance            = rep(j.importance, each = neval/max.q)
#         )
#         return(x2)
# }
# return(render_table)