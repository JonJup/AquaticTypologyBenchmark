# ————————————————————————— #
# ——— Visualize results ——— # 
# ————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 20.01.2025


# load data -------------------------------------------------------------------------
# d <- dir_ls("data/evaluations/")
# d <- lapply(d, readRDS)
# d <- lapply(d, rbindlist)
# d <- rbindlist(d)


d <- readRDS("data/eval_w_dirichlet.rds")

d$types %<>% factor
d$run   %<>% factor
d$rank_sum2 %<>% factor

ggplot(d, aes(x = value, group = rank_sum2)) +
        geom_density(aes(fill = rank_sum2), alpha = 0.5) + 
        facet_wrap(.~metric, scales = "free")


d %>%
        filter(metric == "cs") %>% 
ggplot(aes(y = value, x = env)) +
        geom_density_2d_filled() +
        scale_fill_viridis_d() + 
        facet_wrap(.~rank_sum2)
d %>% 
        filter(metric == "cs") %>%
        ggplot(aes(y = value, x = env_asw, 
                   weight = normalized_dirichlet_density)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = F) + 
        facet_wrap(.~rank_sum2, scales = "free") 

d %>% 
        ggplot(aes(x = rank_sum2, y = value, weight = normalized_dirichlet_density)) + 
        geom_boxplot() + 
        facet_wrap(.~metric, scales = "free")

d %>% 
        ggplot(aes(x = rank_sum2, y = value, weight = normalized_dirichlet_density)) + 
        #geom_jitter(width = .2) + 
        geom_line(stat="smooth",
                  method = "loess",
                  size = 1.5,
                  aes(group = run),
                  alpha = 0.1) +
        geom_line(stat = "smooth", 
                  method = "loess",
                  aes(group = 1),
                  col = "red", 
                  size = 2) + 
        facet_wrap(.~metric, scales = "free") 
d %>% 
        ggplot(aes(x = env_asw, y = value, weight = normalized_dirichlet_density)) + 
        #geom_jitter(width = .2) + 
        geom_line(stat="smooth",
                  method = "loess",
                  size = 1.5,
                  aes(group = run),
                  alpha = 0.1) +
        geom_line(stat = "smooth", 
                  method = "loess",
                  aes(group = 1),
                  col = "red", 
                  size = 2) + 
        facet_wrap(.~metric, scales = "free") 
   
d %>% 
        ggplot(aes(y = env_asw, x = contraction_points)) + 
        geom_point(alpha = 0.5) + 
        geom_smooth()
d %>% 
        ggplot(aes(y = env_asw, x = contraction_centroids)) + 
        geom_point(alpha = 0.5) + 
        geom_smooth()
     

library(ggridges)
ggplot(d, aes(x = value, y = rank_sum2, fill = factor(rank_sum2), 
              weight = normalized_dirichlet_density)) +
        geom_density_ridges(alpha = 0.6, scale = 0.9, aes(fill = rank_sum2)) +
        scale_fill_discrete() +
        theme_ridges() +
        labs(x = "value", y = "") +
        theme(legend.position = "right") + 
        facet_wrap(.~metric, scales = "free")



d[, mean_ntype := mean(x=value), by = "types"]
d[, mean_run   := mean(value), by = "run"]

overall_run_mean <- weighted.mean(x = unique(d, by = "run")$value, w= unique(d, by = "run")$normalized_dirichlet_density)

ggplot(d, aes(y = value, x =types)) +
        geom_violin(draw_quantiles = .5) + 
        geom_jitter(width = 0.4, aes(alpha = normalized_dirichlet_density))



ggplot(d, aes(y = mean_ntype, x = types)) +
        geom_violin(draw_quantiles = .5) + 
        geom_jitter(width = 0.4, aes(alpha = normalized_dirichlet_density))

d %>% unique(by="run") %>%
        mutate(run = fct_reorder(run, value)) %>%
        ggplot(aes(y = value, x = run)) +
        geom_point(aes(color = normalized_dirichlet_density)) + 
        geom_hline(data = overall_run_mean, yintercept = overall_run_mean) + 
        viridis::scale_colour_viridis()


d %>% 
        ggplot(aes(x =scale_max, y =  space)) + 
        geom_point()


# ternary plots  --------------------------------------------------------------------

library(ggtern)

# Create your d
ggtern(data = d, 
       aes(x = env, y = bio, z = space, col = stochastic)) +
        geom_point(size = 2) +
        theme_showarrows()+ 
        #theme_bw() +
        labs(x = "env",
             y = "bio", 
             z = "space") + 
        scale_color_viridis_b()

# newer plots -----------------------------------------------------------------------

d %>%
        ggplot(aes(x = contraction_centroids, y = value)) + 
        geom_jitter(width = 0.1, shape = 21) + 
        geom_smooth(se = FALSE) + 
        facet_wrap(.~metric, scales = "free")
d %>%
        ggplot(aes(x = types, y = value, fill = run)) + 
        geom_boxplot() + 
        facet_wrap(.~metric, scales = "free")
d %>%
        ggplot(aes(x = env_asw, y = value, fill = run)) + 
        #geom_point(alpha = 0.3) + 
        geom_smooth(method = "lm", se = F, alpha = .4) + 
        facet_wrap(.~metric, scales = "free")


# against ENV ------------------------------------------------------------
d %>% unique(by="run") %>%
        ggplot(aes(y = value, x = env)) +
        geom_smooth(method = "lm", se =F) +
        geom_point(aes(color = normalized_dirichlet_density)) + 
        geom_hline(data = overall_run_mean, yintercept = overall_run_mean) + 
        viridis::scale_colour_viridis()
# against SPACE ------------------------------------------------------------
d %>% unique(by="run") %>%
        ggplot(aes(y = value, x = space)) +
        geom_smooth(method = "lm", se =F) + 
        geom_point(aes(color = normalized_dirichlet_density)) + 
        geom_hline(data = overall_run_mean, yintercept = overall_run_mean) + 
        viridis::scale_colour_viridis()
# against BIO ------------------------------------------------------------
d %>% unique(by="run") %>%
        ggplot(aes(y = value, x = bio)) +
        geom_smooth(method = "lm", se =F) + 
        geom_point(aes(color = normalized_dirichlet_density)) + 
        geom_hline(data = overall_run_mean, yintercept = overall_run_mean) + 
        viridis::scale_colour_viridis()
# against STOCHASTICITY ------------------------------------------------------------
d %>% unique(by="run") %>%
        ggplot(aes(y = value, x = stochastic)) +
        geom_smooth(method = "lm", se =F) + 
        geom_point(aes(color = normalized_dirichlet_density)) + 
        geom_hline(data = overall_run_mean, yintercept = overall_run_mean) + 
        viridis::scale_colour_viridis()

