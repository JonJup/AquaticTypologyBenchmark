

d <- readRDS("data/eval_w_dirichlet.rds")
d$types <- factor(d$types)
d$run <- factor(d$run)

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

