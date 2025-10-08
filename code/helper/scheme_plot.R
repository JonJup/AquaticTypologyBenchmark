setwd(rstudioapi::getActiveProject())

library(data.table)

files <- list.files("data/biota/", pattern = "03", full.names = T)
files <- lapply(files, readRDS)
files <- rbindlist(files, fill = T)

library(ggplot2)

taxon_counts <- as.data.frame(table(files$taxon))
colnames(taxon_counts) <- c("taxon", "count")
taxon_labs <- setNames(
        paste0(taxon_counts$taxon, " (n=", taxon_counts$count, ")"),
        taxon_counts$taxon
)

files$organismQuantityType <- factor(files$organismQuantityType)
ggplot(files, aes(x = samples, n_taxa)) + 
        geom_jitter(width = 20, height = 10, 
                    aes(fill = organismQuantityType ), 
                    shape = 21, alpha = 0.7, size = 3)  + 
        facet_wrap(.~taxon, scales = "free", labeller = labeller(taxon = taxon_labs)) + 
        scale_fill_brewer(palette = "Set2") +
        labs(
                x = "Samples",
                y = "Number of Taxa",
                title = "Distribution of Taxa Counts by Sample",
                subtitle = "Faceted by Taxon (n = count in labels)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
                strip.text = element_text(face = "bold", size = 12),
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "bottom",
                legend.title = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_rect(fill = "#fffaf3"),
                panel.background = element_rect(fill = "#fcfbf8ff")
        )
ggsave("output/figures/scheme_plot.tiff")
