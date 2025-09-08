process_geweke_results <- function(geweke_list) {
        # Create an empty list to store our processed results
        processed_results <- list()
        
        # Loop through each chain in the geweke results
        for(chain_num in seq_along(geweke_list)) {
                # Extract the z-scores for this chain
                chain_scores <- unlist(geweke_list[[chain_num]])
                chain_scores <- chain_scores[1:(length(chain_scores)-2)]
                # Create a data frame for this chain's results
                chain_df <- data.table(
                        parameter = names(chain_scores),
                        z_score = unname(chain_scores),
                        status = cut(abs(unname(chain_scores)),
                                     breaks = c(-Inf, 1.5, 2, Inf),
                                     labels = c("Good", "Marginal", "Poor")),
                        chain = paste("Chain", chain_num)
                )
                
                # Add to our processed results
                processed_results[[chain_num]] <- chain_df
        }
        
        # Combine all chains into one data frame
        final_df <- rbindlist(processed_results)
        
        # Add row names
        
        return(final_df)
}

plot_geweke_results <- function(geweke_df) {
        # Create the plot
        p <- ggplot(geweke_df, aes(x = parameter, y = z_score, color = status)) +
                geom_point(size = 3) +
                geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red", alpha = 0.5) +
                geom_hline(yintercept = c(-1.5, 1.5), linetype = "dashed", color = "orange", alpha = 0.5) +
                geom_hline(yintercept = 0, linetype = "solid", color = "gray", alpha = 0.5) +
                facet_wrap(~chain, scales = "free") +
                scale_color_manual(values = c("Good" = "forestgreen", 
                                              "Marginal" = "orange", 
                                              "Poor" = "red")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(title = "Geweke Diagnostic Results",
                     subtitle = "Dashed lines at ±1.5 (orange) and ±2 (red)",
                     x = "Parameter",
                     y = "Z-score",
                     color = "Convergence Status")
        
        return(p)
}
ge2 %>% filter(status != "Good") %>% plot_geweke_results()
