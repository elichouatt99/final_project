library(ggplot2)

# load data
phenotypes <- readRDS("datasets/traits/phenotypes.rds")

# scatter of motor10 vs. decline
phenotypes.names <- colnames(phenotypes)

for (pheno1 in phenotypes.names){
  for (pheno2 in phenotypes.names){
    if (pheno1 < pheno2){
      x_limits <- c(min(phenotypes[[pheno1]],na.rm = TRUE), max(phenotypes[[pheno1]],na.rm = TRUE))
      y_limits <- c(min(phenotypes[[pheno2]],na.rm = TRUE), max(phenotypes[[pheno2]],na.rm = TRUE))
      scatterplot <- ggplot(phenotypes, aes(x=!!sym(pheno1), y=!!sym(pheno2))) +
        geom_point() + 
        lims(x = x_limits, y = y_limits) +
        labs(x = pheno1, y = pheno2) + # Label axes
        ggtitle(paste0("Scatterplot of ",pheno1 ," vs. ", pheno2, "(cor: ", round(cor(phenotypes[[pheno1]],phenotypes[[pheno2]],use = "complete.obs"), 2), ")"))
      
      # Save the scatterplot as a file (e.g., in PNG format)
      ggsave(paste0("plots/transcriptomic_analysis/scatterplots/", pheno1, '_vs_', pheno2, '.png'), plot = scatterplot, width = 6, height = 6, dpi = 300)
      
    }
  }
}

complete_cases <- na.omit(phenotypes[c("cogn_global_lv", "motor10_lv")])
cor.test(complete_cases$cogn_global_lv, complete_cases$motor10_lv)

# Violin plots of different traits: AB, Tau, TDP43
for (pheno in phenotypes.names){
    violinplot <- ggplot(phenotypes, aes(x = pheno, y = !!sym(pheno))) +
      geom_violin(fill="lightblue") +
      labs(x = "Trait", y = " Value") +
      ggtitle(paste("Violin Plot of", pheno))
      ggsave(paste0("plots/transcriptomic_analysis/violinplots/", "violin_plot", pheno, '.png'), plot = violinplot, width = 6, height = 6, dpi = 300)
}
