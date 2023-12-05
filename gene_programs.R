library(dplyr)

dir.create("plots/gene_programs", showWarnings = FALSE)

# number of genes per modules (per tissue ?) comparaison entre l4 et l5


dir.create("plots/gene_programs/histograms", showWarnings = FALSE)


#' compute histogram of the size of the modules
#'
#' @param data MEs
#' @param level level chosen in the Speak Easy algorithm cluster tree
#'
#' @return save histogram
size_hist <- function(data, level){
  dir.create(paste0("plots/gene_programs/histograms/", level), showWarnings = FALSE)

  total_count <- data_frame(module=NULL, Number_of_Genes=NULL)
  
  # per tissue
  tissues <- c("Brain", "SpinalCord", "Muscle")
  for (tissue in tissues){
    module_counts <- data[[tissue]]$genes %>% distinct(gene, module) %>%
      group_by(module) %>%
      summarize(Number_of_Genes = n())
    total_count <- rbind(total_count, module_counts)
    filename <- paste0("plots/gene_programs/histograms/",level, "/histogram_size_", tissue, ".png")
    png(filename, width = 800, height = 600)
    hist(module_counts$Number_of_Genes, main = paste("Histogram of Module Sizes for", tissue, "at", level), xlab = "Number of Genes", ylab = "Frequency")
    dev.off()
  }
  
  # overall
  filename <- paste0("plots/gene_programs/histograms/",level, "/histogram_size_all_tissues.png")
  png(filename, width = 800, height = 600)
  hist(total_count$Number_of_Genes, main = paste("Histogram of Module Sizes at", level), xlab = "Number of Genes", ylab = "Frequency")
  dev.off()
}


#' compute histogram of the correlation between the modules
#'
#' @param data MEs
#' @param level level chosen in the Speak Easy algorithm cluster tree
#'
#' @return save histogram
cor_hist <- function(data, level){
  dir.create(paste0("plots/gene_programs/histograms/", level), showWarnings = FALSE)
  
  total_cor <- data_frame(module1=NULL, module2=NULL, correlation_value=NULL)
  
  # per tissue
  tissues <- c("Brain", "SpinalCord", "Muscle")
  for (tissue in tissues){
    # Créez un data frame avec la matrice de corrélation
    mat.cor <- cor(data[[tissue]]$MEs)
    modules_names <- colnames(data[[tissue]]$MEs)
    correlation_df <- data.frame(module1 = rep(modules_names, each = length(modules_names)),
                                 module2 = rep(modules_names, times = length(modules_names)),
                                 correlation_value = as.vector(mat.cor))
    
    # Filtrez pour avoir module1 < module2
    correlation_df <- correlation_df %>%
      filter(module1 < module2)
    
    total_cor <- rbind(total_cor, correlation_df)
    
    filename <- paste0("plots/gene_programs/histograms/",level, "/histogram_cor_", tissue, ".png")
    png(filename, width = 800, height = 600)
    hist(correlation_df$correlation_value, main = paste("Histogram of Correlation Values for", tissue, "at", level), xlab = "Number of Pairs of different Modules", ylab = "Frequency")
    dev.off()
  }
  # overall
  filename <- paste0("plots/gene_programs/histograms/",level, "/histogram_cor_all_tissues.png")
  png(filename, width = 800, height = 600)
  hist(total_cor$correlation_value, main = paste("Histogram of all Module Correlation at ", level), xlab = "Number of Pairs of different Modules", ylab = "Frequency")
  dev.off()
}

modules_data_l4 <- readRDS("datasets/modules/modules_data1.rds")
size_hist(modules_data_l4, "l4")
cor_hist(modules_data_l4, "l4")
tissues <- c("Brain", "SpinalCord", "Muscle")
for (tissue in tissues){
  print(paste(tissue, ":", length(colnames(modules_data_l4[[tissue]]$MEs))))
}

modules_data_l5<- readRDS("datasets/modules/modules_data_l5.rds")
size_hist(modules_data_l5, "l5")
