#' For a given tissue, compute the Correlation matrix with annotation of the 
#' clusters.
#' 
#'
#' @param tissue tissue to analyse
#' @param num_clusters final number of clusters
#' @param cut_height cut height of the cluster tree (dendogram)
#'
#' @return Save heatmap correlation matrix with annotation of the clusters
cluster_heatmap <- function(tissue,traits=c("cogng_demog_slope"), num_clusters=5, cut_height=4.5){
  MEs <- load.data(tissue)
  matrix_cor <- qgraph::cor_auto(MEs)

  clusters <- get_clusters(matrix_cor)
  cluster_membership <- data.frame(row_names = rownames(matrix_cor), cluster = clusters)
  
  heatmap_object <- Heatmap(
    matrix_cor,
    name = "Correlation",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_column_dend = F
  )
  
  # Generate dynamic colors for clusters
  cluster_colors <- rainbow(length(unique(clusters)))
  
  cluster_color_list <- setNames(cluster_colors, unique(as.character(clusters)))
  
  # Create clusters row annotations
  cluster_ha <- rowAnnotation(clusters = as.character(clusters), 
                                  col=list(clusters=cluster_color_list), border=T)
  
  # Create sig row annotation
  sig_ha <- NULL
  for (trait in traits){
    SE <- readRDS(paste0("datasets/traits/SE_", tissue,"_l4.rds"))
    beta.sign <- ifelse(is.na(SE$cor[colnames(MEs), trait]), 1, sign(SE$cor[colnames(MEs), trait]))
    fdr <- ifelse(is.na(SE$pval.signed[colnames(MEs), trait]), 1, SE$pval.signed[colnames(MEs), trait])
    pval <- 10^(-fdr * beta.sign)
    binary_sig <- ifelse( pval < 0.05, "**", ifelse(pval < 0.1, "*", "not sig"))
    binary_sig[is.na(binary_sig)] <- "not sig"
    row_ha <- rowAnnotation(cogng = binary_sig, col = list(cogng = c("not sig"= "white", "*" = "grey", "**" = "black")), border=T)
    label <- ifelse(trait == 'cogng_demog_slope', 'cogng', trait)
    row_ha@anno_list[["sig"]]@label = label
    sig_ha <- row_ha + sig_ha 
  }

  # Combine the heatmap and row annotations
  heat_annotation_combined <- sig_ha + cluster_ha + heatmap_object
  
  filename <- paste0("plots/link_to_traits/heatmaps/",tissue, ".png")
  png(filename,2048,2048)
  # Plot the combined heatmap with row annotations
  draw(heat_annotation_combined, heatmap_legend_side = "left")
  dev.off()
}


tissues <- c("Brain", "SpinalCord", "Muscle")
for (tissue in tissues){
  cluster_heatmap(tissue)
}


#' For a given tissue, compute the correlation matrix of its modules, its 
#' clusters, and for each one, compute the subnetwork with traits.
#'
#' @param tissue 
#' @param traits 
#' 
#' @return save the subnetworks and the corresponding adjacency matrix 
#' @export
#'
#' @examples
cluster_network <- function(tissue, traits=c("cogng_demog_slope")){
  MEs <- load.data(tissue)
  matrix_cor <- qgraph::cor_auto(MEs)
  clusters <- get_clusters(matrix_cor)
  cluster_membership <- data.frame(row_names = rownames(matrix_cor), cluster = clusters)
  for (cluster_name in unique(clusters)){
    Results <- compute_network(tissue,traits = traits, cluster=cluster_membership[cluster_membership$cluster == cluster_name, ]$row_names, threshold = T, penalize=T)
    dir.create("plots/link_to_traits/network/", showWarnings = FALSE)
    dir.create(paste0("plots/link_to_traits/network/", tissue), showWarnings = FALSE)
    filename <- paste0("plots/link_to_traits/network/", tissue, "/",cluster_name, ".pdf")
    pdf(filename)
    plot(Results$Network)
    draw(Heatmap(Results$Network$graph))
    dev.off()
  }
}

tissues <- c("Brain", "SpinalCord", "Muscle")
for (tissue in tissues){
  cluster_network(tissue)
}
