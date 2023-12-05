library(ggplot2)
library(ComplexHeatmap)
library(tidyr)

# Load the data 
load.dataset <- function(tissue, level){
  infos <- read.csv("datasets/samples_infos.csv")
  path <- paste0("datasets/partial_cor/", tissue, "_modules.rds")
  dataset <- readRDS(path)
  dataset$MEs <- dataset$MEs[,grepl(paste0("^", level), colnames(dataset$MEs))]
  rownames(dataset$MEs) <- plyr::mapvalues(rownames(dataset$MEs), infos$SampleID_BatchID, infos$projid)
  return(dataset)
}

filter.length.MEs <- function(MEs, tissue){
  filtered.MEs <- MEs
  modules <- readRDS(paste0("datasets/partial_cor/",tissue,'_modules.rds'))
  for (module in colnames(MEs)){
    if (length(distinct(modules$genes[modules$genes$module == module,]$gene)) <30){
      filtered.MEs[[module]] <- NULL
    }
  }
  return(filtered.MEs)
}

tissue <- "Brain"
level <- "l4"
MEs <- load.MEs(tissue, level)
phenotypes <- readRDS("datasets/traits/phenotypes.rds")
traits <- c("cogng_demog_slope", "motor10_lv", "tdp_st4", "amylsqrt", "tangsqrt")
common_row_names <- intersect(rownames(MEs), rownames(phenotypes))

readROSMAPInfos <- function(){
  df <- read.csv("datasets/samples_infos.csv") %>% filter(keep.to.wgcna)
  #df$tissue <- unlist(lapply(strsplit(df$ExternalID, "_"), function(annotation){annotation[2]}))
  df$UNMAPPED_READS <- log10(df$UNMAPPED_READS + 1)
  df$PF_READS_ALIGNED <- log10(df$PF_READS_ALIGNED + 1)
  df$CODING_BASES <- log10(df$CODING_BASES + 1)
  df$amylsqrt <- as.numeric(df$amylsqrt)
  df$motor10_lv <- as.numeric(df$motor10_lv)
  df$BatchID <- as.character(paste0("Batch", df$BatchID))
  df$gait_speed_lv <- (as.numeric(replace(df$gait_speed_lv, ".", NA)))[1:nrow(df)]
  df$motor_gait_lv <- (as.numeric(replace(df$motor_gait_lv, ".", NA)))[1:nrow(df)]
  df$cogn_ep_lv <- (as.numeric(replace(df$cogn_ep_lv, ".", NA)))[1:nrow(df)]
  df$cogn_po_lv <- (as.numeric(replace(df$cogn_po_lv, ".", NA)))[1:nrow(df)]
  df$cogn_ps_lv <- (as.numeric(replace(df$cogn_ps_lv, ".", NA)))[1:nrow(df)]
  df$cogn_se_lv <- (as.numeric(replace(df$cogn_se_lv, ".", NA)))[1:nrow(df)]
  all.clinical <- read.csv("datasets/celmod.donors.metadata.csv")
  df$cogng_demog_slope <- -(plyr::mapvalues(df$projid, all.clinical$X, all.clinical$cogng_demog_slope) %>% as.numeric())
  return(df)
}

get_features.to.calculate <- function(vec = NULL){
  if(is.null(vec)){
    vec <- c(features.by.cat(c("motor", "pathologies", "cognitive", "msex", "metabolic")), "age_death")
  }
  return(vec)
}

features.by.cat <- function(categories = c("msex", "motor", "pathologies", "cognitive", "quality", "metabolic")){
  df.infos<- readROSMAPInfos()
  lev.motors <- colnames(df.infos)[endsWith(colnames(df.infos), "lv")]
  all_lst <- list()
  all_lst[["msex"]] <- c("msex")
  all_lst[["motor"]] <- c("parksc_lv", "motor_dexterity_lv",
                          "parkinsonism_yn_lv",
                          "motor10_lv", "tdp_st4",
                          "gait_speed_lv",
                          "motor_handstreng_lv", "motor_gait_lv"
                          )
  all_lst[["pathologies"]] <- c("gpath", 
                                "amylsqrt", "tangsqrt")
  all_lst[["metabolic"]] <- c("glucose", "hba1c", "vasc_risks_sum", "diabetes_sr_rx_ever")
  all_lst[["cognitive"]] <- c("cogn_global_lv", "cogdx", "cogng_demog_slope")
  all_lst[["quality"]] <- c("RQN", "pmi", 
                            "age_death", 
                            "MEDIAN_CV_COVERAGE", "UNMAPPED_READS", "MEDIAN_3PRIME_BIAS","CODING_BASES", "PCT_INTRONIC_BASES"
                            )
  #return(intersect(as.vector(unlist(lapply(categories, function(cat){all_lst[[cat]]}))), colnames(infos)))
  return(as.vector(unlist(lapply(categories, function(cat){all_lst[[cat]]}))))
}

get_phenotype <- function(features, samples, column.id){
  clinical.data <- readROSMAPInfos() %>% dplyr::select(c(column.id, features)) %>% unique() %>% `rownames<-`( NULL ) %>% tibble::column_to_rownames(column.id)
  phenotype <- clinical.data[samples,]
}

# Output a dataframe of the association (lm) of obs.m by traits.m
# Corrected by mediators
# Dataframe output is of the form: trait_, module_, association, pvalue
calculate.corrected.association <- function(obs.m, traits.m, mediators.m){
  vars.m <- cbind(traits.m[rownames(obs.m),], mediators.m[rownames(obs.m),])
  traits_mediator <- colnames(mediators.m)
  mediators <- paste(traits_mediator, collapse = " + ")
  df.m <- cbind(obs.m, vars.m)
  mat <- do.call(rbind,apply(expand.grid(setdiff(colnames(traits.m), colnames(mediators.m)),
                                         colnames(obs.m)), 1, function(params){
                                           trait <- as.vector(params[[1]])
                                           module <- as.vector(params[[2]])
                                           df.m.sub <- df.m[complete.cases(df.m[,setdiff(c(trait, module, traits_mediator), "")]),]
                                           f_mediators <- as.formula(paste(module, paste(trait, mediators, sep='+'), sep='~'))
                                           f_without <- as.formula(paste(module, trait, sep='~'))
                                           x <- lm(f_mediators, data=df.m.sub)$coefficients[[trait]]
                                           y <- lm(f_without, data=df.m.sub)$coefficients[[trait]]
                                           pvalue <- coef(summary(lm(f_mediators, data=df.m.sub)))[trait, "Pr(>|t|)"]
                                           #return(data.frame(trait_ = trait, module_ = module, association=x, not_corrected = y, pvalue = pvalue))
                                           return(data.frame(trait_ = trait, module_ = module, association=x, pvalue = pvalue))
                                         }))
  return(mat)
}

# Output a dataframe of the association (lm) of obs.m by traits.m
# Corrected by mediators
# Dataframe output is of the form: trait_, module_, association, pvalue
# Calculate correlations between clinical traits and MEs - correct by mediators
calculate.corrected.MEs_association <- function(data.df,
                                                traits_mediator,
                                                features.to.calculate = NULL, column.id = "projid"){
  features <- get_features.to.calculate(features.to.calculate)
  phenotype <- get_phenotype( unique(c(traits_mediator, features)), rownames(data.df$MEs), column.id)
  mat <- calculate.corrected.association(data.df$MEs, phenotype[,features], phenotype[,traits_mediator])
  return(mat)
}

pvalue_to_str <- function(vec){
  vec.str <-  ifelse(vec < 0.1, "*", "")
  vec.str <- ifelse(vec < 0.05, "**", vec.str)
  vec.str <- ifelse(vec < 0.01, "***", vec.str)
  return(vec.str)
}

#data.df = list(MEs = x);traits_mediator=c("age_death", "msex"); features.to.calculate =  setdiff(features.by.cat(), features.by.cat("quality")); column.id = "projid";return.df = T;association.by.cor=F
associations.traits.corrected.df.plots <- function(data.df,
                                                   traits_mediator = "",
                                                   set_definition = NULL,
                                                   return.df = F,
                                                   plot.color.pvalue = F, features.to.calculate = NULL, column.id = "projid", association.by.cor = F){
  # Calculate associations between MEs and clinical, correcting by mediators
  if(!association.by.cor){
    mat <- calculate.corrected.MEs_association(data.df, traits_mediator, features.to.calculate = features.to.calculate, column.id = column.id)
    print(mat)
  } else {
    mat <- calculate.correlation(data.df, merge.it = T, method = "pearson", features = features.to.calculate, column.id = column.id) %>% `colnames<-`(c("module_", "trait_", "pvalue", "association"))
  }
  if("infos" %in% names(data.df)){
    mat <- mat %>% merge(., data.df$infos[,c("module", "tissue")], by.x = "module_", by.y = "module")
  } else {
    mat <- mat %>% mutate(tissue = "tissue")
  }
  print(mat)
  # Correct for multiple hypothesis and add str
  mat <- mat %>% group_by(trait_, tissue) %>% mutate(fdr = p.adjust(pvalue, method = "fdr"),
                                                     BH = p.adjust(pvalue, method = "BH"),
                                                     bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% mutate(fdr.str = pvalue_to_str(fdr),
                                                                                                                      fdr.signed = -log10(fdr)*sign(association)) %>% arrange(tissue)%>% ungroup()%>% dplyr::select(-c("tissue"))
  print(mat)
  ms <- c(cor = "association", pval.signed = "fdr.signed", pval = "fdr.str")
  mats <- setNames(lapply(names(ms),function(m_){(mat %>% dplyr::select(c("trait_", "module_", ms[[m_]])) %>% spread(., trait_, ms[[m_]]) %>% tibble::column_to_rownames("module_"))}) , names(ms))
  print(mats)
  if(return.df){
    return(mats)
  }
  modules.annots <- data.df$infos[,c(set_definition, "n.genes", "module")] %>% tibble::column_to_rownames("module")
  colors.cor <- color.pheatmap.neg.pos(mats$cor, colors = c("#003BC0", "#FFFFFF", "#D80D29"))
  ann_colors <-  list(tissue = annotations.colors()$tissue)
  mediators_labs <- paste(traits_mediator, collapse = " + ")
  if (!plot.color.pvalue){
    return(cowplot::plot_grid(pheatmap::pheatmap(mats$cor, display_numbers = mats$pval ,
                                                 annotation_row = modules.annots,
                                                 breaks = colors.cor$breaks,
                                                 color = colors.cor$colors,
                                                 annotation_colors = ann_colors, cluster_rows = F, main = mediators_labs, silent= T)$gtable,
                              pheatmap::pheatmap(mats$cor, display_numbers = mats$pval, annotation_row = modules.annots,
                                                 breaks = colors.cor$breaks, color = colors.cor$colors,
                                                 annotation_colors = ann_colors, main = mediators_labs, silent= T)$gtable))
  } else {
    return(cowplot::plot_grid(pheatmap::pheatmap(cor(data.df$MEs, use = "pairwise.complete.obs"), clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", silent = T,
                                                 breaks = color.pheatmap.corr()$breaks,
                                                 color = color.pheatmap.corr()$colors, legend = F)$gtable,
                              pheatmap::pheatmap(mats$pval.signed, display_numbers = mats$pval , annotation_row = modules.annots,
                                                 breaks = colors.cor$breaks, color = colors.cor$colors,
                                                 annotation_colors = ann_colors, cluster_rows = F, main = mediators_labs, silent= T)$gtable,
                              pheatmap::pheatmap(mats$cor, display_numbers = mats$pval , annotation_row = modules.annots,
                                                 breaks = colors.cor$breaks, color = colors.cor$colors,
                                                 annotation_colors = ann_colors, cluster_rows = F, main = mediators_labs, silent= T, scale = "column")$gtable, nrow = 1))
  }
}

mediators <- c("msex", "age_death")
for (tissue in tissues){
  data.df <- load.dataset(tissue, "l4")
  associations.traits.corrected.df.plots(data.df, mediators, features.to.calculate = traits)
}
