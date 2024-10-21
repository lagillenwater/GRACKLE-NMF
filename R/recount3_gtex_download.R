## Load required libraries 
source("./R/functions.R")
library(recount3)

human_projects <- available_projects() # load all human projects

gtex_projects <- human_projects %>% filter(file_source == "gtex")


project_data <- lapply(1:nrow(gtex_projects), function(x) create_rse(gtex_projects[x,]))


metadata <- lapply(project_data, colData)

gtex_metadata <- lapply(metadata, function(x) {
  x %>% as.data.frame() %>%select(contains("gtex"))
})




gene_metadata <- lapply(project_data, rowRanges)

## filter by protein coding genes
protein_coding_genes <- lapply(gene_metadata, function(x) {
  x %>% as.data.frame() %>% filter(gene_type == "protein_coding")
})


expression <- lapply(project_data,assay)

protein_coding_expression <- lapply(expression, function(x) {
  x %>% as.data.frame() %>% filter(rownames(.) %in% protein_coding_genes[[1]]$gene_id)
})




inner_join_by_rownames <- function(mat_list) {
  common_rownames <- Reduce(intersect, lapply(mat_list, rownames))
  joined_matrix <- do.call(cbind, lapply(mat_list, function(mat) mat[common_rownames, , drop = FALSE]))
  return(joined_matrix)
}

gtex_expression <- inner_join_by_rownames(protein_coding_expression)

gtex_expression_min_max <- apply(gtex_expression,1,min_max_scale)

mean_absolute_deviation <- function(x) {
  mean(abs(x - mean(x)))
}



gtex_mad <- apply(gtex_expression_min_max,1,mean_absolute_deviation)
mad_rank <- gtex_mad[order(gtex_mad, decreasing = TRUE)]



# 
# project <- human_source %>% filter(project == opt$project)
# 
# rse_gene = create_rse(project)
# #assay(rse_gene, "counts") = transform_counts(rse_gene)
# #assays(rse_gene)$RPKM = recount::getRPKM(rse_gene)
# expression =as.data.frame(t(assays(rse_gene)$RPKM))
# metadata =  tryCatch({                               
#   expand_sra_attributes(rse_gene)
# }, error = function(e) {})
# 
# metadata = colData(metadata)
# 
# ## write expression and metadata
# write.csv(expression, file = paste0(opt$out_dir,project$project,"_expression.csv"))
# write.csv(metadata, file = paste0(opt$out_dir, project$project,"_metadata.csv"))
# 
