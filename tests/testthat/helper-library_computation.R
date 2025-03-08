
# load gene expression summary, James' results and size factor
file_path <- system.file("extdata", "gene_expression_summary.rds", package = "perturbplan")
gene_expression <- readRDS(file_path)

# number of treatment and control cell in their resupt
num_element <- 14
num_total_cell <- 56000
treatment_size <- c(50, 75, 100, 150, 250, 500, 750, 1000,
                    1400, 1800, 2500, 4000, 7500, 10000)

# specify the number of gRNAs in each group (15)
K <- 15
trt_mat <- sapply(treatment_size, function(x){
  rep(round(x / K), K)
})
target_cell_df <- data.frame(
  gRNA_id = rep(sprintf("gRNA%d", 1:K), num_element),
  gRNA_target = rep(sprintf("enh%d", 1:num_element), each = K),
  num_cells = as.vector(trt_mat)
)
control_cell_vec <- setNames(as.vector(num_total_cell - treatment_size),
                             sprintf("enh%d", 1:num_element))

# extract the gene expression iformation: mean and size parameter
num_gene <- nrow(gene_expression)
baseline_expression <- setNames(gene_expression$mean, gene_expression$ensembl)
size_parameter <- setNames(1 / gene_expression$dispersion, gene_expression$ensembl)

# provide effect size matrix and effect_size_sd matrix
effect_size_mean <- 0.85
effect_size_sd <- 0.13

############## construct dfs specific to compute_power_posthoc #################
# construct discovery pairs
discovery_pairs <- expand.grid(
  response_id = gene_expression$ensembl,
  grna_target = unique(target_cell_df$gRNA_target)
)

# construct cells_per_grna
cells_per_grna <- target_cell_df |>
  dplyr::rename(grna_target = gRNA_target,
                grna_id = gRNA_id)

# construct baseline_expression_stats
baseline_expression_stats <- data.frame(
  response_id = gene_expression$ensembl,
  expression_mean = baseline_expression,
  expression_size = size_parameter
)
