library(dplyr)
library(testthat)
library(tidyr)

# load simulation results
file_path <- system.file("extdata", "simulation_pvalues.rds", package = "perturbplan")
simulation_pvalues <- readRDS(file_path)

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

# test when there is no signal
test_that("compute_power_posthoc aligns well with the simulation results (no cutoff)!", {

  # specify alpha
  alpha <- 1e-3
  simulation_power <- simulation_pvalues |>
    dplyr::group_by(run_id) |>
    dplyr::mutate(
      rejection = tidyr::replace_na((p.adjust(p_value, method = "BH") <= alpha), FALSE),
      QC = is.na(p_value)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(gene, grna) |>
    dplyr::summarise(
      sceptre_power = mean(rejection),
      QC_sceptre = mean(QC)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(pair = sprintf("%s_%s", grna, gene))


  compute_power_posthoc_results <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change_mean = effect_size_mean,
    fold_change_sd = effect_size_sd,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7,
    side = "left",
    multiple_testing_alpha = alpha
  )

  # merge two dfs
  merged_df <- compute_power_posthoc_results$individual_power |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::left_join(simulation_power, "pair")
  # compute the absolute difference
  abs_diff <- merged_df |>
    dplyr::mutate(abs_diff = abs(power - sceptre_power)) |>
    dplyr::select(abs_diff) |> dplyr::pull() |> unique()
  expect_lt(max(abs_diff),
            0.12)
})

# # plot scatter plot
# merged_df |>
#   ggplot2::ggplot(ggplot2::aes(x = sceptre_power, y = power)) +
#   ggplot2::geom_point() + ggplot2::geom_abline() +
#   ggplot2::theme_bw() +
#   ggplot2::labs(title = "With QC")


# test when there is no signal
test_that("compute_power_posthoc aligns well with the simulation results (with cutoff)!", {

  # specify alpha
  cutoff <- 1e-3
  simulation_power <- simulation_pvalues |>
    dplyr::group_by(run_id) |>
    dplyr::mutate(
      rejection = tidyr::replace_na(p_value <= cutoff, FALSE),
      QC = is.na(p_value)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(gene, grna) |>
    dplyr::summarise(
      sceptre_power = mean(rejection),
      QC_sceptre = mean(QC)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(pair = sprintf("%s_%s", grna, gene))


  compute_power_posthoc_results <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = "complement",
    fold_change_mean = effect_size_mean,
    fold_change_sd = effect_size_sd,
    num_total_cells = num_total_cell,
    n_nonzero_trt_thresh = 7,
    n_nonzero_cntrl_thresh = 7,
    side = "left",
    cutoff = cutoff
  )

  # merge two dfs
  merged_df <- compute_power_posthoc_results$individual_power |>
    dplyr::mutate(pair = paste(grna_target, response_id, sep = "_")) |>
    dplyr::left_join(simulation_power, "pair")
  # compute the absolute difference
  abs_diff <- merged_df |>
    dplyr::mutate(abs_diff = abs(power - sceptre_power)) |>
    dplyr::select(abs_diff) |> dplyr::pull() |> unique()
  expect_lt(max(abs_diff),
            0.09)
})

