# This is a R script generating data from James' simulation setup


#' This is a function generating data from James' DGP
#'
#' @param effect_size Effect size: a scalar
#' @param guide_sd Guide RNA sd
#' @param discovery_pairs Discovery pairs
#' @param sce sce object
#' @param grna_target_data_frame grna target dataframe
#' @param trt_group treatment group
#' @param size_factor A logic value; if TRUE, use the one in sce; false, set to 1
#'
#' @return A lsit of information
#' @export
#' @importFrom dplyr filter pull
#' @importFrom stats setNames


generate_James_data <- function(effect_size, guide_sd, discovery_pairs,
                                sce, grna_target_data_frame, trt_group, size_factor = TRUE){

  # reset the cell size factor if size_factor is false
  if(!size_factor){

    # extract the cell size factor vector
    cell_size_factor <- SummarizedExperiment::colData(sce)[, "size_factors"]

    # create a new size factor of the same length with cell_size_factor but with value 1
    new_size_factor <- stats::setNames(rep(1, length(cell_size_factor)),
                                       names(cell_size_factor))

    # store the new_size_factor to sce
    SummarizedExperiment::colData(sce)[, "size_factors"] <- new_size_factor
  }

  # Initialize the pert object with the given pert and sce
  pert <- trt_group
  pert_object <- pert_input(pert, sce = sce, pert_level = "cre_perts")

  # Get all the guides that target the current `pert`
  pert_guides <- grna_target_data_frame |>
    dplyr::filter(grna_target == pert) |>
    dplyr::pull(grna_id)

  # Get perturbation status and gRNA perturbations for all cells
  pert_status <- SummarizedExperiment::colData(pert_object)$pert
  grna_perts <- SummarizedExperiment::assay(SingleCellExperiment::altExp(pert_object, "grna_perts"), "perts")

  # Convert to a sparse matrix, so the sampling function works in `create_guide_pert_status`
  grna_perts <- as(grna_perts, "CsparseMatrix")
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts, pert_guides = pert_guides)

  # Create effect size matrix (sampled from negative binomial distribution around effect_size or 1)
  effect_sizes <- structure(rep(effect_size, nrow(pert_object)), names = rownames(pert_object))

  # Create and center effect size matrices
  es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                      gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
  es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status, gene_effect_sizes = effect_sizes)
  es_mat_use <- es_mat[, colnames(SummarizedExperiment::assay(pert_object, "counts"))]

  # generate the gene count
  sim_counts <- sim_tapseq_sce(pert_object, effect_size_mat = es_mat_use)
  simulated_response_matrix <- as(SummarizedExperiment::assay(sim_counts, "counts"), "RsparseMatrix")

  # return the output
  return(list(response_matrix = simulated_response_matrix,
              grna_matrix = SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, "grna_perts"), "perts"),
              grna_target_data_frame = grna_target_data_frame,
              discovery_pairs =  discovery_pairs[discovery_pairs$grna_target == pert,]))
}
