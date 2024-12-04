# This is a Rscript computing the power function using score test

#' Power function for CRISPR screen experimental design
#'
#' @param control_cell_vec Vector including number of cell negative control group for each element
#' @param target_cell_mat Matrix including number of cell in group (l,k) (row: L genomic elements; column: K gRNA libraries)
#' @param UMI_s Estimate for UMI count per singlet for the underlying type of cell
#' @param library_size Library size from pilot data or set as a parameter in simulation study
#' @param expression_level Expression level for J genes (values between 0 and 1)
#' @param size_parameter_list Size parameter list including for all J genes
#' @param gRNA_efficiency_mat gRNA efficiency matrix of dimension L by K
#' @param recovery_rate Recovery rate for cells surviving the library preparation
#' @param doublet_rate Doublet rate for droplet containing more than one cell
#' @param doublet_factor atio of averaged UMI count per doublet to UMI per singlet
#' @param planned_read Planned total sequencing reads
#' @param mapping_efficiency Mapping efficiency for sequenced reads
#' @param effect_size Effect size matrix for each element l and gene j (L by J dimension)
#' @param sideness Left, right or both
#' @param correction Multiplicity correction, example including BH, bonferroni
#' @param sig_level False discovery rate level
#' @param QC A logic value; if TRUE, then do quality control
#' @param n_thresh QC threshold
#' @param intermediate_outcome A logic value indicating if only local mean and QC prob is output
#' @param size_factor_list List of size factor with two matrices
#'
#' @return Either power list or power list and discovery size
#' @importFrom dplyr if_else
#' @export

power_function <- function(control_cell_vec, target_cell_mat, UMI_s, library_size = NULL,  # cell-level parameter
                           expression_level, size_parameter_list,                      # gene expression parameters
                           gRNA_efficiency_mat,                                        # gRNA library parameters
                           recovery_rate, doublet_rate, doublet_factor,                # Library prep parameters
                           planned_read, mapping_efficiency,                           # Sequencing parameters
                           effect_size, sideness, correction = "BH", sig_level = 0.1,  # Test related parameters
                           QC = TRUE, n_thresh = 7, intermediate_outcome = FALSE,
                           size_factor_list = list(size_factor = matrix(1, nrow = nrow(target_cell_mat), ncol = 2,
                                                                        dimnames = list(element = NULL,
                                                                                        group = c("trt", "ctl"))),
                                                   size_factor_sq = matrix(1, nrow = nrow(target_cell_mat), ncol = 2,
                                                                           dimnames = list(element = NULL,
                                                                                           group = c("trt", "ctl"))))){

  # either library size or UMI_s has to be provided
  if(all(is.null(library_size), is.null(UMI_s))){
    stop("One of library size or UMI per singlet has to be provided!")
  }

  # compute the list of number of total planned cells
  target_cell_vec <- rowSums(target_cell_mat)
  planned_cell_list <- control_cell_vec + target_cell_vec

  # tell if the element in planned_cell_list is all the same or not
  if(length(unique(planned_cell_list)) > 1){
    stop("Sum of number of control cell and treatment cell should be the same!")
  }

  # compute the total number of cell
  total_cell <- planned_cell_list[1]
  num_element <- length(target_cell_vec)
  num_gene <- length(expression_level)

  # extract information
  size_factor <- size_factor_list$size_factor
  size_factor_sq <- size_factor_list$size_factor_sq

  # impute the size factor matrices
  size_factor_trt <- matrix(rep(size_factor[,"trt"], num_gene), ncol = num_gene)
  size_factor_ctl <- matrix(rep(size_factor[,"ctl"], num_gene), ncol = num_gene)
  size_factor_sq_trt <- matrix(rep(size_factor_sq[,"trt"], num_gene), ncol = num_gene)
  size_factor_sq_ctl <- matrix(rep(size_factor_sq[,"ctl"], num_gene), ncol = num_gene)

  # compute the gene expression related part in power formula
  if(is.null(library_size)){
    ## compute the reads per cell
    read_c <- read_per_cell(planned_read = planned_read,
                            mapping_efficiency = mapping_efficiency,
                            planned_cell = total_cell,
                            recovery_rate = recovery_rate)

    ## compute the averaged library size with read per cell
    library_size <- library_computation(UMI_s =  UMI_s,
                                        read_c = read_c,
                                        doublet_rate = doublet_rate,
                                        doublet_factor = doublet_factor)
  }

  # compute the gRNA efficiency
  gRNA_efficiency_list <- efficiency_computation(pi_mat = gRNA_efficiency_mat,
                                                 cell_mat = target_cell_mat)
  gRNA_mat <- matrix(rep(as.vector(gRNA_efficiency_list), times = num_gene), ncol = num_gene)

  # compute the gene expression related part in power formula
  gene_part <- gene_part_computation(expression_level_list = expression_level,
                                     library_size = library_size,
                                     effect_size = effect_size,
                                     num_control = control_cell_vec,
                                     num_trt = target_cell_vec,
                                     size_factor_trt = size_factor_trt,
                                     size_factor_ctl = size_factor_ctl)

  # compute the size parameter matrix
  size_mat <- matrix(rep(size_parameter_list, num_element),
                     ncol = num_gene, byrow = TRUE)

  # compute the different variance quantities
  pooled_mean <- gene_part$pooled_mat
  trt_mean <- gene_part$trt_mat
  ctl_mean <- gene_part$ctl_mat
  pooled_var <- var_nb(mean = pooled_mean, size = size_mat,
                       size_factor = matrix(1, nrow = num_element, ncol = num_gene),
                       size_factor_sq = matrix(1, nrow = num_element, ncol = num_gene))
  trt_var <- var_nb(mean = trt_mean, size = size_mat,
                    size_factor = size_factor_trt, size_factor_sq = size_factor_sq_trt)
  ctl_var <- var_nb(mean = ctl_mean, size = size_mat,
                    size_factor = size_factor_ctl, size_factor_sq = size_factor_sq_ctl)

  # compute the gene expression part depending if all gRNA is efficient or not
  if(min(gRNA_efficiency_list) == 1){

    # compute the denominator
    denominator <- sweep(pooled_var, 1,
                         (1/control_cell_vec + 1/target_cell_vec), "*")

    # compute the asymptotic mean
    asy_mean <- as.vector(ctl_mean * (exp(effect_size) - 1) / sqrt(denominator))

    # compute the asymptotic sd
    asy_var_1 <- sweep(ctl_var, 1, 1 / control_cell_vec, "*")
    asy_var_2 <- sweep(trt_var, 1, 1 / target_cell_vec, "*")
    asy_var <- (asy_var_1 + asy_var_2) / denominator
    asy_sd <- as.vector(sqrt(asy_var))

  }else{
    stop("Code is not finished when gRNA efficiency is not 1!")
  }

  # next step depends if QC_prob is used or not
  if(QC){

    # compute the QC_prob
    QC_prob <- QC_prob(effect_size = effect_size,
                       baseline_expression = expression_level * library_size,
                       size_parameter = size_parameter_list,
                       num_control = control_cell_vec, num_trt = target_cell_vec,
                       n_thresh = n_thresh)

    # compute the adjusted power and discovery set with QC probability computed
    adjusted_power_list <- adjusted_power(mean_list = asy_mean,
                                          sd_list = asy_sd,
                                          sig_level = sig_level,
                                          correction = correction,
                                          sideness = sideness,
                                          QC_prob = QC_prob)
  }else{

    # add QC probability
    QC_prob <- rep(0, length(asy_mean))

    # compute the adjusted power and discovery set with QC probability being 0
    adjusted_power_list <- adjusted_power(mean_list = asy_mean,
                                          sd_list = asy_sd,
                                          sig_level = sig_level,
                                          correction = correction,
                                          sideness = sideness,
                                          QC_prob = rep(0, length(asy_mean)))
  }

  # extract gRNA and gene names
  gRNA_name <- rownames(target_cell_mat)
  gene_name <- names(expression_level)

  # output depending the intermediate outcome is required or not
  if(intermediate_outcome){

    # create a dataframe with outered names
    result_df <- data.frame(
      pair = as.vector(outer(gRNA_name, gene_name, paste, sep = "_")),
      asy_mean = asy_mean,
      asy_sd = asy_sd,
      QC_prob = QC_prob
    )

    # return the dataframe
    return(result_df)

  }else{

    # create a dataframe with outered names
    result_df <- list(
      individual_power = data.frame(
        pair = as.vector(outer(gRNA_name, gene_name, paste, sep = "_")),
        adjusted_power = adjusted_power_list$adjusted_power
      ),
      num_discovery = adjusted_power_list$discovery_size_estimate
    )

    return(result_df)
  }
}
