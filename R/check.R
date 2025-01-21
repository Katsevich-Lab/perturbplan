#' Input checking function (currently only available for input_check)
#'
#' @inheritParams compute_power
#'
#' @return NULL
input_check <- function(
    ######################## type of perturb-seq experiment ####################
    perturb_type = NULL,

    ######################## power-determining parameters ######################
    target_cell_df = NULL, control_cell_vec = NULL,
    size_parameter = NULL, baseline_expression = NULL,
    effect_size_mean = NULL, effect_size_sd = NULL,

    ######################## QC-related parameters #############################
    n_nonzero_trt_thresh = NULL, n_nonzero_cntrl_thresh = NULL,

    ####################### test-related parameters ############################
    side = NULL, multiple_testing_method = NULL, multiple_testing_alpha = NULL, cutoff = NULL,

    ###################### output parameter ####################################
    intermediate_outcome = NULL
){
  ########################## perturb_type ######################################
  # argument missingness checking
  if(is.null(perturb_type)){
    stop("perturb_type should be specified!")
  }
  # now only CRISPRi can be considered
  if(perturb_type != "CRISPRi"){
    stop("Currently only CRISPRi experiments are supported.")
  }

  ########################## target_cell_df ####################################
  # argument missingness checking
  if(is.null(target_cell_df)){
    stop("target_cell_df should be specified!")
  }
  # data frame checking
  if(!is.data.frame(target_cell_df)){
    stop("target_cell_df should be a data frame!")
  }
  # colnames checking
  if(!dplyr::setequal(colnames(target_cell_df), c("gRNA_id", "gRNA_target", "num_cells"))){
    stop("Column names of target_cell_df should be exactly: gRNA_id, gRNA_target and num_cells!")
  }
  # data type checking
  if(!all(target_cell_df$num_cells %% 1 == 0)){
    stop("Column num_cells in target_cell_df should be integer!")
  }

  ########################## control_cell_vec ##################################
  # argument missingness checking
  if(is.null(control_cell_vec)){
    stop("control_cell_vec should be specified!")
  }
  # vector checking
  if(!is.vector(control_cell_vec)){
    stop("control_cell_vec should be a vector!")
  }
  # names missingness checking
  if(any(is.null(names(control_cell_vec)))){
    stop("Every position in control_cell_vec should be named with a particular genome element of interest!")
  }
  # unique names checking
  if(length(control_cell_vec) > length(unique(names(control_cell_vec)))){
    stop("Repeated control cell size for the same genome element is not allowed!")
  }
  ## interaction checking with target_cell_df
  # set of gRNA_target should be same for control and treatment cell
  if(!dplyr::setequal(names(control_cell_vec), unique(target_cell_df$gRNA_target))){
    stop("gRNA_target in target_cell_df should match the names of the contro_cell_vec if the same set of genome elements are of interest!")
  }

  ########################## size_parameter ####################################
  # argument missingness checking
  if(is.null(size_parameter)){
    stop("Size parameter for genes of interest should be specified!")
  }
  # vector checking
  if(!is.vector(size_parameter)){
    stop("size_parameter should be a vector!")
  }
  # names missingness checking
  if(any(is.null(names(size_parameter)))){
    stop("Every position in size_parameter should be named with a particular gene of interest!")
  }
  # unique names checking
  if(length(size_parameter) > length(unique(names(size_parameter)))){
    stop("Repeated size parameter for the same gene is not allowed!")
  }
  # positivity checking
  if(any(size_parameter <= 0)){
    stop("elements in size_parameter should all be positive!")
  }

  ########################## baseline_expression ###############################
  # argument missingness checking
  if(is.null(baseline_expression)){
    stop("baseline_expression should be specified!")
  }
  # vector checking
  if(!is.vector(baseline_expression)){
    stop("baseline_expression should be a vector!")
  }
  # names missingness checking
  if(any(is.null(names(baseline_expression)))){
    stop("Every position in baseline_expression should be named with a particular gene of interest!")
  }
  # unique names checking
  if(length(baseline_expression) > length(unique(names(baseline_expression)))){
    stop("Repeated baseline expression for the same gene is not allowed!")
  }
  # positivity checking
  if(any(baseline_expression <= 0)){
    stop("Elements in baseline_expression should all be positive!")
  }
  ## interaction checking with size_parameter
  # gene set should be same for size_parameter and baseline_expression
  if(!dplyr::setequal(names(baseline_expression), names(size_parameter))){
    stop("Set of genes provided in baseline_expression and size_parameter should be same!")
  }

  ########################## effect_size_mean ##################################
  # argument missingness checking
  if(is.null(effect_size_mean)){
    stop("Mean effect size matrix/scalar should be specified!")
  }
  # data type checking
  if(is.matrix(effect_size_mean)){
    # if matrix, rownames missingness checking
    if(any(is.null(rownames(effect_size_mean)))){
      stop("Rownames in matrix effect_size_mean should be specified!")
    }
    # colnames missingness checking
    if(any(is.null(colnames(effect_size_mean)))){
      stop("Colnames in matrix effect_size_mean should be specified!")
    }
    ## interaction checking with control_cell_vec
    if(!dplyr::setequal(rownames(effect_size_mean), names(control_cell_vec))){
      stop("Rownames of effect_size_mean should equal to the names of control_cell_vec containing all the genome elements of interest!")
    }
    ## interaction checking with baseline_expression
    if(!dplyr::setequal(colnames(effect_size_mean), names(baseline_expression))){
      stop("Colnames of effect_size_mean should equal to the names of baseline_expression!")
    }
  }else{
    # test if a scalar is inputted
    if(any(!is.numeric(effect_size_mean), (length(effect_size_mean) > 1))){
      stop("Either a scalar or a matrix should be inputted for effect_size_mean!")
    }
  }
  # value should be positive
  if(min(effect_size_mean < 0)){
    stop("The fold changes should be a positive number!")
  }

  ########################## effect_size_sd ####################################
  # argument missingness checking
  if(is.null(effect_size_sd)){
    stop("Sd effect size matrix/scalar should be specified!")
  }
  ## interaction checking with effect_size_mean
  if(!all(class(effect_size_sd) == class(effect_size_mean))){
    stop("Data type of effect_size_sd and effect_size_mean should be same!")
  }
  if(is.matrix(effect_size_sd)){
    # test if effect_size_sd and effect_size_mean have the same set of rownames and colnames
    if(dplyr::setequal(rownames(effect_size_sd), rownames(effect_size_mean))){
      if(!dplyr::setequal(colnames(effect_size_sd), colnames(effect_size_mean))){
        stop("Genes in effect_size_sd should be same as those in effect_size_mean!")
      }
    }else{
      stop("Genome elements in effect_size_sd should be same as those in effect_size_mean!")
    }
  }else{
    # test if both effect_size_sd and effect_size_mean have the same length
    if(length(effect_size_sd) != length(effect_size_mean)){
      stop("Length of effect_size_sd should be same as that of effect_size_mean!")
    }
  }

  #################### test-related parameters #################################
  # argument missingness checking
  if(is.null(cutoff)){
    if(any(is.null(c(multiple_testing_method, multiple_testing_alpha)))){
      stop("Either the cutoff should be specified, or the multiple testing method should be specified!")
    }else{
      # multple testing method should be BH or Bonferroni
      if(!(multiple_testing_method %in% c("BH", "Bonferroni"))){
        stop("Multple testing method should be either BH or Bonferroni!")
      }
    }
  }else{
    # cutoff should be between 0 and 1
    if(any(cutoff >= 1, cutoff <= 0)){
      stop("Input cutoff should be strictly between 0 and 1!")
    }
  }

  return(NULL)
}

