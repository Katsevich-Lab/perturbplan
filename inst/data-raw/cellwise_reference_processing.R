library(data.table)
library(rhdf5)

process_k562_gasperini <- function(path_to_dataset) {
  message("Start processing K562_Gasperini")
  path_to_run1 <- file.path(path_to_dataset, "processed", "at_scale" , "run1")
  path_to_run2 <- file.path(path_to_dataset, "processed", "at_scale", "run2")
  k562_data_1 <- perturbplan::reference_data_preprocessing_10x(path_to_run1)
  response_matrix_1 <- k562_data_1[[1]]
  read_umi_table <- k562_data_1[[2]]
  mapping_efficiency <- k562_data_1$mapping_efficiency
  k562_data_2 <- perturbplan::reference_data_preprocessing_10x(path_to_run2)
  response_matrix_2 <- k562_data_2[[1]]
  response_matrix <- cbind(response_matrix_1, response_matrix_2)
  
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency))
}

process_k562_10x <- function(path_to_dataset) {
  message("Start processing K562_10x")
  path_to_runs <- file.path(path_to_dataset, "processed")
  k562_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  response_matrix <- k562_data[[1]]
  read_umi_table <- k562_data[[2]]
  mapping_efficiency <- k562_data$mapping_efficiency
  
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency))
}

process_thp1_10x <- function(path_to_dataset) {
  message("Start processing THP1_Yao")
  path_to_runs <- file.path(path_to_dataset, "processed")
  dir_srrs <- list.dirs(path_to_runs, recursive = FALSE, full.names = TRUE)
  # Load the Seurat object from the specified path
  seurat_obj <- readRDS(file.path(path_to_runs, "GSM6858447_KO_conventional.rds"))
  
  # Select cells that have exactly one guide and it is either "non-targeting" or "safe-targeting"
  nt_cells <- colnames(seurat_obj[["perturbations"]])[which(
    grepl("non-targeting|safe-targeting", seurat_obj@meta.data$Guides) &
      seurat_obj@meta.data$Total_number_of_guides == 1
  )]
  
  # From these, select only the cells that are labeled with "-1" at the end (usually indicates 'first guide')
  nt_cells_once <- nt_cells[which(sub(".*-(.*)$", "\\1", nt_cells) == 1)]
  
  # Extract the gene expression matrix (raw UMI counts) from the "RNA" assay
  response_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  # Subset the expression matrix to include only the selected 'nt_cells_once' cells
  response_matrix <- response_matrix[, nt_cells]
  
  # check whether there are negative values
  if (any(response_matrix < 0)) {
    stop("Response matrix contains negative values, which is unexpected for UMI counts.")
  }
  
  # Read features.tsv.gz (usually contains 3 columns: gene_id, gene_name, gene_type)
  genes <- fread(file.path(dir_srrs[1], "outs" , "filtered_feature_bc_matrix" , "features.tsv.gz"), header = FALSE)
  
  # Rename columns for clarity (typical 10X Genomics format)
  colnames(genes) <- c("gene_id", "gene_name", "gene_type")
  
  # Create a mapping table from gene_name to gene_id, removing duplicates
  gene_map <- genes[!duplicated(gene_name), .(gene_name, gene_id)]
  
  # Keep only rows whose gene names are present in the mapping
  response_matrix <- response_matrix[rownames(response_matrix) %in% gene_map$gene_name, ]
  
  # Set rownames using data.table join (more robust than match)
  setkey(gene_map, gene_name)  # set key for fast join
  gene_ids <- gene_map[rownames(response_matrix), gene_id]
  
  # Replace rownames with gene IDs
  rownames(response_matrix) <- gene_ids
  
  srr_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  read_umi_table <- srr_data[[2]]
  mapping_efficiency <- srr_data$mapping_efficiency
  
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency))
}

process_t_cd8_10x <- function(path_to_dataset) {
  message("Start processing T_CD8_Shifrut")
  path_to_runs <- file.path(path_to_dataset, "processed")
  t_cd8_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  
  response_matrix <- t_cd8_data[[1]]
  read_umi_table <- t_cd8_data[[2]]
  mapping_efficiency <- t_cd8_data$mapping_efficiency
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency))
}


process_a549_10x <- function(path_to_dataset) {
  message("Start processing A549_10x")
  path_to_runs <- file.path(path_to_dataset, "processed")
  a549_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  response_matrix <- a549_data[[1]]
  read_umi_table <- a549_data[[2]]
  mapping_efficiency <- a549_data$mapping_efficiency
  
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency))
}

process_ipsc_10x <- function(path_to_dataset) {
  message("Start processing iPSC_10x")
  ipsc_data <- perturbplan::reference_data_preprocessing_10x(path_to_dataset)
  response_matrix <- ipsc_data[[1]]
  read_umi_table <- ipsc_data[[2]]
  mapping_efficiency <- ipsc_data$mapping_efficiency
  
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency,
                                                downsample_ratio=0.4))
}

process_ipsc_neuron_10x <- function(path_to_dataset) {
  message("Start processing iPSC_neuron_10x")
  ipsc_neuron_data <- perturbplan::reference_data_preprocessing_10x(path_to_dataset)
  response_matrix <- ipsc_neuron_data[[1]]
  read_umi_table <- ipsc_neuron_data[[2]]
  mapping_efficiency <- ipsc_neuron_data$mapping_efficiency
  
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency,
                                                downsample_ratio= c(0.1, 0.3, 0.5, 0.7)))
}


process_k562_tap <- function(path_to_dataset) {
  message("Start processing K562_TAP")
  path_to_runs <- file.path(path_to_dataset, "processed")
  k562_tap_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  response_matrix <- k562_tap_data[[1]]
  read_umi_table <- k562_tap_data[[2]]
  mapping_efficiency <- k562_tap_data$mapping_efficiency
  # read a gene list from file.path(path_to_dataset, "sample1", "outs", "target_panel.csv")
  # the gene list is in the first column without header (start from 7th row)
  # only keep unique genes
  panel_path <- file.path(path_to_runs, "perturbplan-demo", "outs", "target_panel.csv")
  if (!file.exists(panel_path)) {
    stop("target_panel.csv not found at: ", panel_path)
  }
  
  # helper: safe trim
  .trim <- function(x) gsub("^\\s+|\\s+$", "", x)
  
  gene_list <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    # first column, no header, start from row 7 -> skip first 6 lines
    dt <- data.table::fread(panel_path, header = FALSE, skip = 6, select = 1L, data.table = FALSE)
    gene_list <- unique(.trim(dt[[1]]))
  } else {
    # base R fallback
    con <- file(panel_path, open = "r"); on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    if (length(lines) < 7) stop("target_panel.csv has fewer than 7 lines.")
    first_col <- vapply(strsplit(lines[-(1:6)], ",", fixed = TRUE), function(z) z[[1]], character(1L))
    gene_list <- unique(.trim(first_col))
  }
  
  # drop empties / NAs just in case
  gene_list <- gene_list[!is.na(gene_list) & nzchar(gene_list)]
  
  if (length(gene_list) == 0L) {
    stop("Parsed gene_list is empty from ", panel_path)
  }
  
  return(perturbplan::reference_data_processing(response_matrix=response_matrix,
                                                read_umi_table=read_umi_table,
                                                mapping_efficiency = mapping_efficiency,
                                                gene_list = gene_list,
                                                downsample_ratio= c(0.1, 0.3, 0.5, 0.7)))
}
