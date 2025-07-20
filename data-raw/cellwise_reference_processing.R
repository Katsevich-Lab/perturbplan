library(data.table)
library(rhdf5)
get_nt_only_barcodes <- function(run_dir) {
  h5_file <- file.path(run_dir, "molecule_info.h5")

  feature_idx  <- h5read(h5_file, "feature_idx")
  barcode_idx  <- h5read(h5_file, "barcode_idx")
  barcodes     <- h5read(h5_file, "barcodes")
  features     <- h5read(h5_file, "features/id")
  feature_type <- h5read(h5_file, "features/feature_type")

  # Identify rows corresponding to CRISPR guides
  is_guide <- feature_type == "CRISPR Guide"

  # Identify non-targeting guide names (e.g., containing "non" or "NTC")
  nontargeting_features <- grep("non|NTC", features, ignore.case = TRUE, value = TRUE)

  # Subset rows where features are guides
  guide_rows <- which(is_guide[feature_idx + 1])

  # Create a mapping of barcode to guide
  df <- data.frame(
    barcode = barcodes[barcode_idx[guide_rows] + 1],
    guide   = features[feature_idx[guide_rows] + 1],
    stringsAsFactors = FALSE
  )

  # Select barcodes that are only associated with non-targeting guides
  nt_only_barcodes <- df %>%
    group_by(barcode) %>%
    summarise(all_nt = all(guide %in% nontargeting_features), .groups = "drop") %>%
    filter(all_nt) %>%
    pull(barcode)

  return(nt_only_barcodes)
}

process_k562_10x <- function(path_to_dataset) {
  message("Start processing K562_10x")
  path_to_run1 <- file.path(path_to_dataset, "processed", "at_scale" , "run1")
  path_to_run2 <- file.path(path_to_dataset, "processed", "at_scale", "run2")
  k562_data_1 <- perturbplan::reference_data_preprocessing_10x(path_to_run1)
  response_matrix_1 <- k562_data_1[[1]]
  read_umi_table <- k562_data_1[[2]]
  k562_data_2 <- perturbplan::reference_data_preprocessing_10x(path_to_run2)
  response_matrix_2 <- k562_data_2[[1]]
  response_matrix <- cbind(response_matrix_1, response_matrix_2)

  return(perturbplan::reference_data_preprocessing(response_matrix=response_matrix,
                                                   read_umi_table=read_umi_table))
}

process_thp1_10x <- function(path_to_dataset) {
  message("Start processing THP1_10x")
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


  read_umi_table <- obtain_qc_read_umi_table(dir_srrs[1])
  read_umi_table <- read_umi_table|>
    dplyr::filter(cell_id %in% nt_cells_once)

  return(perturbplan::reference_data_preprocessing(response_matrix=response_matrix,
                                                   read_umi_table=read_umi_table))
}

process_t_cd4_10x <- function(path_to_dataset){
  message("Start processing T_CD4_10x")
  path_to_runs <- file.path(path_to_dataset, "processed")
  t_cd4_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  response_matrix <- t_cd4_data[[1]]
  read_umi_table <- t_cd4_data[[2]]
  nt_only_barcodes <- get_nt_only_barcodes(path_to_runs)
  response_matrix <- response_matrix[, colnames(response_matrix) %in% nt_only_barcodes]
  read_umi_table <- read_umi_table[read_umi_table$barcode %in% nt_only_barcodes, ]

  return(perturbplan::reference_data_preprocessing(response_matrix=response_matrix,
                                                   read_umi_table=read_umi_table))


}

process_t_cd8_10x <- function(path_to_dataset) {
  message("Start processing T_CD8_10x")
  path_to_runs <- file.path(path_to_dataset, "processed")
  t_cd8_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  response_matrix <- t_cd8_data[[1]]
  read_umi_table <- t_cd8_data[[2]]

  return(perturbplan::reference_data_preprocessing(response_matrix=response_matrix,
                                                   read_umi_table=read_umi_table))
}
