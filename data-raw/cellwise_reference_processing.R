process_k562_10x <- function(path_to_dataset) {
  path_to_runs <- file.path(path_to_dataset, "processed")
  k562_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  response_matrix <- k562_data[[1]]
  read_umi_table <- k562_data[[2]]
  h5 <- "molecule_info.h5"
  feature_idx <- h5read(h5, "feature_idx")
  barcode_idx <- h5read(h5, "barcode_idx")
  barcodes <- h5read(h5, "barcodes")
  features <- h5read(h5, "features/id")
  feature_type <- h5read(h5, "features/feature_type")
  
  # === Step 2: identify non-targeting features ===
  is_guide <- feature_type == "CRISPR Guide"
  nontargeting_features <- grep("non|NTC", features, ignore.case = TRUE, value = TRUE)
  
  # === Step 3: map molecules to barcodes and guides ===
  guide_rows <- which(is_guide[feature_idx + 1])
  df <- data.frame(
    barcode = barcodes[barcode_idx[guide_rows] + 1],
    guide   = features[feature_idx[guide_rows] + 1]
  )
  
  barcode_with_only_nt <- df %>%
    group_by(barcode) %>%
    summarise(all_nt = all(guide %in% nontargeting_features)) %>%
    filter(all_nt) %>%
    pull(barcode)
  
  # === Step 4: filter expression matrix ===
  response_matrix_nt <- response_matrix[, colnames(response_matrix) %in% barcode_with_only_nt]
  
  return(perturbplan::reference_data_preprocessing(response_matrix=response_matrix,
                                                   read_umi_table=read_umi_table))
}

process_thp1_10x <- function(path_to_dataset) {
  path_to_runs <- file.path(path_to_dataset, "processed")
  seurat_obj <- readRDS(file.path(path_to_runs, "GSM6858447_KO_conventional.rds"))
  nt_cells <- colnames(seurat_obj[["perturbations"]])[which(grepl("non-targeting|safe-targeting", seurat_obj@meta.data$Guides) & 
                                                              seurat_obj@meta.data$Total_number_of_guides==1)]
  nt_cells_once <- nt_cells[which(sub(".*-(.*)$", "\\1", nt_cells)==1)]
  
  thp1_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  response_matrix <- thp1_data[[1]]
  read_umi_table <- thp1_data[[2]]
  return(perturbplan::reference_data_preprocessing(response_matrix=response_matrix,
                                                   read_umi_table=read_umi_table))
}