# load path to raw data directory
document::load_all()
source("~/.Rprofile")
path_to_raw_data <- paste0(.get_config_path("LOCAL_SAKELLAROPOULOS_2024_RAW_DATA_DIR"), "/processed")

# extract response matrix and read-UMI table
preprocessed_data <- reference_data_preprocessing_10x(path_to_raw_data)

# obtain the necessary information for baseline expression and library size information
preprocessed_results <- reference_data_preprocessing(response_matrix = response_matrix,
                                                     read_umi_table = read_umi_table)

# save the baseline information to data
saveRDS(preprocessed_results, "data/A549_10x.rds")

