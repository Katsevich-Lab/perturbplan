# load path to raw data directory
devtools::load_all()
source("~/.Rprofile")
path_to_raw_data <- paste0(.get_config_path("LOCAL_TIAN_2019_RAW_DATA_DIR"), "/processed/iPSC")

# extract response matrix and read-UMI table
preprocessed_data <- reference_data_preprocessing_10x(path_to_raw_data)

# obtain the necessary information for baseline expression and library size information
preprocessed_results <- reference_data_preprocessing(response_matrix = preprocessed_data$response_matrix,
                                                     read_umi_table = preprocessed_data$read_umi_table)

# save the baseline information to data
usethis::use_data(preprocessed_results,
                  name = "iPSC_10x",  # object name users will type
                  overwrite = TRUE)          # allow reruns
