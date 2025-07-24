# load path to raw data directory
devtools::load_all()
source("~/.Rprofile")
path_to_raw_data <- paste0(.get_config_path("LOCAL_TIAN_2019_RAW_DATA_DIR"), "/processed/iPSC-neuron")

# extract response matrix and read-UMI table
preprocessed_data <- reference_data_preprocessing_10x(path_to_raw_data)

# obtain the necessary information for baseline expression and library size information
preprocessed_results <- reference_data_preprocessing(response_matrix = preprocessed_data$response_matrix,
                                                     read_umi_table = preprocessed_data$read_umi_table,
                                                     D2_rough = 0.1,
                                                     downsample_ratio = c(0.1, 0.3, 0.5, 0.7))

# save the baseline information to data
usethis::use_data(preprocessed_results,
                  name = "iPSC_neuron_10x",  # object name users will type
                  overwrite = TRUE)          # allow reruns
