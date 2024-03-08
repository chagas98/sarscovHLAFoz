#' Concatenate all prediction results from mutations directory
#'
#' @param working_directory Path to the root of the mutations directory
#' @return A strongly typed `arrow::Dataset` with the combined data
#' @import dplyr
#' @import tools
#' @import arrow
concat_mutations <- function(working_directory) {
  
  # Get a list of subdirectories in the working directory
  subdirectories <- list.dirs(working_directory, recursive = FALSE)
  
  # Initialize an empty data frame to store concatenated data
  combined_data <- data.frame()
  
  # Iterate over each subdirectory
  for (subdir in subdirectories) {
    # Get a list of CSV files in the current subdirectory
    tsv_files <- list.files(subdir, pattern = "\\_prediction_result.tsv$", full.names = TRUE)
    print(tsv_files)
    # Iterate over each CSV file and concatenate the data
    for (file in tsv_files) {
      data <- read.csv(file, sep='\t')
      data_clean <- data %>% 
        dplyr::select(-dplyr::contains(c("bam.PL", "bam.GT","AF1", 
                                        "BQBZ", "MQ", "MQ0F", 
                                        "MQBZ", "PV4", "RPBZ", 
                                        "SCBZ", "SGB", "VDB"))) %>% 
        dplyr::mutate(pos = as.integer(pos)) %>% 
        dplyr::mutate_if(is.character,
                list(~case_when(tolower(.)  == "true" ~ TRUE,
                                tolower(.)  == "false" ~ FALSE,
                                TRUE ~ .)))


      # Add the data from the current CSV file to the combined data
      combined_data <- dplyr::bind_rows(data_clean, combined_data)
    }
  }
  
  # Write the combined data to a new CSV file
  write.csv(combined_data, "results/combined_data.csv", row.names = FALSE)
  
  # Print a message indicating the process is complete
  cat("Data concatenation complete. Combined data saved to 'combined_data.csv'.\n")
  str(combined_data)
  type_dict  <-  list(
    "integer" = int64(),
    "string" = string(),
    "double" = float64(),
    "numeric" = float64(),
    "logical" = bool(),
    "Date" = date32(),
    "character" = string(),
    "NA" = NA
  )
  chosen_schema <- schema(
    purrr::map(names(combined_data), ~Field$create(name = .x, type = type_dict[[typeof(combined_data[[.x]])]]))
  )

  # Return the strongly typed combined data
  final_dataset  <- arrow::open_dataset(sources = "results/combined_data.csv",
                                        format = "csv",
                                        schema = chosen_schema, skip_rows = 1)

  return(final_dataset)
}

peptides_dataset <- concat_mutations("results/mutations")

query <- peptides_dataset  %>%  
 filter(length > 10)  %>% 
 collect()

query
