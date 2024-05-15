#!/usr/bin/env Rscript

###############################################################
###############################################################
# Author: Samuel Chagas de Assis
# Institution: Federal University of Latin American Integration (UNILA)
# Research Group: Medical Science Research Group
# Funding: National Council for Scientific and Technological Development (CNPq)
###############################################################
###############################################################

# Retrieve the variant sequences found in Foz do Iguazu and recorded in GISAID database.

#renv::restore() #Use this command to install renv and dependencies to run this code
print(getwd())

#Load R Packages
requiredPackages <- c("renv", "devtools", "readr", "GISAIDR")

for(p in requiredPackages){

  if(!require(p,character.only = TRUE)) {
    if(p == 'GISAIDR'){

      devtools::install_github("Wytamma/GISAIDR")

    } else {

      renv::install(p)
    }

    library(p,character.only = TRUE)
  }
}

# Get command-line arguments
command_line_args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were provided
if (length(command_line_args) == 0) {
  stop("No command-line arguments provided.")
}

# Access the first argument (assuming it's the variable you want)

# Extracting values from command-line arguments
get_arg_value <- function(arg_name) {
  arg_value <- grep(paste0('--', arg_name, '*'), command_line_args, value = TRUE)
  if (length(arg_value) > 0) {
    arg_value <- strsplit(arg_value, split = '=')[[1]][2]
  } else {
    arg_value <- NULL
  }
  return(arg_value)
}

# Extracting values for each parameter
renviron <- get_arg_value('renviron')
start_date <- get_arg_value('start_collection')
end_date <- get_arg_value('end_collection')
city_name <- get_arg_value('location')

print('ok')
print(renviron)
# Credentials
readRenviron(renviron)
username = Sys.getenv("GISAIDR_USERNAME")
password = Sys.getenv("GISAIDR_PASSWORD")

# GISAID Login 
credentials <- login(username = username, 
                     password = password, 
                     database="EpiCoV")

# Query
gisaid_ids <- query(
  credentials = credentials,
  location = city_name,
  from = start_date,
  to = end_date,
  fast = TRUE
)

print(gisaid_ids)
# Output
full_df <- download(credentials = credentials, 
                    list_of_accession_ids = gisaid_ids$accession_id)

write_csv(full_df, 'seq_info_gisaid.csv')

export_fasta(full_df, out_file_name = 'GISAID_sequences.fasta', 
             date_format='%Y-%m-%d', delimiter='|')



