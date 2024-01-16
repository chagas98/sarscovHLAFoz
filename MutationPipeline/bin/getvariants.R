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
#source('renv/activate.R')

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
renviron <- command_line_args[1]

# Credentials
readRenviron(renviron)
username = Sys.getenv("GISAIDR_USERNAME")
password = Sys.getenv("GISAIDR_PASSWORD")

print(username)
print(password)

# GISAID Login 
credentials <- login(username = username, 
                     password = password, 
                     database="EpiCoV")

# Query
gisaid_ids <- query(
  credentials = credentials,
  location = "Foz do Iguacu",
  from_subm = "2020-01-01",
  to_subm = "2023-01-01",
  fast = TRUE
)

# Output
full_df <- download(credentials = credentials, 
                    list_of_accession_ids = gisaid_ids$accession_id)

write_csv(full_df, 'seq_info_gisaid.csv')

export_fasta(full_df, out_file_name = 'GISAID_sequences.fasta', 
             date_format='%Y-%m-%d', delimiter='|')



