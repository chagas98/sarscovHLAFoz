################################################################################
## Robust inferences from the HLA project (Foz do Iguaçu - Brazil):           ##
##                                                                            ##
##                                                                            ##
##                                                                            ##
##                                                                            ##
## sc.assis.2017@aluno.unila.edu.br - Aug 2024                                ##
################################################################################

#rm(list = ls())
library(dplyr)
library(arrow)
library(tools)
library(purrr)
library(gtsummary)
library(tidyverse)

#sessioninfo::package_info()


#' Concatenate all prediction results from mutations directory
#'
#' @param working_directory Path to the root of the mutations directory
#' @return A strongly typed `arrow::Dataset` with the combined data
#' @import dplyr
#' @import tools
#' @import arrow
#' @export
#' @exportType arrow::Dataset
#' @example
#' concat_mutations("results/mutations")
#'
concat_mutations <- function(working_directory) {
  
  # Get a list of subdirectories in the working directory
  subdirectories <- list.dirs(working_directory, recursive = FALSE)
  
  # Initialize an empty data frame to store concatenated data
  combined_data <- data.frame()
  
  # Iterate over each subdirectory
  for (subdir in subdirectories) {
    # Get a list of CSV files in the current subdirectory
    tsv_files <- list.files(path = subdir, 
                            pattern = "\\_prediction_result.tsv$", 
                            full.names = TRUE)
    print(tsv_files)

    rename_columns <- c(
      variant_protein = "variant.details..protein.", 
      variant_genomic = "variant.details..genomic.")

    # Iterate over each CSV file and concatenate the data
    for (file in tsv_files) {
      data <- read.csv(file, sep='\t')

      columns_to_bool <- data %>%
        select_if(function(col) any(tolower(col) %in% c('false', 'true'))) %>%
        names()
      
      data_clean <- data %>% 
        dplyr::select(-dplyr::contains(c("bam.PL", "bam.GT","AF1", 
                                        "BQBZ", "MQ", "MQ0F", 
                                        "MQBZ", "PV4", "RPBZ", 
                                        "SCBZ", "SGB", "VDB"))) %>%
        dplyr::mutate(multiple_mutation = case_when(grepl(",", pos) ~ "double", 
                                                    TRUE ~ "single")) %>%
        dplyr::mutate(pos = as.integer(pos)) %>% 
        dplyr::mutate_at(vars(columns_to_bool), ~ ifelse(tolower(.) == 'true', TRUE, FALSE)) %>% 
        dplyr::rename(all_of(rename_columns))


      # Add the data from the current CSV file to the combined data
      combined_data <- dplyr::bind_rows(data_clean, combined_data)
    }
  }
  
  # Write the combined data to a new CSV file
  write.csv(combined_data, "results/combined_data.csv", row.names = FALSE)
  
  # Print a message indicating the process is complete
  cat("Data concatenation complete. Combined data saved to 'combined_data.csv'.\n")

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

  # Return  combined data
  #final_dataset  <- arrow::open_dataset(sources = "results/combined_data.csv",
  #                                      format = "csv",
  #                                      schema = chosen_schema, skip_rows = 1)

  return(combined_data)
}


#' Deduplicate rows
#' @param data
#' @param columns
#' @import dplyr
#' @export 
#' @exportType data.frame
#' @example
#' deduplicating(df , c("x", "y"))

deduplicating <- function(data, columns) {
  
  data %>%
    dplyr::distinct(pick(all_of(columns)), .keep_all = TRUE)
}


  ##Template - Design Summary Table
template_table <- function(data, div, include_list, translation=NA) {
  
  gtsummary::set_gtsummary_theme(
    theme_gtsummary_compact(set_theme = TRUE)) # https://cran.r-project.org/web/packages/gtsummary/vignettes/themes.html

  tmp <- data %>%
    gtsummary::tbl_summary(by = all_of(div), # stratification variable.
                           include = include_list,
                           percent = "column",
                           statistic = list(
                             all_continuous() ~ "{mean} ({sd})",
                             all_categorical() ~ "{n} ({p}%)",
                             HLA.B.07.02.affinity ~ "{median} ({p25}, {p75})",
                             HLA.B.07.02.affinity_refseq ~ "{median} ({p25}, {p75})"
                             #VentilatorySupport_days ~ "{median} ({p25}, {p75})",
                             #HospitalPeriod_days ~ "{median} ({p25}, {p75})",
                             #ICU_days ~ "{median} ({p25}, {p75})"
                           ),
                           digits = list(all_continuous() ~ c(1, 1),
                                         all_categorical() ~ 1),
                           label = translation) %>%
    bold_labels() %>%
    italicize_levels()

  return(tmp)

}

peptides_dataset <- concat_mutations(
  working_directory = "results/mutations")

  
peptides_deduplicate  <-  peptides_dataset %>%
  dplyr::filter(multiple_mutation == "single") %>%
  tidyr::drop_na(pos) %>%
  deduplicating(data = ., 
                columns = c("sequence", "pos", "gene")) #variant_protein possui grupos de variantes iguais

peptides_deduplicate  %>% 
  group_by(pos, gene,variant_protein)  %>% 
  tally(sort = T) %>%
  print(n = 200)


##Translation
lista_labels = list(
  length = "k-mers",
  `HLA.B.07.02.affinity` = "B07:02 score Mutant",
  `HLA.B.07.02.affinity_refseq` = "B07:02 score WT")
names(lista_labels)
  ##Collect between outcomes
tab_stats <- template_table(data = peptides_deduplicate,
                            div = "gene",
                            include_list = names(lista_labels),
                            translation = lista_labels)

summary_table_final <-
  gtsummary::tbl_merge(list(tab_stats)) %>%
  gtsummary::modify_header(label = "") 

summary_table_final
  #flextable::set_table_properties(width = 1, layout = "autofit")
  ##Collect between years
tab_statsYears <- template_table(data = data_final,
                                 div = "SamplesGroupYear",
                                 translation = lista_labels)

##Between outcomes in 2020
pval20 <- data_final %>%
  dplyr::filter(SamplesGroupYear %in% "2020") %>%
  template_table(., "Outcome_icu", lista_labels) %>%
  gtsummary::add_p(Age ~ "wilcox.test") %>%
  gtsummary::modify_header(p.value ~ "**Alta20 vs. Obito20**") %>%
  gtsummary::modify_column_hide(all_stat_cols())

  ##Between outcomes in 2021
pval21 <- data_final %>%
  dplyr::filter(SamplesGroupYear %in% "2021") %>%
  template_table(., "Outcome_icu", lista_labels) %>%
  gtsummary::add_p(Age ~ "wilcox.test") %>%
  gtsummary::modify_header(p.value ~ "**Alta21 vs. Obito21**") %>%
  gtsummary::modify_column_hide(all_stat_cols())

  ##Between incidence in years
pvalYears <- data_final %>%
  template_table(., "SamplesGroupYear", lista_labels) %>%
  gtsummary::add_p(Age ~ "wilcox.test") %>%
  gtsummary::modify_header(p.value ~ "**2020 vs. 2021**") %>%
  gtsummary::modify_column_hide(all_stat_cols())

  ##Gather all tables
summary_table_final <-
  gtsummary::tbl_merge(list(tab_stats, tab_statsYears,
                            pval20, pval21, pvalYears)) %>%
  gtsummary::modify_header(label = "") %>%
  gtsummary::modify_spanning_header(
    list(
      gtsummary::all_stat_cols() ~ "**Desfechos**",
      gtsummary::starts_with("p.value") ~ "**p-valores**"
    )
  ) %>%
  gtsummary::as_flex_table() %>%
  flextable::set_table_properties(width = 1, layout = "autofit")


peptides_dataset <- concat_mutations("results/mutations")

df <- tibble(
  x = sample(10, 100, rep = TRUE),
  y = sample(10, 100, rep = TRUE)
)

deduplicating(df, c('x', 'y'))
