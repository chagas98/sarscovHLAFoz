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
library(gridExtra)
library(grid)
library(ggplot2)
library(stringr)

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
                             HLA.B.07.02.rank ~ "{median} ({p25}, {p75})",
                             HLA.B.07.02.rank_refseq ~ "{median} ({p25}, {p75})"
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

# concat
peptides_dataset <- concat_mutations(
  working_directory = "results/mutations")

# deduplicate  
peptides_deduplicate  <-  peptides_dataset %>%
  dplyr::filter(multiple_mutation == "single") %>%
  tidyr::drop_na(pos) %>%
  deduplicating(data = ., 
                columns = c("sequence", "length", "pos", "gene")) #variant_protein possui grupos de variantes iguais

nrow(peptides_deduplicate)
colnames(peptides_deduplicate)
peptides_deduplicate  %>% 
  count(gene, sort=TRUE)
###################################################################################
################################### TABLE 1 #######################################
###################################################################################

lista_labels = list(
  length = "k-mers",
  `HLA.B.07.02.rank` = "B07:02 rank Mutant",
  `HLA.B.07.02.rank_refseq` = "B07:02 rank WT")

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

###################################################################################
################################### FIG 1 #########################################
###################################################################################

impact_dataset <- peptides_deduplicate  %>%
  tidyr::pivot_longer(
    cols = contains('.binder'), 
    names_sep = '.binder', 
    names_to = c('alleles', 'seqclass')
  ) %>%
  dplyr::mutate(
    seqclass = case_when(
      grepl("_refseq", seqclass) ~ 'rank_refseq',
      TRUE ~ 'rank_mutant'
    )
  ) %>% 
    tidyr::pivot_wider(
    names_from = seqclass, 
    values_from = value
  ) %>%
  dplyr::mutate(
    alleles = gsub(".binder_refseq|.binder", "", alleles),
    alleles = gsub("HLA.B.", "B*", alleles),
    alleles = gsub("\\.", ":", alleles),
    impact = case_when(
      grepl("NB", rank_refseq) & grepl("SB", rank_mutant) ~ 'Gain',
      grepl("SB", rank_refseq) & grepl("NB", rank_mutant) ~ 'Loss',
      grepl("NB", rank_refseq) & grepl("WB", rank_mutant) ~ 'Weak Gain',
      grepl("SB", rank_refseq) & grepl("WB", rank_mutant) ~ 'Weak Loss',
      TRUE ~ 'No Effect')
  ) 

impact_dataset  %>% 
  count(gene)  %>% 
  print(n=200)

impact_dataset_gain  <- impact_dataset %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(grepl('Gain', impact))

impact_dataset_loss  <- impact_dataset %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(grepl('Loss', impact))

impact_dataset_noeffect  <- impact_dataset %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(grepl('No Effect', impact))

alleles_names  <- unique(impact_dataset$alleles)
alleles_sorted  <- str_sort(alleles_names, numeric = TRUE)

genes_names  <- unique(impact_dataset$gene)
genes_sorted  <- genes_names[order(nchar(genes_names))]


g.mid <- ggplot(impact_dataset,aes(x=1,y=factor(alleles, alleles_sorted)))+
  geom_text(aes(label=alleles))+
  geom_segment(aes(x=0.94,xend=0.96,yend=alleles))+
  geom_segment(aes(x=1.04,xend=1.065,yend=alleles))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))
g.mid
g1 <- ggplot(data = impact_dataset_gain, 
             aes(x = factor(alleles, alleles_sorted), 
             fill = impact)) +
  geom_histogram(stat = "count") + 
  theme(legend.position="left",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_x_discrete(drop=FALSE) +
  scale_y_reverse() + coord_flip()

g2 <- ggplot(data = impact_dataset_loss, 
             aes(x = factor(alleles, alleles_sorted),
             fill = impact)) +
  xlab(NULL) +
  geom_histogram(stat = "count") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  scale_x_discrete(drop=FALSE) +
  coord_flip()

gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

final_plot <- grid.arrange(g1,g.mid,g2, ncol=3,widths=c(4/9,1/9,4/9))
ggsave(file="test.pdf", final_plot,  width = 12, height = 8, dpi = 150) #saves g

###################################################################################
################################### FIG 2 #########################################
###################################################################################


g3 <- ggplot(data = impact_dataset_noeffect, aes(x = alleles)) +
  geom_histogram(stat = "count") + 
  theme(axis.title.x = element_blank(), 
        plot.margin = unit(c(1,-1,1,0), "mm"))
g3


###################################################################################
################################### FIG 3 #########################################
###################################################################################

g.mid2 <- ggplot(impact_dataset,aes(x=1,y=factor(gene, genes_sorted))) +
  geom_text(aes(label=gene))+
  geom_segment(aes(x=0.94,xend=0.96,yend=gene))+
  geom_segment(aes(x=1.04,xend=1.065,yend=gene))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))

g4 <- ggplot(data = impact_dataset_gain, 
             aes(x = factor(gene, genes_sorted), 
             fill = impact)) +
  geom_bar(stat = "count", width = 0.7) + 
  theme(legend.position="left",
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_reverse() + coord_flip()

g5 <- ggplot(data = impact_dataset_loss, aes(x = factor(gene, genes_sorted), fill = impact)) +xlab(NULL)+
  geom_bar(stat = "count", width = 0.7)+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        #axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  scale_x_discrete(drop = FALSE) +
  coord_flip()
g5
gg4 <- ggplot_gtable(ggplot_build(g4))
gg5 <- ggplot_gtable(ggplot_build(g5))
gg.mid2 <- ggplot_gtable(ggplot_build(g.mid2))

grid.arrange(g4,g.mid2,g5,ncol=3,widths=c(4/9,1/9,4/9))

