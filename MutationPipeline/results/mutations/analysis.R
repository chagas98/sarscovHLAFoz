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
library(tools)
library(purrr)
library(arrow)
library(gtsummary)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(lubridate)
library(cowplot)
library(viridis)
library(ggnewscale)
library(grid)
library(ggrepel)

#sessioninfo::package_info()
setwd('../../')
getwd()
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
                             all_categorical() ~ "{n} ({p}%)"
                             #HLA.B.07.02.rank ~ "{median} ({p25}, {p75})",
                             #HLA.B.07.02.rank_refseq ~ "{median} ({p25}, {p75})"
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

# Function to translate mutation annotation
translate_mutation <- function(mutation_string) {
  # Mapping of three-letter amino acid codes to single-letter codes
  aa_mapping <- c("Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Cys" = "C", "Glu" = "E",
                  "Gln" = "Q", "Gly" = "G", "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K",
                  "Met" = "M", "Phe" = "F", "Pro" = "P", "Ser" = "S", "Thr" = "T", "Trp" = "W",
                  "Tyr" = "Y", "Val" = "V")

  # Split the string by comma if present
  parts <- unlist(strsplit(mutation_string, ","))

  # Get the first annotation
  annotation <- parts[1]

  # Extract information
  match <- regexpr("[[:digit:]]+", annotation)  # Find the position
  position <- substr(annotation, match, match + attr(match, "match.length") - 1)
  original_aa <- substr(annotation, 3, match - 1)
  mutated_aa <- substr(annotation, match + attr(match, "match.length"), nchar(annotation))

  # Translate three-letter amino acid codes to single-letter codes
  original_aa_single <- aa_mapping[original_aa]
  mutated_aa_single <- aa_mapping[mutated_aa]

  # Concatenate and return translated annotation
  translated_annotation <- paste0(original_aa_single, position, mutated_aa_single)
  return(translated_annotation)
}

# Function to sort elements in a string separated by commas
sort_string <- function(string) {
  sorted <- sort(unlist(strsplit(string, ",")))
  paste(sorted, collapse = ",")
}

# concat
peptides_dataset <- concat_mutations(
  working_directory = "results/mutations")

# deduplicate  
peptides_single  <-  peptides_dataset %>%
  dplyr::filter(multiple_mutation == "single") %>%
  tidyr::drop_na(pos)  %>% 
  dplyr::mutate(variant_protein_sorted = map(variant_protein, sort_string))  %>% 
  dplyr::mutate(variant_protein_sorted = map(variant_protein, translate_mutation))
  
peptides_deduplicate <- peptides_single %>% 
  deduplicating(data = ., 
                columns = c("sequence", 
                            "length", 
                            "pos", 
                            "gene")
  ) #variant_protein possui grupos de variantes iguais

###################################################################################
################################### METADATA ######################################
###################################################################################
metadata <- read.csv('results/mutations/metadata_info.csv', sep=',')
patients <- read.csv('results/raw_data/patients_hla.csv', sep=',')

metadata_date <- metadata %>%
  group_by(Date) %>% 
  mutate(month_year = format(as.Date(Date), "%m-%Y"),
         year = format(as.Date(Date), "%Y"),
         group_date = cur_group_id()
  ) %>% 
  ungroup()

detection_sorted  <- metadata_date %>%
  group_by(Lineage) %>%
  slice(which.min(group_date)) %>%
  select(Lineage, first_detection = group_date) %>% 
  arrange(first_detection) %>% 
  pull(Lineage) 


lineage_date <- metadata_date %>%
  count(Lineage, month_year) %>%
  group_by(month_year) %>%         
  mutate(prop = prop.table(n)*100)

###################################################################################
################################### MAIN DATA ######################################
###################################################################################

impact_dataset_all <- peptides_single  %>%
  dplyr::select(-dplyr::contains('.rank_changes')) %>%
  tidyr::pivot_longer(
    cols = contains('.rank'), 
    names_sep = '.rank', 
    names_to = c('alleles', 'rank_type')
  ) %>%
  dplyr::mutate(
    rank_type = case_when(
      grepl("_refseq", rank_type) ~ 'rank_values_refseq',
      TRUE ~ 'rank_values_mutant'
    )
  ) %>% 
    tidyr::pivot_wider(
    names_from = rank_type, 
    values_from = value
  ) %>%  
  dplyr::mutate(
    rank_mutant = case_when(
      rank_values_mutant < 0.5 ~ 'SB',
      rank_values_mutant >= 0.5 & rank_values_mutant < 2 ~ 'WB',
      rank_values_mutant >= 2 ~ 'NB',
      TRUE ~ 'NA'
    ),
    rank_refseq = case_when(
      rank_values_refseq < 0.5 ~ 'SB',
      rank_values_refseq >= 0.5 & rank_values_refseq < 2 ~ 'WB',
      rank_values_refseq >= 2 ~ 'NB',
      TRUE ~ 'NA'
    )
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

impact_dataset_deduplicate <- impact_dataset_all %>% 
  deduplicating(data = ., 
                columns = c("sequence", "length",
                            "pos", "gene", "alleles"))


impact_dataset_gain  <- impact_dataset_deduplicate %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(impact %in% c('Gain', 'Weak Gain'))

impact_dataset_loss  <- impact_dataset_deduplicate %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(impact %in% c('Loss', 'Weak Loss'))

impact_dataset_noeffect  <- impact_dataset_deduplicate %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(grepl('No Effect', impact))


impact_only_loss  <- impact_dataset_deduplicate %>% 
  dplyr::filter(impact == 'Loss')


impact_only_gain  <- impact_dataset_deduplicate %>% 
  dplyr::filter(impact == 'Gain')

impact_only_weakloss  <- impact_dataset_deduplicate %>% 
  dplyr::filter(impact == 'Weak Loss')


impact_only_weakgain  <- impact_dataset_deduplicate %>% 
  dplyr::filter(impact == 'Weak Gain')


###################################################################################
################################### TABLE 1 #######################################
###################################################################################

lista_labels = list(
  length = "k-mers",
  gene = 'Proteins')

##Collect between outcomes
tab_stats <- template_table(data = impact_dataset_deduplicate,
                            div = "impact",
                            include_list = names(lista_labels),
                            translation = lista_labels)

summary_table_final <-
  gtsummary::tbl_merge(list(tab_stats)) %>%
  gtsummary::modify_header(label = "") 

summary_table_final



###################################################################################
################################### FIG1C #########################################
###################################################################################

alleles_names  <- unique(impact_dataset_deduplicate$alleles)
alleles_sorted  <- str_sort(alleles_names, numeric = TRUE)

genes_names  <- unique(impact_dataset_deduplicate$gene)
genes_sorted  <- genes_names[order(nchar(genes_names))]

unique_alleles <- impact_dataset_deduplicate %>% 
  select(alleles) %>% 
  distinct()

g.mid <- ggplot(unique_alleles,aes(x=1,y=factor(alleles, alleles_sorted)))+
  geom_text(aes(label=alleles, fontface='plain', family='sans'),
            size=3.5)+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.99,1.01))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(-5,-1, 1,-1), "mm"))

g1 <- ggplot(data = impact_dataset_gain, 
             aes(x = factor(alleles, alleles_sorted), 
             fill = impact)) +
  geom_histogram(stat = "count") + 
  scale_fill_manual(values = c("gray30", "gray80")) +
  labs(x=NULL,
       fill = 'Efeito') +
  theme_classic() +
  theme(legend.position="left",
        legend.title = element_text(face="bold", size=13),
        legend.text =  element_text(face= "italic", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(size=13),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.bottom = element_line(color = "black",
                                          linewidth = 1.5,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "gray",
                                          linewidth = 0.5,
                                          linetype = 1),
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_x_discrete(drop=FALSE) +
  scale_y_reverse() + coord_flip()

g2 <- ggplot(data = impact_dataset_loss, 
             aes(x = factor(alleles, alleles_sorted),
             fill = impact)) +
  geom_histogram(stat = "count") +
  scale_fill_manual(values = c("red", "salmon")) +
  labs(x=NULL,
       fill = 'Efeito') +
  theme_classic() +
  theme(legend.title = element_text(face="bold", size=13),
        legend.text =  element_text(face= "italic", size=12),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x=element_text(size=13),
        axis.ticks.y = element_blank(),
        axis.line.y.left = element_blank(),
        axis.line.x.bottom = element_line(color = "black",
                                          linewidth = 1.5,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "gray",
                                          linewidth = 0.5,
                                          linetype = 1),
        plot.margin = unit(c(1,0,1,-1), "mm")
        ) +
  scale_x_discrete(drop=FALSE) +
  coord_flip()


fig1C <- (g1 + g.mid + g2)  + plot_layout(ncol=3, widths = c(4,1.4,4), guides = "collect")


###################################################################################
################################### FIG 1A ########################################
###################################################################################

unique_gene <- impact_dataset_deduplicate %>% 
  select(gene) %>% 
  distinct()

g.mid2 <- ggplot(unique_gene,aes(x=1,y=factor(gene, genes_sorted))) +
  geom_text(aes(label=gene, fontface='plain', family='sans'),
            size=3.5)+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.99,1.01))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(-5,-1, 1,-1), "mm"))

g4 <- ggplot(data = impact_dataset_gain, 
             aes(x = factor(gene, genes_sorted), 
             fill = impact)) +
  geom_histogram(stat = "count") + 
  scale_fill_manual(values = c("gray30", "gray80")) +
  labs(x=NULL,
       fill = 'Efeito') +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(size=13),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.bottom = element_line(color = "black",
                                          linewidth = 1.5,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "gray",
                                          linewidth = 0.5,
                                          linetype = 1),
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_x_discrete(drop=FALSE) +
  scale_y_reverse() + coord_flip()

g5 <- ggplot(data = impact_dataset_loss, aes(x = factor(gene, genes_sorted), fill = impact)) +xlab(NULL)+
  geom_histogram(stat = "count") +
  scale_fill_manual(values = c("red", "salmon")) +
  labs(x=NULL,
       fill = 'Efeito') +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x=element_text(size=13),
        axis.ticks.y = element_blank(),
        axis.line.y.left = element_blank(),
        axis.line.x.bottom = element_line(color = "black",
                                          linewidth = 1.5,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "gray",
                                          linewidth = 0.5,
                                          linetype = 1),
        plot.margin = unit(c(1,0,1,-1), "mm")
  ) +
  scale_x_discrete(drop=FALSE) +
  coord_flip()

fig1A <- g4  + g.mid2 + g5 + plot_layout(ncol=3, widths = c(4,1.4,4))


###################################################################################
################################### FIG 1B ########################################
###################################################################################


impact_dataset_gain_all <- dplyr::semi_join(
  impact_dataset_all, 
  impact_dataset_gain, 
  by = c("sequence", "alleles", "variant_protein")
)

impact_dataset_loss_all <- dplyr::semi_join(
  impact_dataset_all, 
  impact_dataset_loss, 
  by = c("sequence", "alleles", "variant_protein")
)

lineages_names  <- unique(peptides_single$variant_lineage)
lineages_sorted  <- str_sort(lineages_names)

unique_lineage <- impact_dataset_all %>% 
  select(variant_lineage) %>% 
  distinct()

g.mid3 <- ggplot(unique_lineage,aes(x=1,y=factor(variant_lineage, detection_sorted))) +
  geom_text(aes(label=variant_lineage, fontface='plain', family='sans'),
            size=3.5)+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.99,1.01))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(-5,-1, 1,-1), "mm"))


g6 <- ggplot(data = impact_dataset_gain_all, 
             aes(x = factor(variant_lineage, lineages_sorted), 
             fill = impact)) +
  geom_histogram(stat = "count") + 
  scale_fill_manual(values = c("gray30", "gray80")) +
  labs(x=NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(size=13),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.bottom = element_line(color = "black",
                                          linewidth = 1.5,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "gray",
                                          linewidth = 0.5,
                                          linetype = 1),
        plot.margin = unit(c(0,-0,0,0), "mm")) +
  scale_x_discrete(drop=FALSE) +
  scale_y_reverse() + coord_flip()

g7 <- ggplot(data = impact_dataset_loss_all, aes(x = factor(variant_lineage, lineages_sorted), fill = impact)) +xlab(NULL)+
  geom_histogram(stat = "count") +
  scale_fill_manual(values = c("red", "salmon")) +
  labs(x=NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x=element_text(size=13),
        axis.ticks.y = element_blank(),
        axis.line.y.left = element_blank(),
        axis.line.x.bottom = element_line(color = "black",
                                          linewidth = 1.5,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "gray",
                                          linewidth = 0.5,
                                          linetype = 1),
        plot.margin = unit(c(0,0,0,-0), "mm")
  ) +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  coord_flip()


fig1B <- g6  + g.mid3  + g7 + plot_layout(ncol=3, widths = c(4,1.4,4))



fig1 <- (wrap_elements(fig1A) | plot_spacer() | wrap_elements(fig1B) | plot_spacer() | wrap_elements(fig1C))  + 
  plot_layout(ncol=5, widths = c(1.3/5, 0.0005/5, 1.5/5, 0.0005/5, 1.7/5)) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(plot.tag = element_text(size = 26))

fig1
ggsave('Fig1.png', fig1, height = 14, width = 40, scale = 1,  units = "cm")


###################################################################################
################################### FIG 2 #########################################
###################################################################################

mutation_dataset <- impact_dataset_all %>% 
  deduplicating(data = ., columns = c("variant_protein_sorted", "gene", "variant_lineage")) %>% 
  group_by(pos,gene,variant_protein_sorted) %>% 
  count()


  
# Plot
fig2A <- ggplot(mutation_dataset, aes(x=pos, y=n)) +
  geom_segment(
    aes(x=pos, xend=pos, y=0, yend=n), 
    color=ifelse(mutation_dataset$n > 10, "orange", "grey"), 
    linewidth=ifelse(mutation_dataset$n > 10, 1, 0.7)
  ) +
  geom_point(
    color=ifelse(mutation_dataset$n > 10, "orange", "grey"), 
    size=ifelse(mutation_dataset$n > 10, 4, 2)
  ) +
  geom_point(data=filter(mutation_dataset, variant_protein_sorted %in% impact_only_loss$variant_protein_sorted),
    color="red", 
    size=2
  ) +
  geom_point(data=filter(mutation_dataset, variant_protein_sorted %in% impact_only_gain$variant_protein_sorted),
    color="gray30", 
    size=2
  ) +
#  geom_text_repel(data=filter(mutation_dataset, n>10), 
#                  aes(label=variant_protein_sorted),
#                  max.overlaps = Inf,
#                  min.segment.length = 10, seed = 42, point.padding = 3,
#                  max.time = 1, max.iter = 1e5,
#                  direction = "x",
#                 color='orange'
 # ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16), 
        axis.title.y = element_text(size=14), 
        axis.text.y = element_text(size=13), 
        axis.text.x=element_text(size=13),
        axis.ticks.y = element_blank(),) +
  scale_y_continuous(limits = c(0, 25), expand = c(0,0)) +
  scale_x_continuous(breaks=seq(0, 30000, 5000)) +
  ylab("Frequência de Mutações") +
  xlab("Posição (Nucleotídeos)")

fig2A
ggsave('Fig2A.png', fig2A, height = 8, width = 26, scale = 1,  units = "cm")

mutation_dataset_spike <-  mutation_dataset %>% 
  filter(gene == 'S') %>% 
  mutate(prot_pos = as.numeric(map(variant_protein_sorted, parse_number)))

paste0(mutation_dataset_spike$prot_pos, collapse = ' resid ')
fig2B <- ggplot(mutation_dataset_spike, aes(x=prot_pos, y=n)) +
  geom_segment(
    aes(x=prot_pos, xend=prot_pos, y=0, yend=n), 
    color=ifelse(mutation_dataset_spike$n > 10, "orange", "grey"), 
    size=ifelse(mutation_dataset_spike$n > 10, 1, 0.7)
  ) +
  geom_point(
    color=ifelse(mutation_dataset_spike$n > 10, "orange", "grey"), 
    size=ifelse(mutation_dataset_spike$n > 10, 4, 2)
  ) +
  geom_point(data=filter(mutation_dataset_spike, variant_protein_sorted %in% impact_only_loss$variant_protein_sorted),
             color="red", 
             size=2
  ) +
  geom_point(data=filter(mutation_dataset_spike, variant_protein_sorted %in% impact_only_gain$variant_protein_sorted),
             color="gray30", 
             size=2
  ) +
  geom_text_repel(data=filter(mutation_dataset_spike, n>10 && !(variant_protein_sorted %in% impact_only_loss$variant_protein_sorted)), 
                  aes(label=variant_protein_sorted),
                  max.overlaps = Inf,
                  min.segment.length = 10, seed = 42, point.padding = 3,
                  max.time = 1, max.iter = 1e5,
                  direction = "x",
                  color='orange'
  ) +
  geom_text_repel(data=filter(mutation_dataset_spike, variant_protein_sorted %in% impact_only_loss$variant_protein_sorted), 
                  aes(label=variant_protein_sorted),
                  color="red"
  ) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16), 
        axis.title.y = element_text(size=13), 
        axis.text.y = element_text(size=13), 
        axis.text.x=element_text(size=13),
        plot.margin = unit(c(7,6,5,6), "mm")
  ) +
  scale_y_continuous(limits = c(0, 25), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 1273), breaks=seq(0, 1200, 200), expand = c(0,0)) +
  ylab("Frequência de Mutações") +
  xlab("Resíduos")

fig2B
ggsave('Fig2B.png', fig2B, height = 7, width = 24, scale = 1,  units = "cm")

mutation_dataset_spike
###################################################################################
################################### FIG 5 #########################################
###################################################################################

alleles_binned <- patients %>%
  filter(HospitalPeriod_days < 65)  %>% 
  mutate(month_year = format(as.Date(Hospital_admission), "%m-%Y")) %>%
  select(allele1, allele2, month_year) %>%
  pivot_longer(cols = c(allele1, allele2),
               names_to = 'type',
               values_to = 'alleles') %>%
  count(alleles) %>%
  group_by(alleles) %>%          
  mutate(alleles = ifelse(n < 10, "Outros", alleles)) %>%
  ungroup() 
alleles_binned$alleles

# Select patients
patients_date <- patients %>%
  filter(HospitalPeriod_days < 65) %>% 
  mutate(month_year = format(as.Date(Hospital_admission), "%m-%Y")) %>%
  select(allele1, allele2, month_year) %>%
  pivot_longer(cols = c(allele1, allele2),
               names_to = 'type',
               values_to = 'alleles') %>%
  mutate(alleles = ifelse(alleles %in% alleles_binned$alleles, alleles, "Outros")) %>%
  count(alleles, month_year) %>%
  group_by(month_year) %>%          
  mutate(prop = prop.table(n))

# Get Colors
colourCount_variants = length(unique(lineage_date$Lineage))
getPalette_variants = colorRampPalette(brewer.pal(9, "Set1"))

plot_count_variants <-
  ggplot(lineage_date, aes(x = lubridate::my(month_year), y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  scale_x_date(
    date_breaks = '1 month',
    date_labels = "%b-%Y",
    expand = c(0.01, 0)
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "mm")
  )

plot_variants_month <-
  ggplot(lineage_date, aes(
    x = lubridate::my(month_year),
    y = prop,
    fill = factor(Lineage, detection_sorted)
  )) +
  geom_bar(stat = "identity",  colour = "black") +
  labs(x = "",
       y = "Prevalência Relativa (%)") +
  scale_fill_manual('Variantes \n SARS-CoV-2', 
                    values = getPalette_variants(colourCount_variants)
  ) +
  scale_x_date(
    date_breaks = '1 month',
    date_labels = "%b-%Y",
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.x = element_blank(),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 10),
    legend.title.align = 0.3,
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

# get colors
colourCount_alleles = length(unique(patients_date$alleles))
getPalette_alleles = colorRampPalette(brewer.pal(9, "BrBG"))
# extend x axis based on variants
break.vec <- range(as.Date("01-03-2020", format="%d-%m-%Y"), as.Date("01-06-2022", format="%d-%m-%Y"))


plot_count_alleles <-
  ggplot(patients_date, aes(x = lubridate::my(month_year), y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  scale_x_date(
    limits = break.vec,
    date_breaks = '1 month',
    date_labels = "%b-%Y",
    expand = c(0.025, 0)
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "mm")
  ) 

plot_alleles_month <- ggplot(patients_date, aes(x = lubridate::my(month_year), y = prop, fill=alleles)) +
  geom_bar(stat = "identity", colour='black') +
  labs( x = "",
        y = "Prevalência Relativa (%)") +
  scale_x_date(limits = break.vec, date_breaks = '1 month', date_labels = "%b-%Y", expand = c(0.025,0)) +
  scale_fill_manual('Alelos', values = getPalette_alleles(colourCount_alleles)) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_blank(),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 10),
        legend.title.align = 0.2,
        panel.grid=element_blank(),
        panel.background=element_blank())

fig5 <- plot_count_variants + 
  plot_variants_month + 
  plot_count_alleles +
  plot_alleles_month + 
  plot_layout(ncol=1, heights = c(1.2, 4, 1.2, 4), guides = "auto") +
  plot_annotation(tag_levels = list(c("A", "", "B", "")))

fig5

ggsave('Fig5.png', fig5, height = 22, width = 27, scale = 0.8,  units = "cm")
lineage_date

###################################################################################
################################### FIG 6 #########################################
###################################################################################
#"B*07:02" "B*08:01" "B*14:02" "B*15:01" "B*18:01" "B*27:05" "B*35:01" "B*38:01" "B*39:01"
#"B*40:01" "B*44:03" "B*49:01" "B*51:01" "B*53:01" "B*57:01"

impact_dataset_all$fold_change <- log2(impact_dataset_all$rank_values_mutant / impact_dataset_all$rank_values_refseq)


FC_loss_all <- impact_dataset_all %>%
  dplyr::mutate(fold_change = log2(rank_values_mutant / rank_values_refseq)) %>% 
  dplyr::filter(impact %in% c('Loss', 'Weak Loss'))

muh_grob <- grid::rectGrob(
  x=0, y=1:21, gp=gpar(
    color='white', fill= getPalette_variants(colourCount_variants), alpha=0.8))

colourCount_alleles = length(unique(FC_loss_all$alleles))

fig6 <- ggplot(FC_loss_all)+
  ggnewscale::new_scale_colour()+ 
  geom_point(aes(x = fold_change, y = factor(variant_lineage, detection_sorted), 
                 colour = alleles, shape=factor(length)), 
            size = 3, 
            fill = NA, 
            stroke=1)+ 
  labs(x=expression(phantom()*Log[2]*FC), 
       y = '',
       shape="K-mers") +
  coord_cartesian(clip='off') +
  theme(
    plot.title = element_text(size=14, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=15, face="bold", colour = "black"), 
    axis.text.y = element_text(size=10, face="bold", colour = "white"),
    legend.title = element_text(face = "bold", size=13),
    legend.text =  element_text(size=12),
    legend.title.align = 0.4,
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey70")
  ) +
  scale_colour_manual('Alelos', values = getPalette_alleles(colourCount_alleles)) +
  scale_shape_manual(values = c(15, 18, 16, 17)) +
  annotation_custom(
    grob=muh_grob, ymin = 0, ymax = 1, xmin = -0.55, xmax=0.25
  ) +
  scale_y_discrete(drop = FALSE)

fig6
ggsave('Fig6.png', fig6, height = 15, width = 27, scale = 1,  units = "cm")

###################################################################################
################################### FIG 7 ########################################
###################################################################################
library(ggpubr)
library(rstatix)
impact_dataset_comparative <-  impact_dataset_deduplicate %>% 
  dplyr::filter(str_detect(impact, 'Loss|Gain')) %>% 
  mutate(paired = row_number(),
         impact_geral = case_when(
           str_detect(impact, 'Loss') ~ 'Total Loss',
           str_detect(impact, 'Gain') ~ 'Total Gain')
  ) %>% 
  tidyr::pivot_longer(
    cols = contains('rank_values_'),
    names_prefix = 'rank_values_',
    names_to = 'rank_type') %>% 
  mutate(log2value = log2(value),
         rank_type = case_when(
           rank_type=='mutant' ~ 'Mutação',
           rank_type=='refseq' ~ 'Referência'
         ))

fig7 <- ggplot(impact_dataset_comparative, aes(x=rank_type, y=log2value)) + 
  geom_boxplot(width = 0.3, size = 0.4
  ) +
  geom_signif(comparisons = list(c('Mutação', 'Referência'))
  ) +
  geom_dotplot(
    aes(fill = impact_geral),
    alpha = 0.7,
    dotsize = 2,
    binaxis='y', stackdir='center'
  ) +
  geom_line(aes(group = paired), 
            color = "black",
            alpha = 0.5
  ) +
  scale_fill_manual(values = c("gray30", "red")) +
  labs(y = expression(phantom()*Log[2]*('Afinidade de Ligação')),
       fill = 'Efeito') +
  theme_bw()+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=13),
        legend.text =  element_text(face= "italic", size=12),
        axis.text.x = element_text(face="bold", size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(face="bold", size=14),
        strip.text.x = element_text(size = 13)
  ) +
  ylim(NA, 10) +
  facet_wrap(~alleles)

ggsave('Fig7.png', fig7, height = 17, width = 25, scale = 1,  units = "cm")


###################################################################################
################################### FIG 8 ########################################
###################################################################################


###################################################################################
################################### TABLES ########################################
###################################################################################

# Fig5
patients_date
lineage_date

# Fig6
FC_loss_all

# Fig7
impact_dataset_comparative
