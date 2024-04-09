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
library(patchwork)
library(ggplot2)
library(stringr)
library(RColorBrewer)

library(viridis)

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
  dplyr::filter(grepl('Gain', impact))

impact_dataset_loss  <- impact_dataset_deduplicate %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(grepl('Loss', impact))

impact_dataset_noeffect  <- impact_dataset_deduplicate %>% 
  #dplyr::filter(
  #  if_any(contains("impact"), ~ . != 'No Effect'))
  dplyr::filter(grepl('No Effect', impact))

alleles_names  <- unique(impact_dataset_deduplicate$alleles)
alleles_sorted  <- str_sort(alleles_names, numeric = TRUE)

genes_names  <- unique(impact_dataset_deduplicate$gene)
genes_sorted  <- genes_names[order(nchar(genes_names))]


g.mid <- ggplot(impact_dataset_deduplicate,aes(x=1,y=factor(alleles, alleles_sorted)))+
  geom_text(aes(label=alleles, fontface='plain', family='sans'),
            size=3.5)+
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
        plot.margin = unit(c(-5,-1, 1,-1), "mm"))

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

grid.arrange(gg1,gg.mid,gg2, ncol=3,widths=c(4/9,1/9,4/9))

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

g.mid2 <- ggplot(impact_dataset_deduplicate,aes(x=1,y=factor(gene, genes_sorted))) +
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
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  scale_x_discrete(drop = FALSE) +
  coord_flip()

gg4 <- ggplot_gtable(ggplot_build(g4))
gg5 <- ggplot_gtable(ggplot_build(g5))
gg.mid2 <- ggplot_gtable(ggplot_build(g.mid2))

grid.arrange(g4,g.mid2,g5,ncol=3,widths=c(2/5,1/5,2/5))

###################################################################################
################################### FIG 4 #########################################
###################################################################################

impact_dataset  %>% 
  count(variant_protein_sorted) %>% 
  print(n=300)

peptides_single  %>% 
  count(variant_protein_sorted) 


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

g.mid3 <- ggplot(impact_dataset_all,aes(x=1,y=factor(variant_lineage, lineages_sorted))) +
  geom_text(aes(label=variant_lineage))+
  geom_segment(aes(x=0.94,xend=0.96,yend=variant_lineage))+
  geom_segment(aes(x=1.04,xend=1.065,yend=variant_lineage))+
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

g6 <- ggplot(data = impact_dataset_gain_all, 
             aes(x = factor(variant_lineage, lineages_sorted), 
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

g7 <- ggplot(data = impact_dataset_loss_all, aes(x = factor(variant_lineage, lineages_sorted), fill = impact)) +xlab(NULL)+
  geom_bar(stat = "count", width = 0.7)+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  scale_x_discrete(drop = FALSE) +
  coord_flip()

gg6 <- ggplot_gtable(ggplot_build(g6))
gg7 <- ggplot_gtable(ggplot_build(g7))
gg.mid3 <- ggplot_gtable(ggplot_build(g.mid3))

grid.arrange(g6,g.mid3,g7,ncol=3,widths=c(4/9,1/9,4/9))

###################################################################################
################################### FIG 5 #########################################
###################################################################################
library(lubridate)
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
detection_sorted

lineage_date <- metadata_date %>%
  count(Lineage, month_year) %>%
  group_by(month_year) %>%         
  mutate(prop = prop.table(n)*100)

alleles_binned <- patients %>%
  filter(HospitalPeriod_days < 65)  %>% 
  mutate(month_year = format(as.Date(Hospital_admission), "%m-%Y")) %>%
  select(allele1, allele2, month_year) %>%
  pivot_longer(cols = c(allele1, allele2),
               names_to = 'type',
               values_to = 'alleles') %>%
  count(alleles) %>%
  group_by(alleles) %>%          
  mutate(alleles = ifelse(n < 7, "Others", alleles)) %>%
  ungroup() 
  
patients_date <- patients %>%
  filter(HospitalPeriod_days < 65) %>% 
  mutate(month_year = format(as.Date(Hospital_admission), "%m-%Y")) %>%
  select(allele1, allele2, month_year) %>%
  pivot_longer(cols = c(allele1, allele2),
               names_to = 'type',
               values_to = 'alleles') %>%
  mutate(alleles = ifelse(alleles %in% alleles_binned$alleles, alleles, "Others")) %>%
  count(alleles, month_year) %>%
  group_by(month_year) %>%          
  mutate(prop = prop.table(n))

library(viridis)
colourCount_variants = length(unique(lineage_date$Lineage))
getPalette_variants = colorRampPalette(brewer.pal(9, "Set1"))
magma_edit <- c("#000004FF", "#302e57", "#4f4c85", "#7612b9", "#9d12b9" , "#b464c4",
"#d6147f", "#e77e83", "#ffc7c7", "#94011a",  "#ce1f3c", "#e24c54", "#e03a07", "#ff8800", "#fbb861", "#FECE91FF", "#FDE6A8FF", "#c29e00", "#eeca2a", "#e9de40",
"#c9c9c9")

test <- ggplot(lineage_date, aes(x = lubridate::my(month_year), y = prop, fill = factor(Lineage, detection_sorted))) +
  geom_bar(stat = "identity",  colour="black") +
  #scale_fill_manual(values = c("A" = "blue", "B" = "red", "C" = "green")) +  # Change colors as needed
  labs( x = "",
       y = "Relative Prevalence (%)") +
  scale_fill_manual(values = getPalette_variants(colourCount_variants))+
  scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y", expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_blank(),
        legend.title = element_blank(),
        panel.grid=element_blank(),
        panel.background=element_blank())

test1 <- ggplot(lineage_date, aes(x = lubridate::my(month_year), y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y", expand = c(0.01,0)) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))


lineage_date$month_year
colourCount_alleles = length(unique(patients_date$alleles))
getPalette_alleles = colorRampPalette(brewer.pal(9, "BrBG"))


break.vec <- range(as.Date("01-03-2020", format="%d-%m-%Y"), as.Date("01-06-2022", format="%d-%m-%Y"))

test2 <- ggplot(patients_date, aes(x = lubridate::my(month_year), y = prop, fill=alleles)) +
  geom_bar(stat = "identity", colour='black') +
  labs( x = "",
        y = "Relative Prevalence (%)") +
  scale_x_date(limits = break.vec, date_breaks = '1 month', date_labels = "%b-%Y", expand = c(0.025,0)) +
  scale_fill_manual(values = getPalette_alleles(colourCount_alleles))+
  scale_y_continuous(expand = c(0.01,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_blank(),
        legend.title = element_blank(),
        panel.grid=element_blank(),
        panel.background=element_blank())


test3 <- ggplot(patients_date, aes(x = lubridate::my(month_year), y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  scale_x_date(limits = break.vec, date_breaks = '1 month', date_labels = "%b-%Y", expand = c(0.025,0)) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))

(test1 + test + test3 + test2 + plot_layout(ncol=1, heights = c(1, 4, 1, 4), guides = "collect"))

###################################################################################
################################### FIG 6 #########################################
###################################################################################

impact_dataset_all$fold_change <- log2(impact_dataset_all$rank_values_mutant / impact_dataset_all$rank_values_refseq)


FC_loss_all <- impact_dataset_all %>%
  dplyr::mutate(fold_change = log2(rank_values_mutant / rank_values_refseq)) %>% 
  dplyr::filter(impact %in% c('Loss', 'Weak Loss'))

library(ggnewscale)

cbp1 <- c("#ff0000", "#fd8282", "#ff7301", "#e2a807", "#f3d83f", "#206916", "#138f03", 
          "#0dad70", "#2464c4", "#4f90f1", "#00ccff", "#7c0494", "#e900ca", "#ee70ddc2")

library(grid)

muh_grob <- grid::rectGrob(
  x=0, y=1:21, gp=gpar(
    color='white', fill= getPalette_variants(colourCount_variants), alpha=0.8))

colourCount_alleles = length(unique(FC_loss_all$alleles))

ggplot(FC_loss_all)+
  geom_point(aes(x = fold_change, y = factor(variant_lineage, detection_sorted), colour = factor(length)), size = 5)+ 
  scale_colour_brewer("K-mers", type = "div", palette =  'Greys')+
  ggnewscale::new_scale_colour()+ 
  geom_point(aes(x = fold_change, y = factor(variant_lineage, detection_sorted), colour = alleles), size = 3, fill = NA)+ 
  scale_colour_manual('Alleles', values = getPalette_alleles(colourCount_alleles))+
  coord_cartesian(clip='off') +
  theme(
    plot.title = element_text(size=14, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=10, face="bold", colour = "black"), 
    axis.text.y = element_text(size=10, face="bold", colour = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey70"))+
  annotation_custom(
    grob=muh_grob, ymin = 0, ymax = 1, xmin = -0.55, xmax=0.25
  ) +
  scale_y_discrete(drop = FALSE)

