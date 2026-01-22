if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("readr" %in% rownames(installed.packages()) == FALSE) {install.packages("readr")}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyr")}
library (readr, quietly = T)
library (tidyr, quietly = T)
library(dplyr, quietly = T)

options(dplyr.summarise.inform = FALSE)

args = commandArgs(trailingOnly=TRUE)

#filter read data such that the end (3' or 5') of with the lowest level of read mapping is 1/10th or greater than the higher end of read mapping
sampleList = read_tsv(args[1], col_names = "Sample",
                      show_col_types = FALSE)


# -----------------------------
# READ step1 input safely
# -----------------------------
step1_path <- args[2]

if (!file.exists(step1_path) || file.info(step1_path)$size == 0) {

  message("INFO: ", step1_path, " is missing or empty. Writing empty output.")

  if (args[4] == "contig") {
    out <- tibble(
      sample = character(),
      IS = character(),
      contig = character(),
      Insertion_Position = double(),
      IS_Depth_at_Max_Depth_Site = double(),
      `5prime_IS_depth` = double(),
      `3prime_IS_depth` = double(),
      total_IS_depth = double(),
      non_IS_depth = double(),
      depth_percentage = double(),
      IS_fam = character(),
      ORF = character(),
      intragenic = character()
    )
    write_tsv(out, "final_results/pseudoR_output.contig.tsv")

  } else if (args[4] == "ORF") {
    out <- tibble(
      sample = character(),
      IS = character(),
      contig = character(),
      Insertion_Position = double(),
      IS_Depth_at_Max_Depth_Site = double(),
      `5prime_IS_depth` = double(),
      `3prime_IS_depth` = double(),
      total_IS_depth = double(),
      non_IS_depth = double(),
      depth_percentage = double(),
      IS_fam = character(),
      ORF = character()
    )
    write_tsv(out, "final_results/pseudoR_output.ORF.tsv")
  }

  quit(status = 0)
}

# -----------------------------
# Normal path: read step1 file
# -----------------------------
analysis <- read_tsv(
  step1_path,
  show_col_types = FALSE,
  col_types = cols(
    sample = col_character(),
    IS = col_character(),
    contig = col_character(),
    insertion_pos = col_double(),
    `5prime` = col_double(),
    `3prime` = col_double(),
    max_depth_site_insertion_pos = col_double(),
    max_depth_site_depth = col_double(),
    .default = col_guess()
  )
)

#continue original calculation
analysis = analysis %>%
  rowwise() %>%
  mutate (total_IS_depth = sum(`5prime`+`3prime`)) %>%
  filter (total_IS_depth >= 4) %>%
  mutate (low_end = min (c(`5prime`,`3prime`))) %>%
  mutate (high_end = max (c(`5prime`,`3prime`))) %>%
  filter (low_end >= 0.1*(high_end)) %>%
  distinct()

#process non-IS read mapping to generate non-IS read statistics
nonIS_mapping  = data.frame()
for (i in seq(nrow(sampleList))){
  temp = read_tsv (paste("mapping_files/", sampleList$Sample[i], ".contig.reads_mapped.depth", sep=""),
                   col_names = c("contig","start_pos","max_depth_site_insertion_pos","non_IS_depth"), 
                   show_col_types = FALSE) %>%
    mutate (sample = sampleList$Sample[i]) %>%
    mutate (max_depth_site_insertion_pos = ifelse(start_pos==0, 0, max_depth_site_insertion_pos))
  nonIS_mapping = rbind (nonIS_mapping, temp)
}
nonIS_mapping = nonIS_mapping %>% distinct()

analysis_annot = analysis %>% rowwise () %>%
  left_join (nonIS_mapping, by=c("contig" = "contig", "max_depth_site_insertion_pos" = "max_depth_site_insertion_pos",
                                 "sample"= "sample")) %>%
  mutate (depth_percentage = (100*(sum(`5prime`+`3prime`)/(sum(`5prime`+`3prime`)+non_IS_depth))))

analysis_annot = analysis_annot %>%
  select (sample, IS,contig, max_depth_site_insertion_pos, max_depth_site_depth,`5prime`, `3prime`,
          total_IS_depth, non_IS_depth, depth_percentage) %>%
  rename ("5prime_IS_depth" = `5prime`) %>%
  rename ("3prime_IS_depth" = `3prime`)

#this block removes IS elements false hits (likely) and adds a better naming strategy for the IS elements

IS_fam_df  = read_tsv(paste (args[3],"/IS_fam_annot.txt", sep="")) %>%
  rowwise () %>%
  mutate (IS = gsub("_", ":", IS))


analysis_annot = analysis_annot %>% 
  left_join(IS_fam_df)

#if (args[4] == "contig"){
#  inOrf  = read_tsv ("final_results/IS_hits_in_orfs.txt", 
#                     col_names = c("contig","start_pos","max_depth_site_insertion_pos", "ORF"),
#                     show_col_types = FALSE) %>%
#    select (contig,max_depth_site_insertion_pos,ORF) %>%
#    mutate (intragenic = "Yes") %>% distinct(contig, max_depth_site_insertion_pos, .keep_all = TRUE)
#  analysis_annot = analysis_annot %>%
#    left_join(inOrf, by=c("contig" = "contig", "max_depth_site_insertion_pos" = "max_depth_site_insertion_pos")) %>% 
#    mutate (intragenic =ifelse(is.na(intragenic), "No", intragenic)) %>%
#    distinct() %>%
#    rename ("Insertion_Position" = max_depth_site_insertion_pos) %>%
#    rename ("IS_Depth_at_Max_Depth_Site" = max_depth_site_depth) 
#  write_tsv (analysis_annot, "final_results/pseudoR_output.contig.tsv")
#} else if (args[4] == "ORF"){
#  analysis_annot = analysis_annot %>%
#   mutate ("ORF"=contig) %>%
#    rename ("Insertion_Position" = max_depth_site_insertion_pos) %>%
#    rename ("IS_Depth_at_Max_Depth_Site" = max_depth_site_depth)
#  write_tsv (analysis_annot, "final_results/pseudoR_output.ORF.tsv")
#}
if (args[4] == "contig"){

  inOrf_path <- "final_results/IS_hits_in_orfs.txt"

  if (!file.exists(inOrf_path) || file.info(inOrf_path)$size == 0) {
    inOrf <- tibble::tibble(
      contig = character(),
      max_depth_site_insertion_pos = numeric(),
      ORF = character(),
      intragenic = character()
    )
  } else {
    inOrf <- readr::read_tsv(
      inOrf_path,
      col_names = c("contig","start_pos","max_depth_site_insertion_pos","ORF"),
      col_types = readr::cols(
        contig = readr::col_character(),
        start_pos = readr::col_double(),
        max_depth_site_insertion_pos = readr::col_double(),
        ORF = readr::col_character()
      ),
      show_col_types = FALSE
    ) %>%
      select(contig, max_depth_site_insertion_pos, ORF) %>%
      mutate(intragenic = "Yes") %>%
      distinct(contig, max_depth_site_insertion_pos, .keep_all = TRUE)
  }

  # make sure join keys have the same type
  analysis_annot <- analysis_annot %>%
    mutate(max_depth_site_insertion_pos = as.double(max_depth_site_insertion_pos))

  analysis_annot = analysis_annot %>%
    left_join(inOrf, by=c("contig" = "contig", "max_depth_site_insertion_pos" = "max_depth_site_insertion_pos")) %>%
    mutate(intragenic = ifelse(is.na(intragenic), "No", intragenic)) %>%
    distinct() %>%
    rename("Insertion_Position" = max_depth_site_insertion_pos) %>%
    rename("IS_Depth_at_Max_Depth_Site" = max_depth_site_depth)

  write_tsv(analysis_annot, "final_results/pseudoR_output.contig.tsv")
}
