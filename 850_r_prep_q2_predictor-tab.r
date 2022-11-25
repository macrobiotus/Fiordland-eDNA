# *********************************************************
# * Combine and correct sample descriptors and predictors *
# ********************************************************* 
# 08-Dec-2020; 03-Feb-2021, 25-Nov-2022

# load packages
# =============

library("tidyverse")  # work using tibbles
library("readxl")     # read excel sheets
library("stringr")    # rename column names using dplyr

rm(list = ls(all.names = TRUE))
gc()


# load manifest file 
# ==================

# if only to learn sample id format
qmanif <- read_csv("/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/400_create_qiime_manifest__manifest_v2.txt") %>%
  arrange(desc(`sample-id`))
print(qmanif)

# load big_table to format as metadata file
# =========================================

load("/Users/paul/Documents/OU_eDNA/201028_Robjects/210127_200_r_metadata_management__big_table.Rdata")
print(big_table)

# generate a sample id table matching qiime sample id's in 
#     "/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/400_create_qiime_manifest__manifest.txt"
#   by checking 
#     "/Users/paul/Documents/OU_eDNA/200901_scripts/200_r_metadata_management.R"
#   so that values in column "qsamplid" match Qiime-internal sample names of file:
#     "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-tab.qza"

big_table_trimmed <- big_table %>%  
  mutate("qsamplid" = paste0(gsub("\\.", "-", key) , "-", sub('_', '-', sample_name))) %>%
  relocate("qsamplid") %>% arrange(desc("qsamplid"))

#  Create Qiime2-compatible metadata file
# =======================================

# Check if manifest and metadata sample is match up. 
#   Check if "big_table_trimmed$qsamplid" is the same as qmanif$`sample-id`, as it must become:
#   Some empty samples in "big_table" are not in Qiime manifest, but thats not a problem.

big_table_trimmed$qsamplid [!big_table_trimmed$qsamplid %in% qmanif$`sample-id`]

# Get metadata file: 
#   Widen "big_table_trimmed" so that it only has one unique value for "qsamplid"

metadata <- left_join(
                      big_table_trimmed %>% filter(primer_direction == "f_primer"),
                      big_table_trimmed %>% filter(primer_direction == "r_primer") %>%
                        select("key", "primer_direction", "primer_label", "adapter", "index", "templateprimer", "completeprimer") %>%
                        `colnames<-`(c("key", "primer_direction_rev", "primer_label_rev", "adapter_rev", "index_rev", "templateprimer_rev", "completeprimer_rev")),
                      by = c("key")
                      )
  
# Get metadata file: 
#   Trim to samples contained in manifest

metadata <- metadata[which(metadata$qsamplid %in% qmanif$`sample-id`), ]

# Check if manifest and metadata sample still match up - ok. 

metadata$qsamplid %in% qmanif$`sample-id`
qmanif$`sample-id` %in% metadata$qsamplid


# rename an sort columns for export ...
colnames(metadata) <- c("sampleid", "pool-content", "key", "primer-direction", "primer-label", "adapter-fwd",
  "flupad", "index-fwd", "template-primer-fwd", "complete-primer-fwd", "loc-name", "sample-time",
  "sample-type", "inside-reserve", "sample-name", "vol-l", "lat-dd", "long-dd", "depth-m", "notes",
  "xtr-date", "row", "plate", "col",  "ng-ul", "primer-direction_rev", "primer-label-rev", "adapter-rev", "index-rev", "template-primer-rev", "complete-primer-rev")

# .... more of the same and set library pooling amounts correctly for potential use of package "decontam()" downstream 
metadata <- metadata %>% relocate(c("sampleid", "key", "ng-ul", "lat-dd", "long-dd")) %>% replace_na(list("ng-ul" = 0 ))

metadata$`lat-dd` <- enc2utf8(metadata$`lat-dd`)
metadata$`long-dd` <- enc2utf8(metadata$`long-dd`)

#  Write Qiime2-compatible metadata file
# =======================================

p_compl <- "/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/850_prep_q2_predictor-tab__metadata.tsv"
write_tsv(metadata, p_compl, append = FALSE, col_names = TRUE)
