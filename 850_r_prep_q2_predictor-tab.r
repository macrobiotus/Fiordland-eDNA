# *********************************************************
# * Combine and correct sample descriptors and predictors *
# ********************************************************* 
# 12-Aug-2020

# load packages
# =============
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse")  # work using tibbles
library("readxl")     # read excel sheets
library("stringr")    # rename column names using dplyr

# load sample descriptions
# =========================

# formerly
# p_mtd <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Manifest/200420_18S_MF_corrected.xlsx"
# recently - with DNA concentration values for possiby better cleanup
p_mtd <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Manifest/200727_18S_MF.xlsx"

tib_mtd <- read_excel(p_mtd) %>%
  mutate_at(vars(SampleID, BarcodeSequence, LinkerPrimerSequence, TmplHash,  LibContent, SardiID,  XtrOri, XtrContent, Location, Loci, Description ), factor) 
  
tib_mtd$LongDEC <- as.numeric(tib_mtd$LongDEC)
tib_mtd$LatDEC <- as.numeric(tib_mtd$LatDEC)

# check for duplicates - there should be none
duplicated(tib_mtd$SampleID)

# test 
# print.data.frame(tib_mtd)

# load predictor tables 
# =====================
# also set types and NA strings

p_atp <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Metadata/200419_pcm_atp.csv"
p_che <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Metadata/200419_pcm_csbp.csv"
p_geo <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Metadata/200419_pcm_geog.csv"
p_xrd <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Metadata/200812_pcm_xrd.csv"
p_age <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Metadata/200419_pcm_ages.csv"

tib_atp <- read_csv(p_atp,
  col_types = cols(.default = "d",  "SHA" = col_factor()),
  na = c("", "NA", "<NA>")
  )

tib_che <- read_csv(p_che,
  col_types = cols(.default = "d", "SHA" = col_factor(), "COLR" = col_factor(), "GRAVL" = col_factor()),
  na = c("", "NA", "<NA>")
  )

tib_geo <- read_csv(p_geo,
  col_types = cols(.default = "d", "SHA" = col_factor(), "ASPECT" = col_factor()),
  na = c("", "NA", "<NA>")
  )

tib_xrd <- read_csv(p_xrd,
  col_types = cols(.default = "d", "SHA" = col_factor()),
  na = c("", "NA", "<NA>")
  )

tib_age <- read_csv(p_age,
  col_types = cols(.default = "d",  "SHA" = col_factor()),
  na = c("", "NA", "<NA>")
  )

# merge data tables, using hash keys
# ==================================

# testing for missing data 
print.data.frame(anti_join(tib_mtd, tib_atp, by = c("XtrContent" = "SHA")))
print.data.frame(anti_join(tib_mtd, tib_che, by = c("XtrContent" = "SHA")))
print.data.frame(anti_join(tib_mtd, tib_geo, by = c("XtrContent" = "SHA")))
print.data.frame(anti_join(tib_mtd, tib_xrd, by = c("XtrContent" = "SHA")))
print.data.frame(anti_join(tib_mtd, tib_age, by = c("XtrContent" = "SHA")))

# Chemical data: merge tables and check for duplicates  
tib_compl <- tib_mtd %>% left_join(tib_che, by = c("XtrContent" = "SHA")) 
  
# Chemical data: duplicates removed:
tib_compl  %>% 
  group_by(XtrContent) %>% filter(n() != 1) %>%
  arrange(XtrContent) %>% print.data.frame()
  
# ATP data: merge tables and check for duplicates  
tib_compl <- tib_compl %>% left_join(tib_atp, by = c("XtrContent" = "SHA"))

# ATP data: check for duplicates - looks ok
tib_compl  %>% 
  group_by(XtrContent) %>% filter(n() != 1) %>%
  arrange(XtrContent) %>% print.data.frame()  

# XRD data: merge tables and check for duplicates
tib_compl <- tib_compl %>% left_join(tib_xrd, by = c("XtrContent" = "SHA"))

# XRD data: check for duplicates - looks ok
tib_compl %>%
  group_by(XtrContent) %>% filter(n() != 1) %>%
  arrange(XtrContent) %>% print.data.frame()
    
# Geographic data: merge tables and check for duplicates
tib_compl <- tib_compl %>% left_join(tib_geo, by = c("XtrContent" = "SHA"))

# Geographic data: check for duplicates - looks ok
tib_compl %>%
  group_by(XtrContent) %>% filter(n() != 1) %>%
  arrange(XtrContent) %>% print.data.frame()

# Age data: merge tables and check for duplicates
tib_compl <- tib_compl %>% left_join(tib_age, by = c("XtrContent" = "SHA"))

# Age data: check for duplicates - looks ok
tib_compl %>%
  group_by(XtrContent) %>% filter(n() != 1) %>%
  arrange(XtrContent) %>% print.data.frame()

# correct types
tib_compl <- tib_compl %>%  mutate_at(vars(XtrContent, SampleID), factor)

# save intermediate state after import
save(tib_compl, file = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/R/160_soil-data_raw.Rdata")
save.image(file = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/R/160_soil-data_raw_image.Rdata")

# correcting raw data  
# ====================

# check xrd values
tib_compl %>% 
  mutate(row_sum = rowSums(select(., names(tib_xrd)[2:10]))) %>% 
  dplyr:::select(names(tib_xrd)[2:10], row_sum) %>% print.data.frame()

# re-scale xrd values to yield a column sum of 1, i.e. close data
#  after rounding errors in Excel now just using the absolute value.

tib_compl <- tib_compl %>% 
  mutate_at(vars(names(tib_xrd)[2:10]), abs) %>% 
  mutate(row_sum = rowSums(select(., names(tib_xrd)[2:10]))) %>% 
  mutate_at(vars(names(tib_xrd)[2:10]), ~ ./row_sum) %>% 
  select(-row_sum) %>% arrange(SampleID)

# check xrd values again - looking good
tib_compl %>% 
  mutate(row_sum = rowSums(select(., names(tib_xrd)[2:10]))) %>% 
  select(names(tib_xrd)[2:9], row_sum) %>% print.data.frame()


tib_compl$CHLORITE
# plot all raw data
tib_compl %>% keep(is.numeric) %>%  gather() %>% 
  ggplot2:::ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_histogram() + geom_density() + theme_bw() +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(angle=0)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(angle = 0, hjust = 1,  size = 7), 
        axis.ticks.y = element_blank())

ggsave("200810_predictor_data_of_script_160.pdf", plot = last_plot(), 
       device = "pdf", path = "/Users/paul/Documents/OU_pcm_eukaryotes/Manuscript/200622_display_item_development/",
       scale = 3, width = 75, height = 50, units = c("mm"),
       dpi = 500, limitsize = TRUE)

# save intermediate state after corrections
save(tib_compl, file = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/R/160_soil-data_corrected.Rdata")
save.image(file = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/R/160_soil-data_corrected_image.Rdata")

# write_merged tables - for Qiime import
# ======================================

tib_compl
# formerly:
# p_compl <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Manifest/200501_18S_MF_merged.txt" 
# recently
p_compl <- "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/Manifest/200810_18S_MF_merged.txt"
write_tsv(tib_compl, p_compl, append = FALSE, col_names = TRUE)