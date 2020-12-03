# Fiordland project 03-12-2020
#
# Read cells from 
#   "/Users/paul/Documents/OU_eDNA/200128_lab_work/200907_plate_layouts.xlsx"
#   "/Users/paul/Documents/OU_eDNA/191031_primers/200302_Bunce_et_al_0000_MiFishEmod_single_step_primers.xlsx"  
# output appropriate format for
#   "/Users/paul/Documents/OU_eDNA/200128_lab_work/201028_Otago_Genomics_sample_info.xlsx"
#    and later possibly the mapping file - with format as described here
#      `https://docs.qiime2.org/2020.8/tutorials/metadata/`
# example file at `/Users/paul/Documents/OU_eDNA/201126_script_scratch/sample-metadata_example.tsv`

# Environment
# ===========
library(tidyverse)

rm(list = ls())



# Read in sample data
# ===================

path <- c(
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200907_plate_layouts.xlsx",
  "/Users/paul/Documents/OU_eDNA/191031_primers/200302_Bunce_et_al_0000_MiFishEmod_single_step_primers.xlsx"
  )

# later read in barcodes as well?
 
# pool plate U 1 samples and primers
# ----------------------------------

u1s <- readxl::read_excel(path[1], range = "A3:M11", .name_repair = "universal")
u1pf <- readxl::read_excel(path[1], range = "A46:M54", .name_repair = "universal")
u1pr <- readxl::read_excel(path[1], range = "A56:M64", .name_repair = "universal") 

u1s <- u1s %>% pivot_longer(!pool.plate.U.1) %>% rename(row = pool.plate.U.1, col = name, content = value) %>% add_column(plate = "U.1", .before = TRUE) %>% add_column(type = "amplicon", .before = TRUE)
u1pf <- u1pf %>% pivot_longer(!U.1.fwd) %>% rename(row = U.1.fwd, col = name, content = value) %>% add_column(plate = "U.1", .before = TRUE) %>% add_column(type = "f_primer", .before = TRUE)
u1pr <- u1pr %>% pivot_longer(!U.1.rev) %>% rename(row = U.1.rev, col = name, content = value) %>% add_column(plate = "U.1", .before = TRUE) %>% add_column(type = "r_primer", .before = TRUE)

# pool plate U 2 samples and primers
# ---------------------------------

u2s <- readxl::read_excel(path[1], range = "A13:M21", .name_repair = "universal") 
u2pf <- readxl::read_excel(path[1], range = "A67:M75", .name_repair = "universal")
u2pr <- readxl::read_excel(path[1], range = "A77:M85", .name_repair = "universal")

u2s <- u2s %>% pivot_longer(!pool.plate.U.2) %>% rename(row = pool.plate.U.2, col = name, content = value) %>% add_column(plate = "U.2", .before = TRUE) %>% add_column(type = "amplicon", .before = TRUE)
u2pf <- u2pf %>% pivot_longer(!U.2.fwd) %>% rename(row = U.2.fwd, col = name, content = value) %>% add_column(plate = "U.2", .before = TRUE) %>% add_column(type = "f_primer", .before = TRUE)
u2pr <- u2pr %>% pivot_longer(!U.2.rev) %>% rename(row = U.2.rev, col = name, content = value) %>% add_column(plate = "U.2", .before = TRUE) %>% add_column(type = "r_primer", .before = TRUE)

# pool plate E 1 samples and primers
# ---------------------------------

e1s <- readxl::read_excel(path[1], range = "A24:M32", .name_repair = "universal")
e1pf <- readxl::read_excel(path[1], range = "A88:M96", .name_repair = "universal")
e1pr <- readxl::read_excel(path[1], range = "A98:M106", .name_repair = "universal")

e1s <- e1s %>% pivot_longer(!pool.plate.E.1) %>% rename(row = pool.plate.E.1, col = name, content = value) %>% add_column(plate = "E.1", .before = TRUE) %>% add_column(type = "amplicon", .before = TRUE)
e1pf <- e1pf %>% pivot_longer(!E.1.fwd) %>% rename(row = E.1.fwd, col = name, content = value) %>% add_column(plate = "E.1", .before = TRUE) %>% add_column(type = "f_primer", .before = TRUE)
e1pr <- e1pr %>% pivot_longer(!E.1.rev) %>% rename(row = E.1.rev, col = name, content = value) %>% add_column(plate = "E.1", .before = TRUE) %>% add_column(type = "r_primer", .before = TRUE)

# pool plate E 2 samples and primers
# ----------------------------------

e2s <- readxl::read_excel(path[1], range = "A34:M42", .name_repair = "universal")
e2pf <- readxl::read_excel(path[1], range = "A109:M117", .name_repair = "universal")
e2pr <- readxl::read_excel(path[1], range = "A119:M127", .name_repair = "universal")

e2s <- e2s %>% pivot_longer(!pool.plate.E.2 ) %>% rename(row = pool.plate.E.2, col = name, content = value) %>% add_column(plate = "E.2", .before = TRUE) %>% add_column(type = "amplicon", .before = TRUE) 
e2pf <- e2pf %>% pivot_longer(!E.2.fwd) %>% rename(row = E.2.fwd, col = name, content = value) %>% add_column(plate = "E.2", .before = TRUE) %>% add_column(type = "f_primer", .before = TRUE)
e2pr <- e2pr %>% pivot_longer(!E.2.rev) %>% rename(row = E.2.rev, col = name, content = value) %>% add_column(plate = "E.2", .before = TRUE) %>% add_column(type = "r_primer", .before = TRUE)

# get large data frames for further manipulation
# ---------------------------------------------

amplicon_positions <- bind_rows(u1s, u2s, e1s, e2s) %>% 
  mutate(col = str_remove(col, "...")) %>% 
  mutate(col, col = as.integer(col)) %>% 
  mutate(key = paste0(plate,".",row,".",col)) %>% print(n = Inf)
  
# add and format variable for merging with metadata object `mdata`
#  remove everything after first whitespace
amplicon_positions$sample_name <- sub(" .*", "", amplicon_positions$content)

# remove everything that isn't a sample id
amplicon_positions <- amplicon_positions %>% 
  mutate(sample_name = case_when(grepl("PC", sample_name) ~ sample_name, grepl("AK", sample_name) ~ sample_name)) %>%
  print(n = Inf)
 
save(amplicon_positions, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__amplicon_positions.Rdata")

primer_positions <- bind_rows(u1pf, u1pr, u2pf, u2pr, e1pf, e1pr, e2pf, e2pr) %>% 
  mutate(col = str_remove(col, "...")) %>% 
  mutate(col, col = as.integer(col)) %>% 
  mutate(key = paste0(plate,".",row,".",col))
save(primer_positions, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__primer_positions.Rdata")


# Read in primer data
# ===================

up5 <- readxl::read_excel(path[2], range = "E4:J11", .name_repair = "universal",  col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))
up7 <- readxl::read_excel(path[2], range = "E16:J32", .name_repair = "universal", col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))

ep5 <- readxl::read_excel(path[2], range = "E40:J47", .name_repair = "universal", col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))
ep7 <- readxl::read_excel(path[2], range = "E52:J68", .name_repair = "universal", col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))

# get large data frame for further manipulation
# ---------------------------------------------

primer_sequences <- bind_rows(up5, up7, ep5, ep7)
save(primer_sequences, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__primers.Rdata")


# Read amplicon concentration data prior to pooling 
# =================================================

# read in 

qpath <- c(
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate1_HS_ng-ul.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate2_HS_ng-ul.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate1_HS_ng-ul.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate2_HS_ng-ul.xlsx"
)

e1concs <- readxl::read_excel(qpath[1], range = "A4:M11", .name_repair = "universal", col_names =  c("row", as.character(seq(1,12))), col_types = c("guess", rep("numeric", 12))) %>%
  add_column(plate = "E.1")
e2concs <- readxl::read_excel(qpath[2], range = "A4:M11", .name_repair = "universal", col_names =  c("row", as.character(seq(1,12))), col_types = c("guess", rep("numeric", 12))) %>%
  add_column(plate = "E.2")
u1concs <- readxl::read_excel(qpath[3], range = "A4:M11", .name_repair = "universal", col_names =  c("row", as.character(seq(1,12))), col_types = c("guess", rep("numeric", 12))) %>%
  add_column(plate = "U.1")
u2concs <- readxl::read_excel(qpath[4], range = "A4:M11", .name_repair = "universal", col_names =  c("row", as.character(seq(1,12))), col_types = c("guess", rep("numeric", 12))) %>%
  add_column(plate = "U.2")
  
# reformat for merging with outher tables (such as `amplicon_positions and `primer_positions`)
#   col names:   plate row col key     
#   var format:  U.1   A     1 U.1.A.1 

list_concs_long <- lapply(list(e1concs, e2concs, u1concs, u2concs), function(x) x %>% pivot_longer(cols = !c(row, plate), names_to = "col", values_to = "ng_ul") %>% mutate(col = str_remove(col, "...")) %>% mutate(col, col = as.integer(col))) 
list_concs_long <- do.call("rbind", list_concs_long) 
list_concs_long <- list_concs_long %>% mutate(key = paste0(plate,".",row,".",col)) %>% print(n = Inf)

# save for redundancy 
save(list_concs_long, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__pooling_concs.Rdata")


## Read sample metadata
## ====================

# read in
mdpath <- c("/Users/paul/Documents/OU_eDNA/191213_field_work/201130_sample_overview_updated.xlsx")
mdata <- readxl::read_excel(mdpath[1], range = "A1:L91", .name_repair = "universal")

# correct dates and erase superflous information
lubridate::date(mdata$sample_time) <- lubridate::date(mdata$sample_date)
mdata$sample_date <- NULL

save(mdata, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__mdata.Rdata")


## Merge data for further formatting
## =================================

big_table <- left_join(amplicon_positions,  primer_positions, by=c("key"), copy = TRUE, Keep = TRUE) %>% print(n = Inf)
big_table <- left_join(big_table,  primer_sequences, by=c("content.y" = "Name"), copy = TRUE, Keep = TRUE) %>% print(n = Inf)


## Write Excel Sheet compatible with Otago genomics sample submssion sheet
## =======================================================================

# - Extract unique needed rows and key
smpl <- big_table %>% select("key", "content.x") %>% distinct() %>% rename("Library name" = "content.x")
pr_rev <- big_table %>% filter(type.y == "r_primer") %>% select("key", "content.y", "Index") %>% rename("i7 index ID" = "content.y", "i7 index sequence" = "Index")
pr_fwd <- big_table %>% filter(type.y == "f_primer") %>% select("key", "content.y", "Index") %>% rename("i5 index ID" = "content.y", "i5 index sequence" = "Index")

# - Join needed primer rows by key 
smpl <- left_join(smpl, pr_rev, by = "key", copy = TRUE, keep = TRUE)
smpl <- left_join(smpl, pr_fwd, by = c("key.x" = "key"),copy = TRUE, keep = TRUE)  
smpl <- select(smpl, -c("key.x", "key.y", "key")) %>% filter(complete.cases(.))
smpl %>% print(n = Inf)

# write data to clibboard fr pasting into Excel
clip <- pipe("pbcopy", "w")                       
write.table(smpl, file=clip)                               
close(clip)


# Format big table for Qiime imports or other sequence processing
# ========================================================

# pre-cleaning column types
big_table[c("plate.y", "row.y", "col.y")] <- NULL

# merging sample metadata
# -----------------------

#  check input data
# amplicon_positions %>% print(n = Inf)
# mdata %>% print(n = Inf)

#  get keys in metadata table `mdata`
mdata <- left_join(mdata, amplicon_positions %>% select("key", "sample_name"), by = c("sample_name"), copy = TRUE) 

#  check output data
# mdata %>% print(n = Inf)

# merge metadata with all other data
big_table <- left_join(big_table,  mdata, by = c("key"), copy = TRUE) 
# big_table %>% print(n = Inf)

# merging concentration values
# ----------------------------

big_table <- left_join(big_table,  list_concs_long, by = c("key"), copy = TRUE) 
# big_table %>% print(n = Inf)

# correct column names
# --------------------

colnames(big_table)
big_table$type.x <- NULL
big_table$plate.x <- NULL
identical(big_table$row.x, big_table$row)
big_table$row.x <- NULL
identical(big_table$"content.x", big_table$"content.y")
big_table <- dplyr::rename(big_table, "pool_content" = "content.x" )
big_table <- dplyr::rename(big_table, "primer_label" = "content.y" )
identical(big_table$"col.x", big_table$"col")
big_table$"col.x" <- NULL
identical(big_table$"sample_name.y", big_table$"sample_name.x")
big_table <- dplyr::rename(big_table, "sample_name" = "sample_name.y" )
big_table$"sample_name.x" <- NULL
big_table <- dplyr::rename(big_table, "primer_direction" = "type.y")
colnames(big_table) <- tolower(colnames(big_table)) 

# filter one questionable barcode assignment - from 776 x 24 to 760 x 24 (keeping NA's)
big_table <- filter(big_table, is.na(sample_name) | sample_name != "PC191218-BABB-blk") %>% print(n = Inf)

# add information for negative and positive controls - to "sample_name" and "sample_type" column for down stream compatibility 
unique(big_table$pool_content)

unique(big_table$sample_name)
big_table <- big_table %>% mutate(sample_name = case_when(
                grepl("pttly: empty", pool_content) ~ "empty",
                grepl("ncntrl-pcr ", pool_content) ~ "ncntrl-pcr",
                grepl("ncntrl-pcr", pool_content) ~ "ncntrl-pcr",
                grepl("pttly: ncntrl-pcr", pool_content) ~ "ncntrl-pcr",
                grepl("pcntrol-blnd ", pool_content) ~ "pcntrol-blnd",
                grepl("pcntrol-zebra ", pool_content) ~ "pcntrol-zebra",
                grepl("xtr-blnk", pool_content) ~ "ncntrol-xtr",
                is.na(sample_name) ~ "empty",
                TRUE ~ as.character(sample_name)))
unique(big_table$sample_name)


unique(big_table$sample_type)
big_table <- big_table %>% mutate(sample_type = case_when(
               sample_type == "eDNA" ~ sample_type,
               sample_type == "blank" ~ sample_type,
               is.na(sample_type) ~ sample_name,
               TRUE ~ as.character(sample_type)))
unique(big_table$sample_type)


save(big_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__big_table.Rdata")


# Get file for deconvolution using cutadapt >=3.0
# ===============================================

# selcet relavant columns for clarity
demux_table <- big_table %>% select(key, primer_direction, sample_name, index, templateprimer)

# modify questionable characters for bash compatibility
demux_table <- demux_table %>% mutate(key =  gsub("\\.", "-", key)) %>% mutate(primer_direction =  gsub("\\_", "-", primer_direction))

# unite columns to get correct demultiplexing sequences and unique keys for spreading
demux_table <- demux_table %>% unite(seq_set, c(key,sample_name), remove = FALSE) %>% unite(primer_seq, c(index, templateprimer), sep = "", remove = FALSE)

# keep only columns relavant for cutadapat  and bash script
demux_table <- demux_table %>% select(seq_set, primer_seq, primer_direction) %>% print(n = Inf)

# spread data and rename columns 
demux_table <- demux_table %>% pivot_wider(id_cols = seq_set, names_from = primer_direction , values_from = primer_seq) %>% rename( "f_primer" = `f-primer` ) %>% rename( "r_primer" = `r-primer` )

# remove superflous information 
demux_table <- demux_table %>% filter(f_primer != "NANA" & r_primer != "NANA") %>% print(n = Inf)

#reverse complemnet 3' tagged preimer
strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
demux_table <- demux_table %>% mutate(r_primer, strReverse(chartr("acgtACGT", "tgcaTGCA", r_primer))) %>% rename(r_primer_rc = `strReverse(chartr("acgtACGT", "tgcaTGCA", r_primer))`)


# write four-column-text file for bash parsing
write_delim(demux_table, file = "/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/200_cutadapt_barcode_input.txt", delim = " ", append = FALSE, col_names = FALSE, quote_escape = "none", eol = "\n")

# also save R object in case needed later
save(demux_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__demux_table.Rdata")

