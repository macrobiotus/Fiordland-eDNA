# **********************************************
# * Create, filter, and write Physloseq object *
#            from BRUV observations            * 
# **********************************************
# 26-Aug-2022, 22-Feb-2021, 23-Feb-2021, 25-Feb-2021

# I. Load packages
# ================
rm(list = ls(all.names = TRUE))
gc()

library("readxl")     # read Excel files
library("openxlsx")   # write Excel tables
library("janitor")    # clean column names
library("taxonomizr") # query taxon names
library("stringi")     # random strings
library("tidyverse")   # 

# II. Functions
# =============

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)


# III. Get species data for MH observations
# =========================================

# load data
# ----------

mh_obs_raw_path <- c("/Users/paul/Documents/OU_eDNA/191213_field_work/210225_MH_bruv_data_machine_readable_long.csv")
mh_obs_raw <- read_csv("/Users/paul/Documents/OU_eDNA/191213_field_work/210225_MH_bruv_data_machine_readable_long.csv")

# taxonomy database lookup and inspection 
# ----------------------------------------

# 10-Aug-2020: check external hard drive for readily created database files, if unavailable run 
# prepareDatabase(sqlFile = "accessionTaxa.sql", tmpDir = "/Users/paul/Sequences/References/taxonomizR/", vocal = TRUE) # takes a very long time - avoid by reloading full object from disk

# function for mutate to use taxonomic IDs and add taxonomy strings
get_strng <- function(x) {getTaxonomy(x,"/Volumes/HGST1TB/Users/paul/Sequences/References/taxonomizR/accessionTaxa.sql")}

# look up taxonomy table - takes a long time, needs external database.
tax_table <- as_tibble(get_strng(mh_obs_raw$NCBI.TAXID), rownames = "NCBI.TAXID") %>% mutate(NCBI.TAXID = as.numeric(NCBI.TAXID))
tax_table %>% distinct() %>%  print(n = Inf)

# add missing entries manually - indicate which rows need correction 
tax_table <- tax_table %>% mutate(
  NCBI.TAXDB.INC = case_when(
    NCBI.TAXID %in% c("232418", "2696622", "2696554", "206114", "78394") ~ TRUE,
    NCBI.TAXID %!in% c("232418", "2696622", "2696554", "206114", "78394") ~ FALSE,
    )
  )  %>% print(n = Inf)


# taxonomy database lookup correction 
# -----------------------------------

# in first table get unique rows, to accomodate "rows_upsert" et al. 
tax_table_distinct <- tax_table %>% distinct() %>% print(n = Inf)

# in second table re-define ill-defined rows
corrected_taxa <- tribble(
  ~NCBI.TAXID, ~superkingdom,   ~phylum,           ~class,              ~order,           ~family,             ~genus,                      ~species, ~NCBI.TAXDB.INC,
       206114,   "Eukaryota", "Chordata",    "Actinopteri",     "Blenniiformes",  "Tripterygiidae",     "Forsterygion",      "Forsterygion maryannae",            TRUE,
       232418,   "Eukaryota", "Chordata", "Chondrichthyes", "Carcharhiniformes",  "Scyliorhinidae",  "Cephaloscyllium",   "Cephaloscyllium isabellum",            TRUE,
      2696554,   "Eukaryota", "Chordata",    "Actinopteri",       "Perciformes",      "Serranidae",   "Hypoplectrodes",       "Hypoplectrodes huntii",            TRUE,
      2696622,   "Eukaryota", "Chordata",    "Actinopteri",       "Labriformes",        "Labridae",       "Notolabrus",          "Notolabrus cinctus",            TRUE,
        78394,   "Eukaryota", "Chordata",         "Myxini",      "Myxiniformes",       "Myxinidae",       "Eptatretus",        "Eptatretus cirrhatus",            TRUE,
         7898,   "Eukaryota", "Chordata",    "Actinopteri",                  NA,                NA,                 NA,                             NA,           FALSE,
)


tax_table_distinct_curated <- rows_upsert(tax_table_distinct, corrected_taxa, copy = TRUE, in_place = FALSE)  %>% print(n = Inf)
names(tax_table_distinct_curated) <- toupper(names(tax_table_distinct_curated))

# IV. ADD species data to  MH observations
# =========================================

mh_obs <- right_join(mh_obs_raw, tax_table_distinct_curated, by = "NCBI.TAXID", KEEP = FALSE, copy = TRUE)


# V. Export molten long table for merging with eDNA data 
# ======================================================

# save or load molten state 
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210226_995_r_get_BRUV_long_table__image.Rdata")
save(mh_obs, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210226_995_r_get_BRUV_long_table__mh_bruv_obs.Rds")
save(mh_obs, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210226_995_r_get_BRUV_long_table__mh_bruv_obs.Rds")
