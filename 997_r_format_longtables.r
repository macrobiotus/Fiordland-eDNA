#   **************************************************
#   * Combine, filter, and inspect long tables from  *
#   *   from eDNA data and BRUV observations         * 
#   **************************************************
#   26-Feb-2021

# I. Load packages
# ================
rm(list = ls(all.names = TRUE))
gc()

library("readxl")      # read Excel files
library("openxlsx")    # write Excel tables
library("tidyverse")   # tibbles, pipes, and more

# II. Functions
# =============

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)

# III. Combine data and format
# ===========================

# load data
# ----------

# eDNA data by  "/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_long_table.r"
edna_long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210210_990_r_get_eDNA_phyloseq__clean-marine-eDNA.Rds")

# BRUV data by  "/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_BRUV_long_table.r"
bruv_long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210226_995_r_get_BRUV_long_table__mh_bruv_obs.Rds")

# stack data
# ----------

bind_rows(edna_long_table, bruv_long_table)







mh_obs_raw_path <- c("/Users/paul/Documents/OU_eDNA/191213_field_work/210225_MH_bruv_data_machine_readable_long.csv")
mh_obs_raw <- read_csv("/Users/paul/Documents/OU_eDNA/191213_field_work/210225_MH_bruv_data_machine_readable_long.csv")

