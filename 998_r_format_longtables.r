#   **************************************************
#   * Combine, filter, and inspect long tables from  *
#   *   from eDNA data and BRUV observations         * 
#   **************************************************
#   26-Feb-2021, 1-Mar-2021, 7-Jul-2021, 15-Jul-2021
#   14-Dec-2022

# I. Load packages
# ================
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

rm(list = ls(all.names = TRUE))
gc()

library("tidyverse")   # tibbles, pipes, and more
library("magrittr")    # more pipes

library("robis")       # access OBIS data

library("readxl")      # read Excel files
library("openxlsx")    # write Excel tables

options(tibble.print_max = Inf) 

# II. Functions
# =============

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)

# III. Combine data 
# =================

# load data
# ----------

# eDNA data by  "/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_long_table.r"
edna_long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/221210_991_r_get_eDNA_long_table__eDNA_data.Rds")

# BRUV data by  "/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_BRUV_long_table.r"
bruv_long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/221214_996_r_get_BRUV_long_table__mh_bruv_obs.Rds")

# PUBL data by  "/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_PUBL_long_table.r"
publ_long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/221214_996_r_get_PUBL_long_table__publ_obs.Rds")

# stack data
# ----------
stack_long_table <- bind_rows(edna_long_table, bruv_long_table) 
dim(stack_long_table) #  was: 466 x 67 14.12.2022: 477 71
names(stack_long_table)

# IV. Format data
# ===============

# sorting and inspection for sanity reasons
stack_long_table <- stack_long_table %>% arrange(SET.ID) 

# need to have three or two UNIQ.REP.IDS, otherwise can't analyse data
#   7-Jul-21:
#     practically this is keeping only sets with complete eDNA and Bruv observations
#     also check /Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r, 
#     reused there to set the  UNIQ.REP.IDS, but no re-filtering applied
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% mutate(UNIQ.REP.IDS = n_distinct(REP.ID))
stack_long_table <- stack_long_table %>% filter(UNIQ.REP.IDS %in% c(2,3))

# among three or two UNIQ.REP.IDS needs to be on REP.ID = 3 (BRUV data) otherwise can't analyse
stack_long_table <- stack_long_table %>% filter(any(REP.ID == 3)) 

# sorting and inspection for sanity reasons
stack_long_table <- stack_long_table %>% arrange(desc(UNIQ.REP.IDS), SET.ID, REP.ID)

# works now:
# SET.ID: data set labels 
# REP.ID: replicate identifiers (of which "3", marks lines with BRUV data)

# rearrange columns
stack_long_table <-  stack_long_table %>% relocate(SET.ID,	REP.ID, LOC.NAME, INSIDE.RESERVE,  MH.GPS.LAT,	MH.PPS.LONG,  SUPERKINGDOM,	PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS, SPECIES)

# fill missing values for analysis
stack_long_table %>% group_by(SET.ID) %>% print(n = Inf)
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(LOC.NAME) 
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(INSIDE.RESERVE)
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(MH.GPS.LAT, .direction = c("downup"))
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(MH.PPS.LONG, .direction = c("downup"))

# V. Redefine areas inside and out side marine reserves
# ========================================================
# done using GID externally - check `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/210301_sample_map_overview.pdf`
# see subsequent script for site assignments - but not run here to maintain compatibility with MdL scripts
#  long_table <- long_table %>% mutate( RESERVE.GROUP = case_when(RESERVE.GROUP == "A" ~ "WJ", RESERVE.GROUP == "B" ~ "FF", RESERVE.GROUP == "C" ~ "LS"))
stack_long_table <- stack_long_table %>% mutate(RESERVE.GROUP = 
                                            case_when(SET.ID %in% c(21,22,23,24) ~ "A",
                                                      SET.ID %in% c(26,27,28,29) ~ "A",
                                                      SET.ID %in% c(11,12)       ~ "B",
                                                      SET.ID %in% c(17,18,19)    ~ "B",
                                                      SET.ID %in% c(7,8,9,10)    ~ "C",
                                                      SET.ID %in% c(1,3,4,5)     ~ "C"))

stack_long_table <- stack_long_table %>% mutate(RESERVE.GROUP.INSIDE = 
                                            case_when(SET.ID %in% c(21,22,23,24) ~ TRUE,
                                                      SET.ID %in% c(26,27,28,29) ~ FALSE,
                                                      SET.ID %in% c(11,12)       ~ TRUE,
                                                      SET.ID %in% c(17,18,19)    ~ FALSE,
                                                      SET.ID %in% c(7,8,9,10)    ~ FALSE,
                                                      SET.ID %in% c(1,3,4,5)     ~ TRUE))

long_table <- stack_long_table %>% mutate(SAMPLE.TYPE = case_when(REP.ID %in% c(3)   ~  "BRUV",
                                                                        REP.ID %in% c(1,2) ~  SAMPLE.TYPE))

# correct location name in fulls table
long_table <- long_table %>% mutate( RESERVE.GROUP = case_when(RESERVE.GROUP == "A" ~ "WJ", RESERVE.GROUP == "B" ~ "FF", RESERVE.GROUP == "C" ~ "LS"))

# remove three undetermined fish taxa from BRUV - which don't seem to be there anymore

long_table %<>% filter(ORDER != "NA")

# combine RESERVE.GROUP and RESERVE.GROUP.INSIDE to get six locations RESERVE.GROUP.LOCATION
long_table$RESERVE.GROUP
long_table$RESERVE.GROUP.INSIDE

long_table <- long_table %>% mutate( RESERVE.GROUP.LOCATION = 
  case_when(RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == TRUE  ~ "WJ MR",
            RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == FALSE ~ "WJ CTRL",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == TRUE  ~ "FF MR",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == FALSE ~ "FF CTRL",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == TRUE  ~ "LS MR",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == FALSE ~ "LS CTRL")
            )

long_table$RESERVE.GROUP.LOCATION
long_table %<>% relocate(RESERVE.GROUP.LOCATION)

# set BRUV observations to 1 for downstream generation of ASV presence column
long_table <- long_table %>% mutate(ABUNDANCE = 
   case_when(SAMPLE.TYPE == "BRUV" & is.na(ABUNDANCE) ~ 1,
             TRUE ~ ABUNDANCE)
             ) 

# rearrange columns 
long_table %<>% relocate(SET.ID,	REP.ID, SAMPLE.TYPE, LOC.NAME, MH.GPS.LAT,
  MH.PPS.LONG, RESERVE.GROUP,  RESERVE.GROUP.INSIDE, RESERVE.GROUP.LOCATION, SUPERKINGDOM,	
  PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES)

print(long_table)

# VI. Insert 15-Jul-2021 - add NCBI data in (had been lost)
# ==========================================================

# check relavant columns - NCBI.TAXID needs to be filled (again) for eDNA data
long_table |> select(SAMPLE.TYPE, ASV, NCBI.TAXID)


# load "blast_results_final" to get access to the "tax_id" column looked up previously
load(file="/Users/paul/Documents/OU_eDNA/201028_Robjects/221124_get_q2_tax-tab__blast-noenv_with-ncbi_taxonomy.Rdata")
blast_results_final |> select(iteration_query_def, tax_id)


# fill in missing NCBI tax strings 
long_table <- long_table |> 
  left_join( {blast_results_final |> select(iteration_query_def, tax_id) |> setNames(c("ASV", "NCBI.TAXID"))} , by = c("ASV")) |>
  unite(NCBI.TAXID, c(NCBI.TAXID.x, NCBI.TAXID.y), remove = TRUE, na.rm = TRUE) # |>
  # select(SAMPLE.TYPE, ASV, NCBI.TAXID)

# VIII. Insert 23-Jul-2021 - add Publication data
# ==============================================

# for superseded QGIS mapping in /Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210307_sample_map.qgz
# write.csv(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210301_997_r_format_longtables__analysis_input.csv")
write.csv(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/221214_998_r_format_longtables__analysis_input.csv")

long_table <- bind_rows( {long_table |> mutate(NCBI.TAXID = as.numeric(NCBI.TAXID))}, publ_long_table) 

# IX. write intermediate file 
# ===========================
dim(long_table) # was:        267 x 71 / 330 x 73 with PUBL data
                # 14.12.2022: 278 x 75 / 385 x 77 with PUBL data

# Workspace
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/221214_998_r_format_longtables__analysis_input__image.Rdata")
save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221214_998_r_format_longtables__analysis_input__image.Rdata")

# for verbosity
write.xlsx(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/221214_998_r_format_longtables__analysis_input.xlsx", asTable = FALSE, overwrite = TRUE)

# for previous analysis by MDL
# saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210301_997_r_format_longtables__analysis_input.Rds")
# saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210301_997_r_format_longtables__analysis_input.Rds")
saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/221214_998_r_format_longtables__analysis_input.Rds")
saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/221214_998_r_format_longtables__analysis_input.Rds")

