#   **************************************************
#   * Combine, filter, and inspect long tables from  *
#   *   from eDNA data and BRUV observations         * 
#   **************************************************
#   26-Feb-2021, 1-Mar-2021

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

# III. Combine data 
# =================

# load data
# ----------

# eDNA data by  "/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_long_table.r"
edna_long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210301_990_r_get_eDNA_long_table__eDNA_data.Rds")

# BRUV data by  "/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_BRUV_long_table.r"
bruv_long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210226_995_r_get_BRUV_long_table__mh_bruv_obs.Rds")

# stack data
# ----------
stack_long_table <- bind_rows(edna_long_table, bruv_long_table) #  466 x 67


# IV. Format data  
# ===============

# sorting and inspection for sanity reasons
stack_long_table <- stack_long_table %>% arrange(SET.ID) 

# need to have three or two UNIQ.REP.IDS, otherwise can't analyse data
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
stack_long_table <-  stack_long_table %>% relocate(SET.ID,	REP.ID, LOC.NAME, INSIDE.RESERVE,  MH.GPS.LAT,	MH.PPS.LONG,  SUPERKINGDOM,	PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES)

# fill missing values for analysis
stack_long_table %>% group_by(SET.ID) %>% print(n = Inf)
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(LOC.NAME) %>% print(n = Inf)
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(INSIDE.RESERVE) %>% print(n = Inf)
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(MH.GPS.LAT, .direction = c("downup")) %>% print(n = Inf)
stack_long_table <- stack_long_table %>% group_by(SET.ID) %>% fill(MH.PPS.LONG, .direction = c("downup")) %>% print(n = Inf)


# IV. Redefine areas inside and out side marine reserves
# ========================================================
# done using GID externally - check `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/210301_sample_map_overview.pdf`

stack_long_table <- stack_long_table %>% mutate(RESERVE.GROUP = 
                                            case_when(SET.ID %in% c(21,22,23,24) ~ "A",
                                                      SET.ID %in% c(26,27,28,29) ~ "A",
                                                      SET.ID %in% c(11,12)       ~ "B",
                                                      SET.ID %in% c(17,18,19)    ~ "B",
                                                      SET.ID %in% c(7,8,9,10)    ~ "C",
                                                      SET.ID %in% c(1,3,4,5)     ~ "C"))

stack_long_table <- stack_long_table %>% mutate(RESERVE.GROUP.INSIDE = 
                                            case_when(SET.ID %in% c(21,22,23,24) ~ TRUE,
                                                      SET.ID %in% c(26,27,28,29)    ~ FALSE,
                                                      SET.ID %in% c(11,12)       ~ TRUE,
                                                      SET.ID %in% c(17,18,19)    ~ FALSE,
                                                      SET.ID %in% c(7,8,9,10)    ~ FALSE,
                                                      SET.ID %in% c(1,3,4,5)     ~ TRUE))

stack_long_table <- stack_long_table %>% mutate(SAMPLE.TYPE = case_when(REP.ID %in% c(3)   ~  "BRUV",
                                                                        REP.ID %in% c(1,2) ~  SAMPLE.TYPE))


# rearrange columns 
stack_long_table <-  stack_long_table %>% relocate(SET.ID,	REP.ID, SAMPLE.TYPE, LOC.NAME, MH.GPS.LAT,	MH.PPS.LONG, RESERVE.GROUP,  RESERVE.GROUP.INSIDE,  SUPERKINGDOM,	PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES)

# V. write final file 
# ===================

save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210301_997_r_format_longtables__analysis_input__image.Rdata")
saveRDS(stack_long_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210301_997_r_format_longtables__analysis_input.Rds")
saveRDS(stack_long_table, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210301_997_r_format_longtables__analysis_input.Rds")
write.xlsx(stack_long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210301_997_r_format_longtables__analysis_input.xlsx", asTable = FALSE)






