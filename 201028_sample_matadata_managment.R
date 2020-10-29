# Fiordldand project 28-10-2020
#
# Read cells from 
#   "/Users/paul/Documents/OU_eDNA/200128_lab_work/200907_plate_layouts.xlsx"
#   "/Users/paul/Documents/OU_eDNA/191031_primers/200302_Bunce_et_al_0000_MiFishEmod_single_step_primers.xlsx"  
# output appropriate format for
#   "/Users/paul/Documents/OU_eDNA/200128_lab_work/201028_Otago_Genomics_sample_info.xlsx"
#    and later possibly the mapping file


## Environment
## ===========

rm(list = ls())

library(tidyverse)


## Read in sample data
## ===================

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
save(amplicon_positions, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__amplicon_positions")

primer_positions <- bind_rows(u1pf, u1pr, u2pf, u2pr, e1pf, e1pr, e2pf, e2pr) %>% 
  mutate(col = str_remove(col, "...")) %>% 
  mutate(col, col = as.integer(col)) %>% 
  mutate(key = paste0(plate,".",row,".",col))
save(primer_positions, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__primer_positions")

## Read in primer data
## ===================

up5 <- readxl::read_excel(path[2], range = "E4:J11", .name_repair = "universal",  col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))
up7 <- readxl::read_excel(path[2], range = "E16:J32", .name_repair = "universal", col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))

ep5 <- readxl::read_excel(path[2], range = "E40:J47", .name_repair = "universal", col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))
ep7 <- readxl::read_excel(path[2], range = "E52:J68", .name_repair = "universal", col_names = c("Name", "Adapter",	"FluPad", "Index", "TemplatePrimer", "CompletePrimer"))

# get large data frame for further manipulation
# ---------------------------------------------

primer_sequences <- bind_rows(up5, up7, ep5, ep7)
save(primer_sequences, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__primers")


## Read amplicon concentration data
## ================================

# **pending**

## Read sample metadata
## ================================

# **pending**


## Merge data for further fomatting
## =======================================================================

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


## Write mapping file for Qiime imports or other sequenc processing
## ================================================================

# **pending**



primer_plates %>% filter(type %in% c("f_primer","fr_primer"))

primer_plates

pool_plates$type