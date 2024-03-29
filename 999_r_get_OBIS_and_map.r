#   **************************************
#   * Get a nice map and add OBIS data   *
#   *                                    *  
#   **************************************
#   

# I. Load packages
# ================
rm(list = ls(all.names = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
gc()

library("tidyverse")   # tibbles, pipes, and more
library("magrittr") # more pipes

library("robis")       # access OBIS data

library("rgdal")       # create buffer around sampling areas
library("sf")          # create buffer around sampling areas
library("sp")

library("taxonomizr")  # get NCBI tax ids

library("taxize") # better then "taxonomizr" - may use and recode
Sys.setenv(ENTREZ_KEY="YourAPIKeyHere")
Sys.getenv("ENTREZ_KEY")

library("readxl")      # read Excel files
library("openxlsx")    # write Excel tables

options(tibble.print_max = Inf) 

# II. Functions
# =============

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)

# III. Load data 
# ==============

long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/221214_998_r_format_longtables__analysis_input.Rds")

# IV. prepare OBIS data and a map
# ================================
# use CRS 4326
# keep in mind 
#  long_table_dt_map <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("MH.GPS.LAT", "MH.PPS.LONG", "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]
# of
#  /Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r

# get a table with relevant columns for OBIS lookup
lt_obis_lookup <- long_table %>% 
  select("SET.ID", "MH.GPS.LAT", "MH.PPS.LONG",  "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION") %>%
  distinct() %>% arrange(SET.ID) 

# modify SET.ID 
#  from "98" Publication data set 
#  to "99" Obis data set to be added below (to each line)
lt_obis_lookup  <- lt_obis_lookup |> mutate(SET.ID = ifelse(SET.ID == 98, 99, SET.ID)) 

# add missing values for newly added publication data **use those below, again**
# use code to adjust great circle position
left_deg <- c( 0.150)
up_deg   <- c(-0.035)

lt_obis_lookup <- lt_obis_lookup |>  mutate(MH.GPS.LAT = ifelse(is.na(MH.GPS.LAT), (mean(na.omit(lt_obis_lookup$MH.GPS.LAT) + up_deg)), MH.GPS.LAT))
lt_obis_lookup <- lt_obis_lookup |>  mutate(MH.PPS.LONG = ifelse(is.na(MH.PPS.LONG), (mean(na.omit(lt_obis_lookup$MH.PPS.LONG) + left_deg) ), MH.PPS.LONG))
lt_obis_lookup <- lt_obis_lookup |>  mutate(RESERVE.GROUP = ifelse(is.na(RESERVE.GROUP), "FI", RESERVE.GROUP))
lt_obis_lookup <- lt_obis_lookup |>  mutate(RESERVE.GROUP.INSIDE = ifelse(is.na(RESERVE.GROUP.INSIDE), FALSE, RESERVE.GROUP.INSIDE))
lt_obis_lookup <- lt_obis_lookup |>  mutate(RESERVE.GROUP.LOCATION = ifelse(RESERVE.GROUP == "FI", "FI CTRL", RESERVE.GROUP.LOCATION))
lt_obis_lookup %>% print(n= Inf)

# get clean spatial data in degree units 
# ---------------------------------------
#  
# https://gis.stackexchange.com/questions/292327/creating-buffers-around-points-and-merging-with-spatialpolygonsdataframe-to-crea
# https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/understand-epsg-wkt-and-other-crs-definition-file-types/
# check cordinates sytems
EPSG <- rgdal::make_EPSG()
EPSG %>% filter(code %in% c(4326, 2193))
head(EPSG, 4)

# for buffer calculation - get sample coordinates as sf object with coordinates in EPSG 4326 / WGS84
lt_obis_lookup_sf <- st_as_sf(lt_obis_lookup,coords=c("MH.PPS.LONG","MH.GPS.LAT")) 
lt_obis_lookup_sf %<>% st_set_crs(4326) # set CRS to WGS84 

# for (test) maps  - get background layer(s) coordinates as sf object with coordinates in EPSG 4326 / WGS84
nzshp_hires = read_sf("/Users/paul/GIS/NZ_coast/NZ_Coast_isl.shp")
nzshp_hires_WGS84 <- st_transform(nzshp_hires, crs = 4326)
#   getting a lo-resolution map - https://gis.stackexchange.com/questions/243569/simplify-polygons-of-sf-object
nzshp_lores_WGS84 <- rmapshaper::ms_simplify(input = as(nzshp_hires_WGS84, 'Spatial')) %>% st_as_sf()

# create inset map for publication / define a bounding box around the field work area
#   https://geocompr.github.io/post/2019/ggplot2-inset-maps/
bb_fwork <- st_as_sfc(st_bbox(c(xmin = (166.5-0.1), xmax = (167.0+0.1), ymax = (-46.04-0.1), ymin = (-45.52+0.1)), crs = st_crs(4326)))
map_inset <- ggplot(data = nzshp_lores_WGS84) +
    geom_sf(fill = "grey93", color = "red", lwd = 0.5) +
    geom_sf(data = bb_fwork, fill = NA, color = "darkred", size = 1) +
    theme_void()

# get bounding box around sample groups for publication map
#   https://stackoverflow.com/questions/54696440/create-polygons-representing-bounding-boxes-for-subgroups-using-sf
calc_angle <- function(lon,lat) {
  cent_lon <- mean(lon)
  cent_lat <- mean(lat)
  ang <- atan2(lat - cent_lat, lon - cent_lon)
  return(ang)
}

# bounding boxes for all data but Fiordland data
bbox <- lt_obis_lookup %>%
  group_by(RESERVE.GROUP.LOCATION) %>% filter(RESERVE.GROUP.LOCATION != "FI CTRL") %>%
  summarise(xmin = min(MH.PPS.LONG) -0.01 ,ymin = min(MH.GPS.LAT) -0.01, xmax=max(MH.PPS.LONG) + 0.01,  ymax = max(MH.GPS.LAT) +0.01) %>%
  gather(x,lon,c('xmin','xmax')) %>%
  gather(y,lat,c('ymin','ymax')) %>%
  st_as_sf(coords=c('lon','lat'),crs=4326,remove=F) %>%
  group_by(RESERVE.GROUP.LOCATION) %>%
  mutate(angle = calc_angle(lon,lat)) %>%
  arrange(angle) %>%
  summarise(do_union=FALSE) %>%
  st_cast('POLYGON')

# check sf objects - looking ok so far
ggplot() +
    geom_sf(data = nzshp_hires_WGS84, fill = "lightgrey") + 
    geom_sf(data = lt_obis_lookup_sf, colour = "red") +
    geom_sf(data = bbox, fill = NA, color = "red") +   
    coord_sf(xlim = c((166.5-0.1), (167.0+0.1)), ylim = c((-46.04-0.1),(-45.52+0.1)), expand = FALSE) +
    theme_bw()

# get clean spatial data in local distance units 
# ----------------------------------------------

# for buffer generation re-project objects to local kms and check again
lt_obis_lookup_sf_loc <- lt_obis_lookup_sf %>% st_transform(crs = st_crs(2193))
nzshp_hires_WGS84_loc <- nzshp_hires_WGS84 %>% st_transform(crs = st_crs(2193))
bbox_loc <- bbox %>% st_transform(crs = st_crs(2193))

# calculate 2.5 km buffers - for inside cols
# calculate 2.5 km buffers - for outside cols - needs to be one great circle
# dosdy code - check results carefully - last column should have coordinates for a large circle

rows_inside  <- c(1:nrow(lt_obis_lookup_sf_loc)-1)
rows_outside <- c(nrow(lt_obis_lookup_sf_loc))

lt_obis_lookup_sf_buffer_loc <- st_buffer(lt_obis_lookup_sf_loc[rows_inside,  ], 2500)
lt_obis_lookup_sf_buffer_loc <- rbind(lt_obis_lookup_sf_buffer_loc, st_buffer(lt_obis_lookup_sf_loc[rows_outside, ], 38000))

# central coordinate for paper: 166.8948 -45.80689
lt_obis_lookup_sf_loc[rows_outside, ] %>% st_transform(crs = st_crs(2193))

# map to check object and for publication - unit is in km
#  check sf objects  - bounding box as defined per lt_obis_lookup_sf and 10 km in addition
#  inset grob in degrees, but positioned in kilometers
# 19-Sep-2023 - uncommnet all lines to get un-obfuscated locations
map_main <- ggplot() +
    geom_sf(data = nzshp_hires_WGS84_loc, fill = "lightgrey") + 
    geom_sf(data = lt_obis_lookup_sf_buffer_loc, fill = NA, colour = "red") + 
    # geom_sf(data = filter(lt_obis_lookup_sf_loc, RESERVE.GROUP != "FI"), fill = NA, colour = "darkred") + 
    geom_sf(data = bbox_loc, fill = NA, color = "darkred") +
    # stat_sf_coordinates(data = filter(lt_obis_lookup_sf_loc, RESERVE.GROUP != "FI"), aes(shape = RESERVE.GROUP), color = "darkred", size = 3) +
    # stat_sf_coordinates(data = filter(lt_obis_lookup_sf_loc, RESERVE.GROUP != "FI"), aes(shape = RESERVE.GROUP), color = "red", size = 1) +
    geom_sf_label(data=bbox_loc, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 4500, nudge_y = 4000) + 
    coord_sf( xlim = c((1125643-40000), (1125643+5000)), ylim = c((4909254-30000),(4909254+35000)), expand = FALSE) +
    annotation_custom(ggplotGrob(map_inset), xmin = 1125643-40000, xmax = 1125643-25000, ymin = 4909254+5000, ymax = 4909254+45000) + 
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.position=c(.9,.4), 
          legend.background = element_blank(), 
          legend.key=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

# see filename - map saving originally implemented in 
#  /Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r


# save image for combination with other plots
# saveRDS(map_main, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/999_r_get_OBIS_and_map__mapggplot.Rds")

# 19-Sep-2023save redacvted image for combination with other plots
saveRDS(map_main, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/999_r_get_OBIS_and_map__mapggplot_redacted.Rds")


ggsave("221219_999_r_summarize_results_fig1_draft.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 125, height = 175, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# back transform buffers to degree and wkt format for obis lookup
lt_obis_lookup_sf_buffer_d <- lt_obis_lookup_sf_buffer_loc %>% st_transform(crs = 4326)

# add column with buffers in Well-Known-Text format for OBIS lookup
lt_obis_lookup_sf_buffer_d %<>% mutate(BUFFER.WKT = st_as_text(geometry))

# simplify object to make life easier for later analyses
lt_obis_lookup <- st_drop_geometry(lt_obis_lookup_sf_buffer_d)

# V. fetch OBIS data 
#   (Wed Jul 24 10:45:50 NZST 2021)
#   (Mon Dec 19 17:57:47 CET 2022)
# ==================================================

# create nested data frame for OBIS lookup -  WKT polygons will go to $data slot
# renaming $data slot to BUFFER.WKT.NST as it is also a keyword
lt_obis_lookup %<>% select(SET.ID,BUFFER.WKT)
lt_obis_lookup %<>% group_by(SET.ID) %>% nest() %>% rename(BUFFER.WKT.NST = data)

# OBIS lookup function for each WKT POLYGON  - retrive all Chordata (taxon 1821)
get_obis <- function(BUFFER.WKT) occurrence(taxon = 1821, geometry = BUFFER.WKT)

# fill nested data frame with OBIS data
lt_obis_lookup %<>% mutate(OBIS = map(BUFFER.WKT.NST, get_obis)) # OBIS access may take some time

# Retrieved 5278 records of approximately 5278 (100%) - Mon Dec 19 17:57:47 CET 2022


lt_obis_lookup %>% print(n = Inf)                                # continue here after 24-Jul-2021
                                                                 # work space saved 

# create full table of raw results from nested table for cleaning and merging with other data
lt_obis_results <- unnest(lt_obis_lookup, cols = c(BUFFER.WKT.NST, OBIS))

# save table to save lookup time
saveRDS(lt_obis_results, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/221219_999_r_map_and_add_obis__lt_obis_results.Rds")
saveRDS(lt_obis_results, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/221219_999_r_map_and_add_obis__lt_obis_results.Rds")

lt_obis_results <- readRDS("/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/221219_999_r_map_and_add_obis__lt_obis_results.Rds")

# VI. format OBIS data (NCBI taxonomy addition)
# =============================================
# get clean NCBI conform definitions for SUPERKINGDOM PHYLUM CLASS ORDER FAMILY, GENUS, SPECIES, 
# get a "1" for each ABUNDANCE or count in buffer? / per SET.ID ?
# keep SET.ID for merging
# add columns for merging
# UNIQ.REP.IDS = ?
# REP.IDS = 4

# 50 species found (Wed Jul 24 10:45:50 NZST 2021)
# 55 species found (Mon Dec 19 17:57:47 CET 2022)
spc_in <- lt_obis_results |> ungroup()  |> select(species) |> filter(!is.na(species)) |> distinct() |> pull(species)

# Download NCBI annotations (to match taxonomy records of other surveys)
# ----------------------------------------------------------------------
# same as in /Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_PUBL_long_table.r

# sort strings
spc <- spc_in  |> sort() |> unique()   # for assembly of object 
gen <- gsub( " .*$", "", spc)          # for assembly of object 
gen_uniq <- gen |> sort() |> unique()  # for reporting only (below)


# download annotations - list of data frames
gen_ls <- classification(gen, db = "ncbi")
spc_ls <- classification(spc, db = "ncbi")
gen_uniq_list <- classification(gen_uniq, db = "ncbi")

# 38 of 55 species found in NCBI (and 22 not found in NCBI excluded from eDNA)
# ---------------------------------------------------------------------------

spc_in_nf <- tibble( SPECIES = c(
  "Acanthoclinus matti",
  "Aplidium coronum",
  "Aplidium phortax",
  "Aplidium powelli",
  "Botryllus stewartensis",
  "Cnemidocarpa bicornuta",
  "Cnemidocarpa nisiotis",
  "Didemnum inveteratum",
  "Diplosoma velatum",
  "Eudistoma circumvallatum",
  "Hemerocoetes monopterygius",
  "Lissoclinum notti",
  "Ritterella sigillinoides",
  "Synoicum kuranui",
  "Synoicum occidentalis",
  "Synoicum stewartense",
  "Trididemnum shawi"
  ))

gen_in_nf <- c("Ritterella")

# saving workspace manually
save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221219_999_r_map_and_add_obis__post_downloads.Rdata")
load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221219_999_r_map_and_add_obis__post_downloads.Rdata")

# get basic species lists
# -----------------------

# get list of tibbles and strip class attributes for row binding
spc_tibl_ls <- sapply(spc_ls, as_tibble)
gen_tibl_ls <- sapply(gen_uniq_list, as_tibble)

length(gen_tibl_ls) # 41 genera, 45 (20.12.2022)
length(spc_tibl_ls) # 50 spcies, 55 (22.12.2022)

# unnest objects - but keep list names (also species name duplicates) 
#  to later separate rows
spc_tibl <- bind_rows(spc_tibl_ls, .id = "column_label") |> select(-c(value))
gen_tibl <- bind_rows(gen_tibl_ls, .id = "column_label") |> select(-c(value))

# keep only needed taxonomic values
nrnk <- c("SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES")
spc_tibl <- spc_tibl |> mutate(rank = toupper(rank)) |> filter(rank %in% !! nrnk)
gen_tibl <- gen_tibl |> mutate(rank = toupper(rank)) |> filter(rank %in% !! nrnk)

# widen taxonomic information for downstream compatibility
spc <- spc_tibl |> pivot_wider(id_cols = column_label, names_from = rank,  values_from = name, names_sort = FALSE)
gen <- gen_tibl |> pivot_wider(id_cols = column_label, names_from = rank,  values_from = name, names_sort = FALSE)

# add future column NCBI.ID, drop now unneeded "column_label"
spc <- left_join(spc, spc_tibl |> filter(rank == "SPECIES") |> select(column_label, id)) |> select(-c(column_label)) |> rename("NCBI.TAXID" = "id") 
gen <- left_join(gen, gen_tibl |> filter(rank == "GENUS") |> select(column_label, id)) |> select(-c(column_label)) |> rename("NCBI.TAXID" = "id") 
  
# select genera that haven't been found on species level already - not needed
unique(spc$GENUS); unique(spc$SPECIES) # 34 Genera, 38 species fully resolved
gen_not_in_spc <- gen |> filter(GENUS %!in% unique(spc$GENUS))
# add species not found and later fill higher taxonomy columns up 
spc <- bind_rows(spc, spc_in_nf, gen) 

#  add missing species manually for which genus information is available
# ----------------------------------------------------------------------

# sort stuff
col_order <- c("SUPERKINGDOM", "PHYLUM",  "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")                     
spc <- spc |> mutate(GENUS = case_when( is.na(GENUS) == TRUE ~  gsub( " .*$", "", SPECIES), TRUE ~ as.character(GENUS)))
spc <- spc |> relocate(all_of(col_order)) |> arrange(across( rev(col_order[1:(length(col_order)-1)]) ))

# fill NA's that can be filled
spc <- spc |> group_by(GENUS)  |> fill(FAMILY, .direction = c("updown"))
spc <- spc |> group_by(FAMILY) |> fill(ORDER, .direction = c("updown"))
spc <- spc |> group_by(ORDER)  |> fill(CLASS, .direction = c("updown")) 
spc <- spc |> group_by(CLASS)  |> fill(PHYLUM, .direction = c("updown"))
spc <- spc |> group_by(PHYLUM) |> fill(SUPERKINGDOM, .direction = c("updown"))

# remove holes that can be
spc <- spc |> filter(!is.na(SPECIES)) |> distinct()

# plug holes manually 
spc <- spc |> mutate(FAMILY = ifelse(GENUS == "Ritterella", "Salpidae", FAMILY))
spc <- spc |> mutate(ORDER = ifelse(FAMILY == "Salpidae", "Salpida", ORDER))
spc <- spc |> mutate(CLASS = ifelse(ORDER == "Salpida", "Thaliacea", CLASS))

spc <- spc |> mutate(CLASS = ifelse(FAMILY == "Plesiopidae", "Actinopteri", CLASS))
spc <- spc |> mutate(ORDER = ifelse(FAMILY == "Plesiopidae", "Ovalentaria", ORDER))

# remove holes that can be
spc <- spc |> filter(!is.na(SPECIES)) |> distinct()
 

# add other variables for downstream compatibility
# -----------------------------------------------

spc <- spc |> mutate(NCBI.TAXID = ifelse(!is.na(NCBI.TAXID), NCBI.TAXID, as.character("0"))) 
spc <- spc |> mutate(NCBI.TAXID = as.numeric(NCBI.TAXID))
spc <- spc |> mutate(NCBI.TAXID.INC = ifelse(NCBI.TAXID == 0, TRUE, FALSE)) 
spc <- spc |> mutate(SAMPLE.TYPE = "OBIS") |> mutate(ABUNDANCE = 1)

save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221220_999_r_get_OBIS_and_map.Rdata")
load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221220_999_r_get_OBIS_and_map.Rdata")

# combine columns 
# ----------------

# left_join(lt_obis_results spc
lt_obis_truncated <- left_join(
  (lt_obis_results |> 
    select(SET.ID, id, species, depth) |> filter(!is.na(species)) |>
    rename(DEPTH.M = depth) |> mutate(DEPTH.M = ifelse(DEPTH.M < 0, NA,DEPTH.M)) |> 
    rename(ASV = id) |> rename(SPECIES = species)),
  spc)

head(lt_obis_truncated)

lt_obis_truncated <- lt_obis_truncated |> mutate(OBIS.OBS.PRES = 1)   

head(lt_obis_truncated)

lt_obis_truncated %<>% mutate(RESERVE.GROUP = 
  case_when(SET.ID %in% c(98, 99) ~ "FI",
            SET.ID %in% c(21,22,23,24) ~ "WJ",
            SET.ID %in% c(26,27,28,29) ~ "WJ",
            SET.ID %in% c(11,12)       ~ "FF",
            SET.ID %in% c(17,18,19)    ~ "FF",
            SET.ID %in% c(7,8,9,10)    ~ "LS",
            SET.ID %in% c(1,3,4,5)     ~ "LS")
            )

lt_obis_truncated |> select (SET.ID, RESERVE.GROUP) |> distinct()

lt_obis_truncated %<>% mutate(RESERVE.GROUP.INSIDE = 
  case_when(SET.ID %in% c(98, 99) ~ FALSE,
            SET.ID %in% c(21,22,23,24) ~ TRUE,
            SET.ID %in% c(26,27,28,29) ~ FALSE,
            SET.ID %in% c(11,12)       ~ TRUE,
            SET.ID %in% c(17,18,19)    ~ FALSE,
            SET.ID %in% c(7,8,9,10)    ~ FALSE,
            SET.ID %in% c(1,3,4,5)     ~ TRUE))

lt_obis_truncated %<>% mutate(RESERVE.GROUP.LOCATION = 
  case_when(RESERVE.GROUP == "FI" & RESERVE.GROUP.INSIDE == FALSE ~ "FI CTRL",
            RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == TRUE  ~ "WJ MR",
            RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == FALSE ~ "WJ CTRL",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == TRUE  ~ "FF MR",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == FALSE ~ "FF CTRL",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == TRUE  ~ "LS MR",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == FALSE ~ "LS CTRL")
            )

# stack data for subsequent analysis - check dimensions
dim(long_table) # 330 x 73, now 385 x 77

# stack data for subsequent analysis - correct type in previous data for succesful stacking
long_table %<>% mutate(NCBI.TAXID = as.numeric(NCBI.TAXID))
long_table %<>% mutate(DEPTH.M = as.numeric(DEPTH.M))

lt_obis_truncated %<>% mutate(NCBI.TAXID = as.numeric(NCBI.TAXID))
lt_obis_truncated %<>% mutate(DEPTH.M = as.numeric(DEPTH.M))

# stack data for subsequent analysis
long_table %<>% bind_rows(long_table, lt_obis_truncated)
dim(long_table) # 6973   74, now 7101   78

long_table %<>% mutate(RESERVE.GROUP = 
  case_when(SET.ID %in% c(98, 99) ~ "FI",
            SET.ID %in% c(21,22,23,24) ~ "WJ",
            SET.ID %in% c(26,27,28,29) ~ "WJ",
            SET.ID %in% c(11,12)       ~ "FF",
            SET.ID %in% c(17,18,19)    ~ "FF",
            SET.ID %in% c(7,8,9,10)    ~ "LS",
            SET.ID %in% c(1,3,4,5)     ~ "LS")
            )

long_table |> select (SET.ID, RESERVE.GROUP) |> distinct()

long_table %<>% mutate(RESERVE.GROUP.INSIDE = 
  case_when(SET.ID %in% c(98, 99) ~ FALSE,
            SET.ID %in% c(21,22,23,24) ~ TRUE,
            SET.ID %in% c(26,27,28,29) ~ FALSE,
            SET.ID %in% c(11,12)       ~ TRUE,
            SET.ID %in% c(17,18,19)    ~ FALSE,
            SET.ID %in% c(7,8,9,10)    ~ FALSE,
            SET.ID %in% c(1,3,4,5)     ~ TRUE))

long_table |> select (SET.ID, RESERVE.GROUP, RESERVE.GROUP.INSIDE) |> distinct()

long_table %<>% mutate(RESERVE.GROUP.LOCATION = 
  case_when(RESERVE.GROUP == "FI" & RESERVE.GROUP.INSIDE == FALSE ~ "FI CTRL",
            RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == TRUE  ~ "WJ MR",
            RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == FALSE ~ "WJ CTRL",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == TRUE  ~ "FF MR",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == FALSE ~ "FF CTRL",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == TRUE  ~ "LS MR",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == FALSE ~ "LS CTRL")
            )

long_table |> select (SET.ID, RESERVE.GROUP.LOCATION, RESERVE.GROUP, RESERVE.GROUP.INSIDE) |> distinct()

# define unique observations by technique
long_table %<>% mutate(BRUV.OBS.PRES = case_when(SAMPLE.TYPE == "BRUV" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table %<>% mutate(EDNA.OBS.PRES = case_when(SAMPLE.TYPE == "eDNA" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table %<>% mutate(OBIS.OBS.PRES = case_when(SAMPLE.TYPE == "OBIS" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table %<>% mutate(PUBL.OBS.PRES = case_when(SAMPLE.TYPE == "PUBL" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))

long_table |> select (SET.ID, SAMPLE.TYPE, LOC.NAME, RESERVE.GROUP.LOCATION, RESERVE.GROUP, RESERVE.GROUP.INSIDE) |>
  arrange(SET.ID) |> distinct()

long_table <- long_table |> mutate(LOC.NAME = ifelse(SET.ID == 99, "Fiordland", LOC.NAME)) 

long_table |> select (SET.ID, SAMPLE.TYPE, LOC.NAME, RESERVE.GROUP.LOCATION, RESERVE.GROUP, RESERVE.GROUP.INSIDE) |>
  arrange(SET.ID) |> distinct()

# fill missing values for analysis
long_table %>% group_by(SET.ID) %>% print(n = Inf)
long_table %<>% group_by(SET.ID) %>% fill(LOC.NAME)
long_table %<>% group_by(SET.ID) %>% fill(INSIDE.RESERVE)
long_table %<>% group_by(SET.ID) %>% fill(MH.GPS.LAT, .direction = c("downup"))
long_table %<>% group_by(SET.ID) %>% fill(MH.PPS.LONG, .direction = c("downup"))
long_table %<>% group_by(SPECIES) %>%  fill(NCBI.TAXID.INC,  .direction = c("downup"))
long_table %<>% group_by(SPECIES) %>%  fill(NCBI.TAXID, .direction = c("downup"))
long_table %<>% group_by(SPECIES) %>%  fill(NCBI.LEVEL,  .direction = c("downup"))
long_table %<>% ungroup(SET.ID)

# rearrange columns as in previous data combination
long_table %<>% relocate(SET.ID,REP.ID, SAMPLE.TYPE, LOC.NAME, MH.GPS.LAT,
  MH.PPS.LONG, RESERVE.GROUP,  RESERVE.GROUP.INSIDE, RESERVE.GROUP.LOCATION, SUPERKINGDOM,	
  PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES, NCBI.TAXID, NCBI.TAXID.INC, NCBI.LEVEL)

# check again
long_table %>% select(SET.ID,REP.ID, SAMPLE.TYPE, LOC.NAME, MH.GPS.LAT,
  MH.PPS.LONG, RESERVE.GROUP,  RESERVE.GROUP.INSIDE, RESERVE.GROUP.LOCATION, SUPERKINGDOM,	
  PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES, NCBI.TAXID, NCBI.TAXID.INC, NCBI.LEVEL) %>%
  distinct()

long_table %<>% mutate(NCBI.LEVEL = ifelse(is.na(NCBI.LEVEL), "species", NCBI.LEVEL))
long_table %<>% mutate(NCBI.TAXID.INC = ifelse(is.na(NCBI.TAXID.INC), FALSE, NCBI.TAXID.INC))

#   7-Jul-21: 
#     in ~/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r
#     needed to to have 2 two or 3 UNIQ.REP.IDS
#     practically this was keeping only sets with complete eDNA and Bruv observations
#     redefine UNIQ.REP.IDS here 
#     but can't be used for filtering in the old fashion anymore 
long_table %<>% group_by(SET.ID) %>% mutate(UNIQ.REP.IDS = n_distinct(REP.ID))

# VIII. Check data completeness and citations
# ==========================================

long_table %<>% mutate(LOC.NAME = ifelse(RESERVE.GROUP == "FI", "Fiordland", LOC.NAME))

# check table
long_table |> select(SET.ID, SAMPLE.TYPE, LOC.NAME, MH.GPS.LAT, MH.PPS.LONG,
  RESERVE.GROUP, RESERVE.GROUP.INSIDE, RESERVE.GROUP.LOCATION, SUPERKINGDOM,	
  PHYLUM,	CLASS, ORDER,	FAMILY,	GENUS,	SPECIES) |> distinct()

# plug holes

long_table <- long_table |> mutate(SUPERKINGDOM = ifelse(SPECIES == "Puffinus griseus", "Eukaryota", SUPERKINGDOM))
long_table <- long_table |> mutate(PHYLUM       = ifelse(SPECIES == "Puffinus griseus", "Chordata", PHYLUM))
long_table <- long_table |> mutate(CLASS        = ifelse(SPECIES == "Puffinus griseus", "Aves", CLASS))
long_table <- long_table |> mutate(ORDER        = ifelse(SPECIES == "Puffinus griseus", "Procellariiformes", ORDER))
long_table <- long_table |> mutate(FAMILY       = ifelse(SPECIES == "Puffinus griseus", "Ardenna", FAMILY))
long_table <- long_table |> mutate(GENUS        = ifelse(SPECIES == "Puffinus griseus", "Ardenna", GENUS))
long_table <- long_table |> mutate(SPECIES      = ifelse(SPECIES == "Ardenna grisea", "Ardenna", SPECIES))

glimpse(long_table)
unique(long_table$SAMPLE.TYPE)

long_table <- long_table |> mutate(SAMPLE.TYPE = as.factor(SAMPLE.TYPE))
unique(long_table$SAMPLE.TYPE)
unique(long_table$SET.ID)
glimpse(long_table)

save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221220_999_r_get_OBIS_and_map_bug_chase.Rdata")
load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221220_999_r_get_OBIS_and_map_bug_chase.Rdata")

# this line needs to be correct so that sampl types don't get tunred into integers
long_table <- long_table |> mutate(SAMPLE.TYPE = case_when(
  SET.ID == 99 ~ as.factor("OBIS"),
  SET.ID != 99 ~ SAMPLE.TYPE)
)

unique(long_table$SAMPLE.TYPE)

# check table
long_table |> select(SET.ID, SAMPLE.TYPE, LOC.NAME, MH.GPS.LAT, MH.PPS.LONG,
  RESERVE.GROUP, RESERVE.GROUP.INSIDE, RESERVE.GROUP.LOCATION, SUPERKINGDOM,	
  PHYLUM,	CLASS, ORDER,	FAMILY,	GENUS,	SPECIES) |> distinct()

long_table %<>% filter(!is.na(SPECIES)) 

# safe table with citation info
data_citations <- lt_obis_results %>% ungroup() %>% select(bibliographicCitation) %>% filter(!is.na(bibliographicCitation)) %>% 
   distinct() %>% arrange(bibliographicCitation)#  %>% print(n = Inf) 

# write.xlsx(data_citations, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210707_OBIS_data_citations.xlsx", asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(data_citations, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/221220_OBIS_data_citations.xlsx", asTable = TRUE, overwrite = TRUE)


# check data completeness - preformatting
left_join(ungroup(select(lt_obis_results, id)), ungroup(select(long_table, ASV)), by = c("id" = "ASV"))
OBIS_records <- ungroup(lt_obis_results) %>% select(id) %>% distinct() %>% rename(OBIS_record = id)
OBIS_records %<>%  mutate(used = case_when(OBIS_record %in%  long_table$ASV ~ TRUE, TRUE ~ FALSE))

# and numerical summaries for manuscript
nrow(OBIS_records)     # 5263, now 5278
sum(OBIS_records$used) # 5186, now 5199
sum(OBIS_records$used) / nrow(OBIS_records) # 0.9853696, now 0.9850322

# XI. Save results and work space
# ==============================

# for verbosity
# write.xlsx(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/998_r_map_and_add_obis__full_data_raw.xlsx", asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/999_r_map_and_add_obis__full_data_raw.xlsx", asTable = TRUE, overwrite = TRUE)

# for superseded QGIS mapping in /Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210307_sample_map.qgz
# write.csv(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/998_r_map_and_add_obis__full_data_raw.csv")
write.csv(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/999_r_map_and_add_obis__full_data_raw.csv")

# saving work space manually
# save.image("/Users/paul/Documents/OU_eDNA/260705_r_workspaces/210705_998_r_map_and_add_obis__end.Rdata")
# save.image("/Users/paul/Documents/OU_eDNA/201028_Robjects/210705_998_r_map_and_add_obis__end.Rdata")
save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221220_999_r_map_and_add_obis__end.Rdata")
save.image("/Users/paul/Documents/OU_eDNA/201028_Robjects/221220_999_r_map_and_add_obis__end.Rdata")

# for subsequent analyses
# saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_map_and_add_obiss__full_data_raw.Rds")
# saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/998_r_map_and_add_obiss__full_data_raw.Rds")

saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/999_r_map_and_add_obis__full_data_raw.Rds")
saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/999_r_map_and_add_obis__full_data_raw.Rds")
