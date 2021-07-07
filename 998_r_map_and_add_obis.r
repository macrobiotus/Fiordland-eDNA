#   **************************************
#   * Get a nice map and add OBIS data   *
#   *                                    *  
#   **************************************
#   

# I. Load packages
# ================
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse")   # tibbles, pipes, and more
library("magrittr") # more pipes

library("robis")       # access OBIS data

library("sf")          # create buffer around sampling areas
library("sp")

library("taxonomizr")  # get NCBI tax ids


library("readxl")      # read Excel files
library("openxlsx")    # write Excel tables


# II. Functions
# =============

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)

# III. Load data 
# ==============

long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210301_997_r_format_longtables__analysis_input.Rds")

# IV. prepare OBIS data and a map
# ================================
# use CRS 4326
# keep in mind 
#  long_table_dt_map <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("MH.GPS.LAT", "MH.PPS.LONG", "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]
# of
#  /Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r

# get a table with relevant columns  for OBIS lookup
lt_obis_lookup <- long_table %>% 
  select("SET.ID", "MH.GPS.LAT", "MH.PPS.LONG",  "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION") %>%
  distinct() %>% arrange(SET.ID) 
lt_obis_lookup %>% print(n= Inf)

# get clean spatial data in degree units 
# ---------------------------------------
#  
# https://gis.stackexchange.com/questions/292327/creating-buffers-around-points-and-merging-with-spatialpolygonsdataframe-to-crea
# https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/understand-epsg-wkt-and-other-crs-definition-file-types/
# check cordinates sytems
EPSG <- rgdal::make_EPSG()
EPSG %>% filter(code == 4326)
head(EPSG, 3)

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

bbox <- lt_obis_lookup %>%
  group_by(RESERVE.GROUP.LOCATION) %>%
  summarise(xmin = min(MH.PPS.LONG) -0.01 ,ymin = min(MH.GPS.LAT) -0.01, xmax=max(MH.PPS.LONG) +0.01,  ymax = max(MH.GPS.LAT) +0.01) %>%
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
# ---------------------------------------

# for buffer generation re-project objects to local kms and check again
lt_obis_lookup_sf_loc <- lt_obis_lookup_sf %>% st_transform(crs = st_crs("+proj=utm +zone=58G +datum=WGS84 +units=km"))
nzshp_hires_WGS84_loc <- nzshp_hires_WGS84 %>% st_transform(crs = st_crs("+proj=utm +zone=58G +datum=WGS84 +units=km"))
bbox_loc <- bbox %>% st_transform(crs = st_crs("+proj=utm +zone=58G +datum=WGS84 +units=km"))

# calculate 2.5 km buffers
lt_obis_lookup_sf_buffer_loc <- st_buffer(lt_obis_lookup_sf_loc, 2.5)

# map to check object and for publication - unit is in km
#  check sf objects  - bounding box as defined per lt_obis_lookup_sf and 10 km in addition
#  inset grob in degrees, but positioned in kilometers

map_main <- ggplot() +
    geom_sf(data = nzshp_hires_WGS84_loc, fill = "lightgrey") + 
    geom_sf(data = lt_obis_lookup_sf_buffer_loc, fill = NA, colour = "red") + 
    geom_sf(data = lt_obis_lookup_sf_loc, fill = NA, colour = "darkred") + 
    geom_sf(data = bbox_loc, fill = NA, color = "darkred") +
    stat_sf_coordinates(data = lt_obis_lookup_sf_loc, aes(shape = RESERVE.GROUP), color = "darkred", size = 3) +
    stat_sf_coordinates(data = lt_obis_lookup_sf_loc, aes(shape = RESERVE.GROUP), color = "red", size = 1) +
    geom_sf_label(data=bbox_loc, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6) + 
    coord_sf( xlim = c((619.6011-10), (653.8977+10)), ylim = c((-5100.241-10),(-5042.894+10)) , expand = FALSE) +
    annotation_custom(ggplotGrob(map_inset), xmin = 610, xmax = 625, ymin = -5065, ymax = -5025) + 
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.position=c(.9,.1), 
          legend.background = element_blank(), 
          legend.key=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

# see filename - map saving originally implemented in 
#  /Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r

ggsave("210401_998_r_summarize_results_fig1_draft.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 125, height = 175, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# back transform buffers to degree and wkt format for obis lookup
lt_obis_lookup_sf_buffer_d <- lt_obis_lookup_sf_buffer_loc %>% st_transform(crs = 4326)

# add column with buffers in Well-Known-Text format for OBIS lookup
lt_obis_lookup_sf_buffer_d %<>% mutate(BUFFER.WKT = st_as_text(geometry))

# simplify object to make life easier for later analyses
lt_obis_lookup <- st_drop_geometry(lt_obis_lookup_sf_buffer_d)

# V. fetch OBIS data (Wed Jul  7 10:45:50 NZST 2021)
# ==================================================

# create nested data frame for OBIS lookup -  WKT polygons will go to $data slot
# renaming $data slot to BUFFER.WKT.NST as it is also a keyword
lt_obis_lookup %<>% select(SET.ID,BUFFER.WKT)
lt_obis_lookup %<>% group_by(SET.ID) %>% nest() %>% rename(BUFFER.WKT.NST = data)

# OBIS lookup function for each WKT POLYGON 
get_obis <- function(BUFFER.WKT) occurrence(taxon = 2, geometry = BUFFER.WKT)

# fill nested data frame with OBIS data
lt_obis_lookup %<>% mutate(OBIS = map(BUFFER.WKT.NST, get_obis) )
lt_obis_lookup %>% print(n = Inf)

# create full table of raw results from nested table for cleaning and merging with other data
lt_obis_results <- unnest(lt_obis_lookup, cols = c(BUFFER.WKT.NST, OBIS))

# save table to save lookup time
saveRDS(lt_obis_results, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210705_998_r_map_and_add_obis__lt_obis_results.Rds")
saveRDS(lt_obis_results, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210705_998_r_map_and_add_obis__lt_obis_results.Rds")

# VI. format OBIS data (NCBI taxonomy addition)
# =============================================
# get clean NCBI conform definitions for SUPERKINGDOM PHYLUM CLASS ORDER FAMILY, GENUS, SPECIES, 
# get a "1" for each ABUNDANCE or count in buffer? / per SET.ID ?
# keep SET.ID for merging
# add columns for merging
# UNIQ.REP.IDS = ?
# REP.IDS = 4

# getting clean taxonomy - look up NCBI taxonomy ids
#  - look up NCBI strings as for other datasets 
#  - use first taxonomy strings, if more then one dicovered
lt_obis_results %<>% mutate(NCBI.TAXID = getId(scientificName ,"/Volumes/HGST1TB/Users/paul/Sequences/References/taxonomizR/accessionTaxa.sql", onlyScientific = TRUE)) 
lt_obis_results %<>% mutate(NCBI.TAXID = as.numeric(gsub(",.*$", "", NCBI.TAXID)))

# get NCBI taxonomy strings
get_strng <- function(x) {getTaxonomy(x,"/Volumes/HGST1TB/Users/paul/Sequences/References/taxonomizR/accessionTaxa.sql")}
ncbi_strings <- as_tibble(get_strng(unique(lt_obis_results$NCBI.TAXID)), rownames = "NCBI.TAXID") %>% 
                   mutate(NCBI.TAXID= as.numeric(NCBI.TAXID)) %>% filter(!is.na(NCBI.TAXID)) %>% 
                   rename_all(toupper) 
# ncbi_strings %>% print(n=Inf)

# add taxonomy strings to OBIS results
lt_obis_results %<>% left_join(ncbi_strings)


# VII. format OBIS data (match with previous tales and stack)
# ===========================================================

# save results for merging with other data -
# set values to match Excel table
#  "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210301_997_r_format_longtables__analysis_input.xlsx"
lt_obis_truncated <- lt_obis_results %>% 
  select(SET.ID, NCBI.TAXID, SUPERKINGDOM, PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES, id, depth ) %>% 
  filter(rowSums(across(c(SUPERKINGDOM, PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES), ~ !is.na(.))) > 0) %>% 
  rename(DEPTH.M = depth) %>% mutate(DEPTH.M = ifelse(DEPTH.M < 0, NA,DEPTH.M)) %>% 
  rename(ASV = id) %>% 
  add_column(UNIQ.REP.IDS = 1) %>% 
  add_column(ABUNDANCE = 1) %>% 
  add_column(REP.ID = 4) %>% 
  add_column(SAMPLE.TYPE = "OBIS") 

lt_obis_truncated %<>% mutate(RESERVE.GROUP = 
  case_when(SET.ID %in% c(21,22,23,24) ~ "WJ",
            SET.ID %in% c(26,27,28,29) ~ "WJ",
            SET.ID %in% c(11,12)       ~ "FF",
            SET.ID %in% c(17,18,19)    ~ "FF",
            SET.ID %in% c(7,8,9,10)    ~ "LS",
            SET.ID %in% c(1,3,4,5)     ~ "LS")
            )

lt_obis_truncated %<>% mutate(RESERVE.GROUP.INSIDE = 
  case_when(SET.ID %in% c(21,22,23,24) ~ TRUE,
            SET.ID %in% c(26,27,28,29) ~ FALSE,
            SET.ID %in% c(11,12)       ~ TRUE,
            SET.ID %in% c(17,18,19)    ~ FALSE,
            SET.ID %in% c(7,8,9,10)    ~ FALSE,
            SET.ID %in% c(1,3,4,5)     ~ TRUE))

lt_obis_truncated %<>% mutate(RESERVE.GROUP.LOCATION = 
  case_when(RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == TRUE  ~ "WJ MR",
            RESERVE.GROUP == "WJ" & RESERVE.GROUP.INSIDE == FALSE ~ "WJ CTRL",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == TRUE  ~ "FF MR",
            RESERVE.GROUP == "FF" & RESERVE.GROUP.INSIDE == FALSE ~ "FF CTRL",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == TRUE  ~ "LS MR",
            RESERVE.GROUP == "LS" & RESERVE.GROUP.INSIDE == FALSE ~ "LS CTRL")
            )

# stack data for subsequent analysis - check dimensions
dim(long_table) # 267 x 71

# stack data for subsequent analysis -correct type in previous data for succesful stacking
long_table %<>% mutate(DEPTH.M = as.numeric(DEPTH.M))

# stack data for subsequent analysis
long_table %<>% bind_rows(long_table, lt_obis_truncated)
dim(long_table) # 1828   71

# define BRUV.PRES and set EDNA.PRES and BOTH.PRES so as to identify complete BRUV/EDNA data
# define ALL.PRES so as to identify complete BRUV/EDNA/OBIS data 

long_table %<>% mutate( BRUV.PRES = case_when(SAMPLE.TYPE == "BRUV" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table %<>% mutate( EDNA.PRES = case_when(SAMPLE.TYPE == "eDNA" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table %<>% mutate( BOTH.PRES = case_when(BRUV.PRES == 1 | EDNA.PRES == 1 ~ 1, TRUE ~ 0))

long_table %<>% mutate( OBIS.PRES = case_when(SAMPLE.TYPE == "OBIS" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table %<>% mutate(  ALL.PRES = case_when(BRUV.PRES == 1 | EDNA.PRES == 1 | OBIS.PRES == 1 ~ 1, TRUE ~ 0))


# rearrange columns as in previous data combination
long_table %<>% relocate(SET.ID,	REP.ID, SAMPLE.TYPE, LOC.NAME, MH.GPS.LAT,
  MH.PPS.LONG, RESERVE.GROUP,  RESERVE.GROUP.INSIDE, RESERVE.GROUP.LOCATION, SUPERKINGDOM,	
  PHYLUM,	CLASS,	ORDER,	FAMILY,	GENUS,	SPECIES)

# VIII. Check data completness and citations
# ==========================================

# safe table with citation info
data_citataions <- lt_obis_results %>% ungroup() %>% select(bibliographicCitation) %>% filter(!is.na(bibliographicCitation)) %>% 
   distinct() %>% arrange(bibliographicCitation) %>% print(n = Inf) 

write.xlsx(data_citataions, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210707_OBIS_data_citations.xlsx", asTable = TRUE, overwrite = FALSE)


# check data completeness - preformatting
left_join(ungroup(select(lt_obis_results, id)), ungroup(select(long_table, ASV)), by = c("id" = "ASV"))
OBIS_records <- ungroup(lt_obis_results) %>% select(id) %>% distinct() %>% rename(OBIS_record = id)
OBIS_records %<>%  mutate(used = case_when(OBIS_record %in%  long_table$ASV ~ TRUE, TRUE ~ FALSE))

# and numerical summaries for manuscript
nrow(OBIS_records)
sum(OBIS_records$used)
sum(OBIS_records$used) / nrow(OBIS_records)

# XI. Save results and workspace
# ==============================

# for verbosity
write.xlsx(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/998_r_map_and_add_obis__full_data_raw.xlsx", asTable = FALSE)
# for superseded QGIS mapping in /Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210307_sample_map.qgz
write.csv(long_table, "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/998_r_map_and_add_obis__full_data_raw.csv")

# saving workspace manually
save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/210705_998_r_map_and_add_obis.Rdata")
save.image("/Users/paul/Documents/OU_eDNA/201028_Robjects/210705_998_r_map_and_add_obis.Rdata")

# for subsequent analyses
saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_map_and_add_obiss__full_data_raw.Rds")
saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/998_r_map_and_add_obiss__full_data_raw.Rds")

