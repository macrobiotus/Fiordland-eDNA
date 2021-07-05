#   **************************************
#   * Get a nice map and add OBIS data   *
#   *                                    *  
#   **************************************
#   26-Feb-2021, 1-Mar-2021

# I. Load packages
# ================
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse")   # tibbles, pipes, and more
library("magrittr") # more pipes

library("robis")       # access OBIS data

library("sf")          # create buffer around sampling areas
library("sp")


library("readxl")      # read Excel files
library("openxlsx")    # write Excel tables


# II. Functions
# =============

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)

# III. Load data 
# ==============

readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210301_997_r_format_longtables__analysis_input.Rds")

# VII. add OBIS data and get a map
# ================================
# use CRS 4326
# keep in mind 
#  long_table_dt_map <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("MH.GPS.LAT", "MH.PPS.LONG", "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]
# of
#  /Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r

# get clean raw data
# -------------------

# define BRUV.PRES and set EDNA.PRES and BOTH.PRES so as to identify complete BRUV/EDNA data 
long_table <- long_table %>% mutate( BRUV.PRES = case_when(SAMPLE.TYPE == "BRUV" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table <- long_table %>% mutate( EDNA.PRES = case_when(SAMPLE.TYPE == "eDNA" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table <- long_table %>% mutate( BOTH.PRES = case_when(BRUV.PRES == 1 | EDNA.PRES == 1 ~ 1, TRUE ~ 0))

# get a table with relevant columns  for OBIS lookup
lt_obis_lookup <- long_table %>% 
  select("SET.ID", "MH.GPS.LAT", "MH.PPS.LONG",  "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION") %>%
  distinct() %>% arrange(SET.ID) 
lt_obis_lookup %>% print(n= Inf)


# get clean spatial data
# ----------------------
#  
# https://gis.stackexchange.com/questions/292327/creating-buffers-around-points-and-merging-with-spatialpolygonsdataframe-to-crea
# https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/understand-epsg-wkt-and-other-crs-definition-file-types/
# check cordinates sytems
EPSG <- rgdal::make_EPSG()
EPSG %>% filter(code == 4326)
head(EPSG, 20)

# for buffer calculation - get sample coordinates as sf object with coordinates in EPSG 4326 / WGS84
lt_obis_lookup_sf <- st_as_sf(lt_obis_lookup,coords=c("MH.PPS.LONG","MH.GPS.LAT")) 
lt_obis_lookup_sf %<>% st_set_crs(4326) # set CRS to WGS84 

# for test maps  - get background layer coordinates as sf object with coordinates in EPSG 4326 / WGS84
nzshp_hires = read_sf("/Users/paul/GIS/NZ_coast/NZ_Coast_isl.shp")
nzshp_hires_WGS84 <- st_transform(nzshp_hires, crs = 4326)

# check sf objects - looking ok so far
ggplot() +
    geom_sf(data = nzshp_hires_WGS84, fill = "lightgrey") + 
    geom_sf(data = lt_obis_lookup_sf, colour = "red") + 
    coord_sf( xlim = c((166.5-0.1), (167.0+0.1)), ylim = c((-46.04-0.1),(-45.52+0.1)), expand = FALSE) +
    theme_bw()


# calculate buffers
# ------------------

# for buffer generation re-project objects to local kms and check again
lt_obis_lookup_sf_loc <- lt_obis_lookup_sf %<>% st_transform(crs = st_crs("+proj=utm +zone=58G +datum=WGS84 +units=km"))
nzshp_hires_WGS84_loc <- nzshp_hires_WGS84 %<>% st_transform(crs = st_crs("+proj=utm +zone=58G +datum=WGS84 +units=km"))

# calculate 1 km buffers 
lt_obis_lookup_sf_buffer <- st_buffer(lt_obis_lookup_sf, 1)

# check sf objects - still looking good so far - bounding box as defined per lt_obis_lookup_sf and 10 km in addition
#  buffers working as well
ggplot() +
    geom_sf(data = nzshp_hires_WGS84, fill = "lightgrey") + 
    geom_sf(data = lt_obis_lookup_sf_buffer, colour = "red") + 
    geom_sf(data = lt_obis_lookup_sf, colour = "red") + 
    coord_sf( xlim = c((619.6011-10), (653.8977+10)), ylim = c((-5100.241-10),(-5042.894+10)) , expand = FALSE) +
    theme_bw()

# back transform buffers to degree and wkt format for obis lookup
lt_obis_lookup_sf_buffer_d <- lt_obis_lookup_sf_buffer %>% st_transform(crs = 4326)

# add column with buffers in Well-Known-Text format for OBIS lookup
lt_obis_lookup_sf_buffer_d %<>% mutate(wkt_bu = st_as_text(geometry))

# fetch OBIS data
# ---------------

# quick and dirty - get data for AphiaID "2" = animalia
obis_data_raw <- lapply(lt_obis_lookup_sf_buffer_d$wkt_bu, function (x)  occurrence(taxon = 2, geometry = x))


# get data for AphiaID "2" = animalia - will get a sf dataframe with nested occurrence data 
lt_obis_lookup <-  lt_obis_lookup_sf_buffer_d %>% mutate(obis = occurrence(taxon = 141438, geometry = wkt_bu)) 












map_dfr()

get_obis_data = function(df_wkt){

print(wkt_bu)

}

lapply(lt_obis_lookup_sf_buffer_d, get_obis_data)


occurrence(taxonid = 2, geometry = foo, startdepth = 0)
?occurrence()


records <- occurrence(taxon = 141438, geometry = "POLYGON ((0 0, 0 45, 45 45, 45 0, 0 0))")

names(records)

# later filter as in eDNA script
# /Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_long_table.r
# psob_molten_clean_marine <- psob_molten_clean_chordates_top %>% filter(
#   CLASS %in% c("Actinopteri", "Chondrichthyes") |
#   ORDER %in% c("Cetacea") | 
#   FAMILY %in% c("Otariidae")) %>% filter(!(GENUS %in% c("Sardinops")))



?occurrence()



















