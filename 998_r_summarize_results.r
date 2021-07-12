#  **********************************
#  *                                *        
#  *  Get results from long tables  *
#  *                                *
#  **********************************


lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
rm(list = ls(all.names = TRUE))
gc()

# I. Load packages and define functions
# =====================================

library("tidyverse")   # because we can't stop using it anymore
library("magrittr")    # get the %<>% pipe

library("ggpubr") # combine plots -  http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/


library("sf")           # simple feature objects


# library("ggrepel")     # to improve plot labels
# 
# library("future.apply") # faster handling of large tables
# library("data.table")   # faster handling of large tables
# 

# library("rmapshaper")   # simplify shape file layers
# library("ggsflabel")    # label simple feature in ggplot  https://github.com/yutannihilation/ggsflabel - possibly inluded in ggplot
# 
# library("eulerr")       # to compare BRIUV and eDNA
# library("ggplotify")    # base R to Ggplot
# 
# library("vegan")        # for NMDS 
# library("indicspecies") # indicator species  - see citation below
# 
# library("FactoMineR") # MCA
# library("explor")     # check MCA results in browser
# library("factoextra") # get MCA results summaries
# 
# library("jpeg")   # read in jpeg images - see line ~840

# library("nVennR")
# library("UpSetR")    # Conway, J. R., Lex, A. & Gehlenborg, N. 2017 UpSetR: an R package for the
#                      # visualization of intersecting sets and their properties. Bioinformatics 33,
#                      # 2938?2940. (doi:10.1093/bioinformatics/btx364)
#                      # 
#                      # documentation at https://rdrr.io/cran/UpSetR/man/upset.html - hard to follow

# get Euler objects for plotting
# ------------------------------
get_euler_object = function(level, tibl){
  require("eulerr")
  require("tidyverse")
  require("magrittr")
  
  # check if needed columns are in the input data
  stopifnot(c("BRUV.OBS.PRES", "EDNA.OBS.PRES", "OBIS.OBS.PRES") %in% names(tibl))
  stopifnot(level %in% c("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"))
  
  # isolate realvant columns for summary
  tibl %<>% select(SET.ID, BRUV.OBS.PRES, EDNA.OBS.PRES, OBIS.OBS.PRES, RESERVE.GROUP, RESERVE.GROUP.LOCATION, SUPERKINGDOM,  PHYLUM,  CLASS,  ORDER,  FAMILY,  GENUS, SPECIES) %>% distinct()
  
  # sum up unique presences fr Euler plot
  tibl %<>% group_by(get(level)) %>% summarise(eDNA = as.logical(sum(EDNA.OBS.PRES)),
                                          BRUV = as.logical(sum(BRUV.OBS.PRES)),
                                          OBIS = as.logical(sum(OBIS.OBS.PRES))
                                          )
  return(euler(tibl[ , 2:4]))
}

# get Euler Ggplots 
# -----------------
get_euler_ggplot = function(level, euler_ob, plot_label = TRUE){
  require("tidyverse")
  require ("ggplotify")

  # sanitize input
  stopifnot( class(euler_ob)[1] == "euler")
  stopifnot(level %in% c("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"))
  
  euler_ggplot <- as.ggplot(
    plot(euler_ob, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels = list(font=1, cex=0.8))
    ) + {if(plot_label == TRUE) labs(subtitle = str_to_sentence(level))}
  
  return(euler_ggplot)
}

# get a table with relevant columns for mapping
# ---------------------------------------------
#   (from full_biodiv or fish_biodiv )
get_sf_biodiv =  function(tibl){
  require("tidyverse")
  require("magrittr")
  require("sf")
  
  # define columns for mapping add input verification
  cols <- c("SET.ID", "MH.GPS.LAT", "MH.PPS.LONG",  "RESERVE.GROUP", "RESERVE.GROUP.INSIDE",
            "RESERVE.GROUP.LOCATION", "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS",
            "SPECIES", "ASV", "ABUNDANCE", "SAMPLE.TYPE", "BRUV.OBS.PRES", "EDNA.OBS.PRES", "OBIS.OBS.PRES")
  stopifnot(cols %in% names(tibl))
  
  # select relavant data fro mapping
  tibl %<>% ungroup %>% select(all_of(cols)) %>% arrange(SET.ID) 
  
  # get simple feature df for mapping, define coordinates as WGS84 (degrees)  
  tibl %<>% st_as_sf(coords=c("MH.PPS.LONG","MH.GPS.LAT")) %>% st_set_crs(4326) 

}

# get bounding box around an area defined by a variable (here default: RESERVE.GROUP.LOCATION)  
# ------------------------------------------------------
# from https://stackoverflow.com/questions/54696440/create-polygons-representing-bounding-boxes-for-subgroups-using-sf

get_bbox_anyloc <- function(tibl, location = c("RESERVE.GROUP.LOCATION")){
  require("tidyverse")
  require("magrittr")
  require("sf")
  
  # sanitize input and 
  stopifnot( c(location, "MH.PPS.LONG", "MH.GPS.LAT") %in% names(tibl))
  
  # helper function 
  calc_angle <- function(lon,lat) {
    cent_lon <- mean(lon)
    cent_lat <- mean(lat)
    ang <- atan2(lat - cent_lat, lon - cent_lon)
    return(ang)
  }
  
  # calculate bounding box
  bbox <- tibl %>% group_by(across(location)) %>%
  summarise(xmin = min(MH.PPS.LONG) -0.01 ,ymin = min(MH.GPS.LAT) -0.01, xmax=max(MH.PPS.LONG) +0.01,  ymax = max(MH.GPS.LAT) +0.01) %>%
  gather(x,lon,c('xmin','xmax')) %>% gather(y,lat,c('ymin','ymax')) %>%
  st_as_sf(coords=c('lon','lat'),crs=4326,remove=F) %>%
  group_by(across(location)) %>% mutate(angle = calc_angle(lon,lat)) %>%
  arrange(angle) %>% summarise(do_union=FALSE) %>% st_cast('POLYGON')
  
  return(bbox)
}

# get hacked data frame for heat-map plotting 
# --------------------------------------------
#  - isolate coordinates into seperate columns 
#  - ass two crows matching mapping extend to extend plot
get_plot_df = function(sf_df, show_var = NULL) {
  require("sf")
  require("purrr")
  require("magrittr")
  stopifnot("sf" %in% class(sf_df))
   
  sf_df %<>% mutate(lat = unlist(map(sf_df$geometry,2)),
                    lon = unlist(map(sf_df$geometry,1))
                    ) %>% st_drop_geometry %>% 
                    {if(!is.null(show_var)) filter(.,SAMPLE.TYPE == show_var) else .} %>%
                    add_row(tibble_row(lon = 600, lat = -5200, SAMPLE.TYPE = "eDNA")) %>%
                    add_row(tibble_row(lon = 700, lat = -5000, SAMPLE.TYPE = "eDNA")) 
                    
  return(sf_df)
}

# aggregate discrete observation of either method ("BOTH.PRES") per sampling area (e.g.: "RESERVE.GROUP.LOCATION") on "GENUS" or species  level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
#   used to get distance matrices in vegan and for numerical summaries
# get_taxon_matrix <- function(long_dt = long_table_dt , group_var = "RESERVE.GROUP.LOCATION", level = "GENUS") {
# 
#   if (level == "GENUS") {
#     # aggregate dt for provided grouping variable
#     long_table_dt_agg_group_var_level <- long_dt[, lapply(.SD, sum, na.rm=TRUE), by=c(group_var, "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]
#   } else if (level == "SPECIES") {
#     # aggregate dt for provided grouping variable
#     long_table_dt_agg_group_var_level <- long_dt[, lapply(.SD, sum, na.rm=TRUE), by=c(group_var, "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"), .SDcols=c("BOTH.PRES") ]
#   } else {
#     stop("Level needs to be set to either \"GENUS\" or \"SPECIES\"")
#   }
#   
#   # cast matrix
#   taxon_matrix <- as.matrix(data.table::dcast(setDT(long_table_dt_agg_group_var_level), get(group_var)~get(level), value.var="BOTH.PRES", sum, fill=0), rownames=TRUE)
#   return(taxon_matrix)
# 
# }


# II. Read in data
# ================

# check input data of previous script
system("open -a \"Microsoft Excel\" \"/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/998_r_map_and_add_obis__full_data_raw.xlsx\"")

long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_map_and_add_obiss__full_data_raw.Rds")


# III. Read in and format data  
# ============================
# - mark non-NZ species  **(possibly needs to be re-worked)**
# - split "fish" and "full" data
# - filter for data completeness **(possibly needs to be re-worked)**

# Mark non-NZ species  **(possibly needs to be re-worked)**
# ---------------------------------------------------------
#  16-Mar-2021 add asterisks ("*") to non-NZ species, and ("**") to non-fish (mammals and crustaceans)
#  after checking with list 
#  Roberts, C., Stewart, A., Struthers, C., Barker, J. & Kortet, S. 2019 Checklist of the Fishes of New Zealand. 

nonnz_fish <- c("Asterropteryx", "Banjos", "Benitochromis", "Bostrychus", "Bovichtus", "Caprodon", "Coptodon", "Engraulis", "Gobiesox", "Gymnoscopelus", "Helcogramma", "Microcanthus", "Opistognathus", "Phoxinus", "Sander", "Scobinichthys")
nonnz_othr <- c("Macroctopus", "Jasus", "Arctocephalus", "Balaenoptera", "Tursiops")

long_table %<>% mutate(GENUS = 
                         case_when(GENUS %in% nonnz_fish ~ paste0(GENUS, "*"),
                                   GENUS %in% nonnz_othr ~ paste0(GENUS, "**"),
                                                    TRUE ~ GENUS)
                                                    )

# Filter for data completeness **(possibly needs to be re-worked)**
# ------------------------------------------------------------------

# - not done yet -


# continue her after 7-Jul-2021: 

# possibly get equivalent of BOTH.PRES
# ------------------------------
# function needs to 
#   get presence / absence on a {taxonomic level} (SPECIES)
#   per a {location} (SET.ID  RESERVE.GROUP.LOCATION )
#   ? check data completeness  {all "1" in BRUV.OBS.PRES	EDNA.OBS.PRES	OBIS.OBS.PRES}

# Split "fish" and "full" data
# ----------------------------
full_biodiv <- long_table %>% distinct()
fish_biodiv <- long_table %>% distinct() %>% filter(CLASS %in% c("Actinopteri", "Chondrichthyes")) %>% filter(!(GENUS %in% c("Sardinops")))

# III. Get Euler plots
# ====================

# get euler analysis results for plotting / plot_label = TRUE shrinks plots a lot
euler_obs_full_bio <- lapply(list("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"), get_euler_object, full_biodiv)
euler_ggp_full_bio <- mapply(get_euler_ggplot, list("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"),  euler_obs_full_bio, plot_label = FALSE, SIMPLIFY = FALSE)

# plot euler analysis results
euler_obs_fish_bio <- lapply(list("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"), get_euler_object, fish_biodiv)
euler_ggp_fish_bio <- mapply(get_euler_ggplot, list("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"),  euler_obs_fish_bio, plot_label = FALSE, SIMPLIFY = FALSE)

# create compound plot with better labels then with plot_label = TRUE above
ggarrange(
  ggarrange(plotlist = euler_ggp_full_bio,  ncol = 1, nrow = 6,
    labels = str_to_sentence(c("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES")),
    font.label = list(size = 12, color = "black", face = "bold.italic", family = NULL),
    vjust = 4.5
    ),
  ggarrange(plotlist = euler_ggp_fish_bio,  ncol = 1, nrow = 6), ncol = 2,
    labels = c("a","b")
  )

# save compound plot with better labels then with plot_label = TRUE above
ggsave("210712_998_r_summarize_results__euler_edna_bruv_obis.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1.5, width = 75, height = 175, units = c("mm"),
         dpi = 500, limitsize = TRUE)  

# IV. get geographical maps with heat overlays
# =============================================

# compare script  ~/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r

# data preparation
# ----------------

# for mapping: get column-subset sf's with WGS 84 in degrees 
full_biodiv_sf <- get_sf_biodiv(full_biodiv)
fish_biodiv_sf <- get_sf_biodiv(fish_biodiv)

# for mapping: get map layers 
nzshp_hires_WGS84_sf <- read_sf("/Users/paul/GIS/NZ_coast/NZ_Coast_isl.shp") %>% st_transform(crs = 4326)
nzshp_lores_WGS84_sf <- rmapshaper::ms_simplify(input = as(nzshp_hires_WGS84_sf, 'Spatial')) %>% st_as_sf

# for mapping: define bounding boxes as in map in previous script 
#  field work area & sample groups
bbox_fwork <- st_as_sfc(st_bbox(c(xmin = (166.5-0.1), xmax = (167.0+0.1), ymax = (-46.04-0.1), ymin = (-45.52+0.1)), crs = st_crs(4326)))
#  boxes around default value RESERVE.GROUP.LOCATION  
bbox_rgl_full_biodiv <- get_bbox_anyloc(full_biodiv) # must use original object, not sf 
bbox_rgl_fish_biodiv <- get_bbox_anyloc(fish_biodiv) # must use original object, not sf

# for mapping and buffer calculations at correct scale: re-project all sf's to local km  
get_reprojection <- function(sf) st_transform(sf, crs = st_crs("+proj=utm +zone=58G +datum=WGS84 +units=km"))

full_biodiv_sf_km <- get_reprojection(full_biodiv_sf)
fish_biodiv_sf_km <- get_reprojection(fish_biodiv_sf)

nzshp_hires_WGS84_sf_km <- get_reprojection(nzshp_hires_WGS84_sf)
nzshp_lores_WGS84_sf_km <- get_reprojection(nzshp_lores_WGS84_sf)

bbox_fwork_km <- get_reprojection(bbox_fwork)

bbox_rgl_full_biodiv_km <- get_reprojection(bbox_rgl_full_biodiv)
bbox_rgl_fish_biodiv_km <- get_reprojection(bbox_rgl_fish_biodiv)

# calculate 2.5 km buffers
full_biodiv_sf_km_sid_buff <- full_biodiv_sf_km %>% select("SET.ID") %>% distinct %>% st_buffer(2.5)
fish_biodiv_sf_km_sid_buff <- fish_biodiv_sf_km %>% select("SET.ID") %>% distinct %>% st_buffer(2.5)

# get dataframes suitable for plotting with below functions - write as function
full_biodiv_df_edna <- get_plot_df(full_biodiv_sf_km, "eDNA")
full_biodiv_df_bruv <- get_plot_df(full_biodiv_sf_km, "BRUV")
full_biodiv_df_obis <- get_plot_df(full_biodiv_sf_km, "OBIS")

fish_biodiv_df_edna <- get_plot_df(fish_biodiv_sf_km, "eDNA")
fish_biodiv_df_edna <- get_plot_df(fish_biodiv_sf_km, "BRUV")
fish_biodiv_df_edna <- get_plot_df(fish_biodiv_sf_km, "OBIS")

get_ggeom_density(get_plot_df(fish_biodiv_sf_km))
get_ggeom_density(get_plot_df(fish_biodiv_sf_km, c("eDNA")))
get_ggeom_density(get_plot_df(fish_biodiv_sf_km, c("BRUV")))
get_ggeom_density(get_plot_df(fish_biodiv_sf_km, c("OBIS")))


# mapping
# --------

# inset map for subsequent main map (see https://geocompr.github.io/post/2019/ggplot2-inset-maps/)
#  using map in degrees, as scale may be too large for kms
map_inset <-  ggplot(data = nzshp_lores_WGS84_sf) + geom_sf(fill = "grey93", color = "red", lwd = 0.5) +
    geom_sf(data = bbox_fwork, fill = NA, color = "darkred", size = 1) + theme_void()

plot_full_biodiv <- ggplot() +
      geom_density_2d_filled(data = get_plot_df(full_biodiv_sf_km), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
      facet_grid(. ~ SAMPLE.TYPE) +
      geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
      # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") + 
      stat_sf_coordinates(data = full_biodiv_sf_km, aes(shape = RESERVE.GROUP), color = "grey20", size = 2) +
      stat_sf_coordinates(data = full_biodiv_sf_km, aes(shape = RESERVE.GROUP), color = "white", size = 1) +
      coord_sf(xlim = c((619.6011-10), (653.8977+10)), ylim = c((-5100.241-10),(-5042.894+10)) , expand = FALSE) +
      theme_bw() + 
      theme(legend.position= "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            )

plot_fish_biodiv <- ggplot() +
      geom_density_2d_filled(data = get_plot_df(fish_biodiv_sf_km), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
      facet_grid(. ~ SAMPLE.TYPE) +
      geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
      # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") + 
      stat_sf_coordinates(data = fish_biodiv_sf_km, aes(shape = RESERVE.GROUP), color = "grey20", size = 2) +
      stat_sf_coordinates(data = fish_biodiv_sf_km, aes(shape = RESERVE.GROUP), color = "white", size = 1) +
      coord_sf(xlim = c((619.6011-10), (653.8977+10)), ylim = c((-5100.241-10),(-5042.894+10)) , expand = FALSE) +
      theme_bw() + 
      theme(legend.position= "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            )

ggarrange( plot_full_biodiv, plot_fish_biodiv,  
  ncol = 1, nrow = 2, labels = c("a","b") )

# save compound plot with better labels then with plot_label = TRUE above
ggsave("210712_998_r_summarize_results__geoheat_edna_bruv_obis.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 2, width = 85, height = 85, units = c("mm"),
         dpi = 500, limitsize = TRUE)  

# V. Get a table
# ==============

# --- just drafting ----


full_biodiv  %>% select("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES") %>% distinct()
fish_biodiv  %>% select("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES") %>% distinct()



full_biodiv_sf_km %>% filter(SAMPLE.TYPE == "OBIS"") 






# --- unrevised code below ----






# IV. get barplots
# ================
# - not done yet -

get_default_ocurrence_plot = function (psob_molten, taxlev, taxlev_fill = taxlev, ptitl, pxlab, pylab, facet_var = "LOC.NAME" ){
  
  require(ggplot2)
  require(data.table)

  # melting Phyloseq object to data table for merging and speed
  psob_molten <- data.table(psob_molten)
 
  # set sorting key properly
  setkey(psob_molten,ASV)
 
  
  # aggregate on chosen level
  #   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
  psob_molten <- psob_molten[, lapply(.SD, sum, na.rm=TRUE), by=c(facet_var, "ASV", "SAMPLE.TYPE", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS",  "SPECIES"), .SDcols=c("ABUNDANCE") ]

  #  resort for clarity
  keycol <-c("ASV",facet_var)
  setorderv(psob_molten, keycol)
  
  # debugging
  # print(head(psob_molten))

  # add presence-absence abundance column
  psob_molten <- psob_molten[ , ASVPRESENT :=  fifelse(ABUNDANCE == 0 , 0, 1, na=NA)]
  
  # debugging
  # print(head(psob_molten))

  ggplot(psob_molten, aes_string(x = taxlev, y = "ASVPRESENT", fill = taxlev_fill)) +
    geom_bar(stat = "identity", position = "stack", colour = NA, size=0) +
    facet_grid(get(facet_var) ~ SAMPLE.TYPE, shrink = TRUE, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    theme(strip.text.y = element_text(angle=0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7), 
          axis.ticks.y = element_blank()) +
    labs( title = ptitl) + xlab(pxlab) + ylab(pylab)
  
}

get_default_ocurrence_plot (fish_biodiv, facet_var = "RESERVE.GROUP.LOCATION", taxlev = "FAMILY", ptitl = "Phyla across all sample types before filtering", pxlab = "phyla (NCBI taxonomy) ", pylab =  "ASV counts in sample category (y scales fixed)")



foo <- fish_biodiv %>% group_by(RESERVE.GROUP.LOCATION, SAMPLE.TYPE, SPECIES) %>% summarise(SUMABU = case_when(ABUNDANCE >= 1  ~ 1))


foo %>% print(n = Inf)

ggplot(foo, aes(fill=SAMPLE.TYPE, y=SUMABU, x=SPECIES)) + 
    geom_bar(position="stack", stat="identity") +
    facet_grid(rows = vars(RESERVE.GROUP.LOCATION), scales="free") + 
    theme_bw() +
    theme(strip.text.y = element_text(angle=0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7), 
          axis.ticks.y = element_blank()) +
    labs( title = "ptitl") + xlab("pxlab") + ylab("pylab")






















# VI. Present raw observations 
# =============================

# 1.) Get a simple barplot
# ------------------------

# aggregate discrete observation of wither method ("BOTH.PRES") per sampling area (RESERVE.GROUP.LOCATION) on GENUS level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
long_table_dt_agg_gen <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]

RESERVE.GROUP.LOCATION.LABS <- list(
  "LS MR"   = "LS MR\n (n = 4)",
  "LS CTRL" = "LS CTRL\n (n = 4)",
  "FF MR"   = "FF MR\n (n = 2)",
  "FF CTRL" = "FF CTRL\n (n = 3)",
  "WJ MR"   = "WJ MR\n (n = 4)",
  "WJ CTRL" = "WJ CTRL\n (n = 4)"
)

get_label <- function(variable,value){
  return(RESERVE.GROUP.LOCATION.LABS[value])
}

p_barobs <- ggplot(long_table_dt_agg_gen, aes_string(x =  "BOTH.PRES", y = reorder(long_table_dt_agg_gen$GENUS, desc(long_table_dt_agg_gen$GENUS)), fill = "BOTH.PRES")) +
    geom_bar(stat = "identity", position = "stack", colour = NA, size=0) +
    scale_fill_gradient(low="black", high="red") +
    facet_grid(.~RESERVE.GROUP.LOCATION, shrink = TRUE, scales = "fixed", labeller=get_label) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(strip.text.y = element_text(angle=0)) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7, face = "italic"), 
          axis.ticks.y = element_blank()) +
    labs(title = "Sampling Locations") + xlab("Unique Observations") + ylab("Genus")

ggsave("210407_998_r_summarize_results_observations.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 120, height = 150, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# 2.) Get numerical summaries for main text
# -----------------------------------------

# summary plain numbers
# ---------------------
# summary of observations (added 21.04.2021)

nrow(long_table_dt) # number of all obesrvations

long_table_dt %>% select(EDNA.PRES, SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS) %>% 
  filter(EDNA.PRES ==1)
  
long_table_dt %>% select(BRUV.PRES, SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS) %>% 
  filter(BRUV.PRES ==1) 
  

long_table_dt %>% select(BRUV.PRES, SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS) %>% 
  filter(BRUV.PRES ==1) %>% distinct()

long_table_dt %>% select(EDNA.PRES, SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS) %>% 
  filter(EDNA.PRES ==1) %>% distinct()
  
long_table_dt %>% select(EDNA.PRES, SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS) %>% print(n = Inf)


# summary of BRUV depth
summary(long_table_dt$MH.BRUV.DEPTH[!is.na(long_table_dt$MH.BRUV.DEPTH)])
sd(long_table_dt$MH.BRUV.DEPTH[!is.na(long_table_dt$MH.BRUV.DEPTH)])

# summary of eDNA sampling depth 
summary(as.numeric(long_table_dt$DEPTH.M[!is.na(long_table_dt$DEPTH.M)]))
sd(as.numeric(long_table_dt$DEPTH.M[!is.na(long_table_dt$DEPTH.M)]))

# summery of genus observations
mat_rgl_spc <- get_taxon_matrix(long_table_dt, "RESERVE.GROUP.LOCATION", "SPECIES")
mat_rgl_gen <- get_taxon_matrix(long_table_dt, "RESERVE.GROUP.LOCATION", "GENUS")

summary(colSums(mat_rgl_spc))
summary(colSums(mat_rgl_gen))

# genus names with maximum observations in each location 
apply(mat_rgl_gen, 1, which.max)
colnames(mat_rgl_gen)[apply(mat_rgl_gen, 1, which.max)]

# species names with maximum observations in each location 
apply(mat_rgl_spc, 1, which.max)
colnames(mat_rgl_spc)[apply(mat_rgl_spc, 1, which.max)]

# most frequent genus observation
which(mat_rgl_gen == max(mat_rgl_gen), arr.ind = TRUE)
colnames(mat_rgl_gen)[which(mat_rgl_gen == max(mat_rgl_gen), arr.ind = TRUE)[2]]
mat_rgl_gen[which(mat_rgl_gen == max(mat_rgl_gen), arr.ind = TRUE)]

# most frequent species observation: what, where, how many
which(mat_rgl_spc == max(mat_rgl_spc), arr.ind = TRUE)
colnames(mat_rgl_spc)[which(mat_rgl_spc == max(mat_rgl_spc), arr.ind = TRUE)[2]]
mat_rgl_spc[which(mat_rgl_spc == max(mat_rgl_spc), arr.ind = TRUE)]

# heat maps
# ---------

# base heat maps 
#   https://cran.r-project.org/web/packages/plot.matrix/vignettes/plot.matrix.html
#   https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet.html

library(plot.matrix)
par(mar=c(5.1, 10, 4.1, 4.1))
plot(t(mat_rgl_gen), axis.col=list(side=1, las=1), axis.row = list(side=2, las=1), ann = FALSE, digits = 1, fmt.cell='%.0f')


# ggplot heat map with margin totals
#   as per check https://stackoverflow.com/questions/55787412/adding-marginal-totals-to-ggplot-heatmap-in-r

# get margin sums 
h_total <- long_table_dt_agg_gen %>% 
  group_by(GENUS) %>% 
  summarise(BOTH.PRES = sum(BOTH.PRES)) %>% 
  mutate(RESERVE.GROUP.LOCATION = 'TOTAL')

v_total <- long_table_dt_agg_gen %>% 
  group_by(RESERVE.GROUP.LOCATION) %>% 
  summarise(BOTH.PRES = sum(BOTH.PRES)) %>% 
  mutate(GENUS = 'TOTAL')
  
# add margin sums column for merging
v_total <- v_total %>% mutate(GENUS.SUM = BOTH.PRES) %>% select(RESERVE.GROUP.LOCATION, GENUS.SUM)
h_total <- h_total %>% mutate(RESERVE.GROUP.LOCATION.SUM = BOTH.PRES) %>% select(GENUS, RESERVE.GROUP.LOCATION.SUM)

long_table_dt_agg_gen_plot <- long_table_dt_agg_gen
long_table_dt_agg_gen_plot <- left_join(long_table_dt_agg_gen_plot, v_total, by = c("RESERVE.GROUP.LOCATION"), keep = FALSE  )
long_table_dt_agg_gen_plot <- left_join(long_table_dt_agg_gen_plot, h_total, by = c("GENUS"), keep = FALSE )

long_table_dt_agg_gen_plot <- long_table_dt_agg_gen_plot %>% mutate( RESERVE.GROUP.LOCATION = paste0(RESERVE.GROUP.LOCATION, " (", GENUS.SUM, ")"))
long_table_dt_agg_gen_plot <- long_table_dt_agg_gen_plot %>% mutate( GENUS = paste0(GENUS, " (", RESERVE.GROUP.LOCATION.SUM, ")"))


# expand margin to to be compatible with long dt used below - doesn't work 
# 
# h_total_long <- as_tibble(long_table_dt_agg_gen$GENUS) %>% arrange(value) %>% rename(GENUS = value)
# h_total_long <- left_join(h_total_long, h_total, by = "GENUS")
# 
# v_total_long <- as_tibble(long_table_dt_agg_gen$RESERVE.GROUP.LOCATION) %>% rename(RESERVE.GROUP.LOCATION = value)
# v_total_long <- left_join(v_total_long, v_total, by = "RESERVE.GROUP.LOCATION")
# get heat map

# add margin sums to new version of long data frame - looks crap
# long_table_dt_agg_gen_plot <- bind_rows(long_table_dt_agg_gen, h_total) 
# long_table_dt_agg_gen_plot <- bind_rows(long_table_dt_agg_gen_plot, v_total) 

p_heatobs <- ggplot(long_table_dt_agg_gen_plot, aes_string(x = "RESERVE.GROUP.LOCATION", y = reorder(long_table_dt_agg_gen_plot$GENUS, desc(long_table_dt_agg_gen_plot$GENUS)))) +
    geom_tile(aes(fill = BOTH.PRES) ) +
    scale_fill_gradient(low="black", high="red") +
    geom_text(size = 2, aes(label = BOTH.PRES, color = "white")) +
    theme_bw() +
    theme(legend.position = "none", 
          strip.text.y = element_text(angle=0), 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7, face = "italic"), 
          axis.ticks.y = element_blank()
          ) +
    xlab("Sampling Locations") + ylab("Genus Observations")

ggsave("210407_998_r_summarize_results_observations_heat.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 75, height = 165, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# summary corrected for sampling effort
# -------------------------------------

# ...(not done yet)...


# numerical summaries from margin totals
# --------------------------------------

# combined Genus observations
h_total <- long_table_dt_agg_gen %>% 
  group_by(GENUS) %>% 
  summarise(BOTH.PRES = sum(BOTH.PRES)) %>% 
  mutate(RESERVE.GROUP.LOCATION = 'TOTAL')

summary(h_total$BOTH.PRES)
# summary(h_total$BOTH.PRES)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.000   1.000   2.000   4.377   5.000  31.000 


h_total %>% filter(BOTH.PRES == max(BOTH.PRES))
# Parapercis - sand perches


# combined location observations
v_total <- long_table_dt_agg_gen %>% 
  group_by(RESERVE.GROUP.LOCATION) %>% 
  summarise(BOTH.PRES = sum(BOTH.PRES)) %>% 
  mutate(GENUS = 'TOTAL')

summary(v_total$BOTH.PRES)

#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    22.0    39.5    51.0    44.5    52.0    55.0 

v_total %>% filter(BOTH.PRES == max(BOTH.PRES))
# WJ MR                         55 TOTAL



# describe ASV  yield per primer
smpl_eff <- long_table %>% ungroup() %>% select(SET.ID, REP.ID, SAMPLE.TYPE, PRIMER.LABEL, RESERVE.GROUP.LOCATION) %>% filter(SAMPLE.TYPE == "eDNA") %>% arrange(SET.ID, RESERVE.GROUP.LOCATION) %>% print(n = Inf)
smpl_eff_grp <- smpl_eff %>% mutate(PRIMER.LABEL  = gsub(".*Mi", "", PRIMER.LABEL)) %>% group_by(RESERVE.GROUP.LOCATION, PRIMER.LABEL) %>% summarise(n = n()) %>% ungroup()
smpl_eff_grp %>% arrange(PRIMER.LABEL, RESERVE.GROUP.LOCATION)


# **** revisit this above - need species observations and genus observations *****
# sum genus observations for each factor
dt_genussum_rgl <- long_table_dt[,.(RESERVE.GROUP.LOCATION.GENUS.SUM=sum(BOTH.PRES)),.(RESERVE.GROUP.LOCATION)]
dt_genussum_rg  <- long_table_dt[,.(RESERVE.GROUP.GENUS.SUM=sum(BOTH.PRES)),.(RESERVE.GROUP)]

# correct genus observation effort for unequal sample effort
dt_genussum_rgl$RESERVE.GROUP.LOCATION.GENUS.SUM.PS <- dt_genussum_rgl$RESERVE.GROUP.LOCATION.GENUS.SUM / c(4,4,2,3,4,4)
dt_genussum_rgl
summary(dt_genussum_rgl$RESERVE.GROUP.LOCATION.GENUS.SUM.PS)

dt_genussum_rg$RESERVE.GROUP.LOCATION.SUM.PS <- dt_genussum_rg$RESERVE.GROUP.GENUS.SUM / c(8,5,8)
dt_genussum_rg
summary(dt_genussum_rg$RESERVE.GROUP.LOCATION.SUM.PS)



# VII. Show RESERVE.GROUP.LOCATION similarity based on GENUS overlap 
# ===================================================================
 
# reshape to observation matrix digestible by Vegan, discrete observations will be summed per genus
long_table_dt_agg_gen_mat <- get_taxon_matrix(long_table_dt, "RESERVE.GROUP.LOCATION", "GENUS")

jacc_matrix <- vegdist(long_table_dt_agg_gen_mat, distance="jaccard" )
summary(jacc_matrix)

# get a Jaccard distance matrix (distance define by overlap between sites)
#   see https://rpubs.com/CPEL/NMDS
#   see https://peat-clark.github.io/BIO381/veganTutorial.html
#   see https://fromthebottomoftheheap.net/2013/01/12/decluttering-ordination-plots-in-vegan-part-1-ordilabel/
#   see https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo


# run metaMDS - at pressence does just use the presence absence matrix (to keep genus scrores), but distance matrix is possible as well
long_table_dt_agg_gen_mat_jacc_NMS <-  metaMDS(long_table_dt_agg_gen_mat, distance="jaccard",  noshare = FALSE, k = 2, maxit = 5000,  trymax = 5000, wascores = TRUE, shrink = FALSE)
long_table_dt_agg_gen_mat_jacc_NMS
stressplot(long_table_dt_agg_gen_mat_jacc_NMS)

# basic plots - variant A
ordiplot(long_table_dt_agg_gen_mat_jacc_NMS, type = "none") 
orditorp(long_table_dt_agg_gen_mat_jacc_NMS, display =  "sites", cex = 1.25, air = 0.01)

# improve plot as shown here
#  https://jkzorz.github.io/2019/06/06/NMDS.html

#extract NMDS scores (x and y coordinates)
long_table_dt_agg_gen_mat_jacc_NMS.scores <- as_tibble(scores(long_table_dt_agg_gen_mat_jacc_NMS), rownames = "RESERVE.GROUP.LOCATION")
long_table_dt_agg_gen_mat_jacc_NMS.genus <- as_tibble(scores(long_table_dt_agg_gen_mat_jacc_NMS, "species"), rownames = "GENUS")


p_nmds <- ggplot(long_table_dt_agg_gen_mat_jacc_NMS.scores, aes(x = NMDS1, y = NMDS2)) +
   geom_point(data=long_table_dt_agg_gen_mat_jacc_NMS.genus, aes(x=NMDS1,y=NMDS2), size=3) +
   geom_point(size = 6, colour = "darkred", shape = c(16,16,17,17,15,15)) +
   geom_point(size = 5, colour = "red", shape = c(16,16,17,17,15,15)) +
   geom_label_repel(aes(label=RESERVE.GROUP.LOCATION), point.padding = 0.5) +
   coord_flip(xlim =c(-0.5, +0.5), ylim = c(-0.5, +0.5)) +
   theme_bw()
   
ggsave("210312_998_r_summarize_results_jaccard.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 75, height = 75, units = c("mm"),
         dpi = 500, limitsize = TRUE)     



# VIII. ANOSIM to test wether or not genus composition based on factors are significantly different
# ==========================================================================================

# On Anosim: 
# 1. CLARKE, K. R. 1993 Non-parametric multivariate analyses of changes in
# community structure. Austral Ecol. 18, 117?143.
# (doi:10.1111/j.1442-9993.1993.tb00438.x)
# "The ANalysis Of SIMilarity (ANOSIM) test has some similarity to an ANOVA-like
# hypothesis test, however, it is used to evaluate a dissimilarity matrix rather
# than raw data (Clarke, 1993). Further, raw (dis)similarities are often ranked 
# prior to performing an ANOSIM."
# "...the higher the R value, the more dissimilar [...] groups are in terms of [...] community composition."

# aggregate discrete observation of wither method ("BOTH.PRES") per sampling area (RESERVE.GROUP.LOCATION) on GENUS level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
long_table_dt_agg_gen_sets <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("SET.ID", "RESERVE.GROUP.LOCATION", "INSIDE.RESERVE", "RESERVE.GROUP", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]

# rename SET.ID to circumvent naming snafu with package data.table
setnames(long_table_dt_agg_gen_sets, "SET.ID", "SET_ID")

# reshape to observation matrix digestible by Vegan, discrete observations will be summed per genus
long_table_dt_agg_gen_mat_sets <- as.matrix(data.table::dcast(setDT(long_table_dt_agg_gen_sets), SET_ID~GENUS, value.var="BOTH.PRES", sum, fill=0), rownames=TRUE)


# A. Test if RESERVE.GROUP.LOCATIONs are significantly different 
# -----------------------------------------------------------------
# https://jkzorz.github.io/2019/06/11/ANOSIM-test.html
# - To test if there is a statistical difference between the fish communities of two or more groups of samples.
# - Null Hypothesis: there is no difference between the microbial communities of your groups of samples.


# get grouping variable of  RESERVE.GROUP.LOCATION
groupings <- as_tibble(long_table_dt_agg_gen_sets %>% select(SET_ID, RESERVE.GROUP.LOCATION)) %>% distinct()

long_table_dt_agg_gen_mat_sets_ano <-  anosim(long_table_dt_agg_gen_mat_sets, groupings$RESERVE.GROUP.LOCATION, distance = "jaccard", permutations = 9999)
summary(long_table_dt_agg_gen_mat_sets_ano)

# Dissimilarity: jaccard 
# 
# ANOSIM statistic R: 0.1992 
#       Significance: 0.0256 
# 
# Permutation: free
# Number of permutations: 9999

# B. Test if sample locations (RESERVE.GROUP) are significantly different 
# -----------------------------------------------------------------

# get grouping variable of  RESERVE.GROUP.LOCATION
groupings <- as_tibble(long_table_dt_agg_gen_sets %>% select(SET_ID, RESERVE.GROUP)) %>% distinct()

long_table_dt_agg_gen_mat_sets_ano <-  anosim(long_table_dt_agg_gen_mat_sets, groupings$RESERVE.GROUP, distance = "jaccard", permutations = 9999)
summary(long_table_dt_agg_gen_mat_sets_ano)

# Call:
# anosim(x = long_table_dt_agg_gen_mat_sets, grouping = groupings$RESERVE.GROUP,      permutations = 9999, distance = "jaccard") 
# Dissimilarity: jaccard 
# 
# ANOSIM statistic R: 0.1045 
#       Significance: 0.0956 
# 
# Permutation: free
# Number of permutations: 9999


# C. Test if inside/outside MR (INSIDE.RESERVE) are significantly different 
# -----------------------------------------------------------------

# get grouping variable of  RESERVE.GROUP.LOCATION
groupings <- as_tibble(long_table_dt_agg_gen_sets %>% select(SET_ID, INSIDE.RESERVE)) %>% distinct()

long_table_dt_agg_gen_mat_sets_ano <-  anosim(long_table_dt_agg_gen_mat_sets, groupings$INSIDE.RESERVE, distance = "jaccard", permutations = 9999)

# Call:
# anosim(x = long_table_dt_agg_gen_mat_sets, grouping = groupings$INSIDE.RESERVE,      permutations = 9999, distance = "jaccard") 
# Dissimilarity: jaccard 
# 
# ANOSIM statistic R: -0.06009 
#       Significance: 0.829 
# 
# Permutation: free
# Number of permutations: 9999


# IX. Show indicator genera for inside/outside each/all reserve(s)
# ============================================================================
# indicator species analysis
# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

# Using package indicspecies
# De C?ceres, M., Legendre, P. & Moretti, M. 2010 Improving indicator species
# analysis by combining groups of sites. Oikos 119, 1674?1684.
# (doi:10.1111/j.1600-0706.2010.18334.x)

long_table_dt_agg_gen_sets <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("SET.ID", "INSIDE.RESERVE", "RESERVE.GROUP.LOCATION", "RESERVE.GROUP", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]

# rename SET.ID to circumvent naming snafu with package data.table
setnames(long_table_dt_agg_gen_sets, "SET.ID", "SET_ID")

# reshape to observation matrix digestible by Vegan, discrete observations will be summed per genus
long_table_dt_agg_gen_mat_sets <- as.matrix(data.table::dcast(setDT(long_table_dt_agg_gen_sets), SET_ID~GENUS, value.var="BOTH.PRES", sum, fill=0), rownames=TRUE)

# define grouping vectors 
group.INSIDE.RESERVE <- as_tibble(long_table_dt_agg_gen_sets %>% select( SET_ID, INSIDE.RESERVE)) %>% distinct() %>% pull(INSIDE.RESERVE)
group.RESERVE.GROUP <- as_tibble(long_table_dt_agg_gen_sets %>% select(SET_ID, RESERVE.GROUP)) %>% distinct() %>% pull(RESERVE.GROUP)
group.RESERVE.GROUP.LOCATION <- as_tibble(long_table_dt_agg_gen_sets %>% select(SET_ID, RESERVE.GROUP.LOCATION)) %>% distinct() %>% pull(RESERVE.GROUP.LOCATION)

# A. Find indicator species INSIDE.RESERVE
# ----------------------------------------

ind_ir = multipatt(long_table_dt_agg_gen_mat_sets, group.INSIDE.RESERVE, func = "r.g", control = how(nperm=9999))
summary(ind_ir)


# Multilevel pattern analysis
#  ---------------------------
# 
#  Association function: r.g
#  Significance level (alpha): 0.05
# 
#  Total number of species: 61
#  Selected number of species: 1 
#  Number of species associated to 1 group: 1 
# 
#  List of species associated to each combination: 
# 
#  Group TRUE  #sps.  1 
#                 stat p.value  
# Opistognathus* 0.442  0.0365 *


# B. Find indicator species at each RESERVE.GROUP
# -----------------------------------------------

ind_rgl = multipatt(long_table_dt_agg_gen_mat_sets, group.RESERVE.GROUP, func = "r.g", control = how(nperm=9999))
summary(ind_rgl)

# Multilevel pattern analysis
#  ---------------------------
# 
#  Association function: r.g
#  Significance level (alpha): 0.05
# 
#  Total number of species: 61
#  Selected number of species: 2 
#  Number of species associated to 1 group: 2 
#  Number of species associated to 2 groups: 0 
# 
#  List of species associated to each combination: 
# 
#  Group FF  #sps.  1 
#              stat p.value  
# Galeorhinus 0.555  0.0433 *
# 
#  Group WJ  #sps.  1 
#               stat p.value  
# Forsterygion 0.632  0.0239 *
# ---
# Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1


# C. Find indicator species at each RESERVE.GROUP.LOCATION
# ---------------------------------------------------------

ind_rgl = multipatt(long_table_dt_agg_gen_mat_sets, group.RESERVE.GROUP.LOCATION, func = "r.g", control = how(nperm=9999))
summary(ind_rgl)

# Multilevel pattern analysis
#  ---------------------------
# 
#  Association function: r.g
#  Significance level (alpha): 0.05
# 
#  Total number of species: 61
#  Selected number of species: 3 
#  Number of species associated to 1 group: 2 
#  Number of species associated to 2 groups: 1 
#  Number of species associated to 3 groups: 0 
#  Number of species associated to 4 groups: 0 
#  Number of species associated to 5 groups: 0 
# 
#  List of species associated to each combination: 
# 
#  Group FF MR  #sps.  1 
#                 stat p.value   
# Opistognathus* 0.934  0.0036 **
# 
#  Group LS MR  #sps.  1 
#          stat p.value  
# Banjos* 0.775  0.0213 *
# 
#  Group FF MR+LS CTRL  #sps.  1 
#                 stat p.value  
# Scobinichthys* 0.676  0.0462 *


# X. Multiple Correspondence analysis
# =======================================

# as per https://rpubs.com/gaston/MCA

# select pertinent variables for further analysis
# long_table_dt_sfct <- long_table_dt[, c("SAMPLE.TYPE", "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION", "INSIDE.RESERVE", "BRUV.PRES", "EDNA.PRES", "CLASS", "ORDER", "FAMILY", "GENUS")] %>% unique()
long_table_dt_sfct <- long_table_dt[, c("RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION", "GENUS")] %>% unique()

# get genus names in Italics for later plotting 
# long_table_dt_sfct$GENUS <-  paste0("italic(", long_table_dt_sfct$GENUS,")")

# convert all variables to factors
long_table_dt_sfct <- long_table_dt_sfct[, lapply(.SD, as.factor)]

# number of categories per variable
cats <- apply(long_table_dt_sfct, 2, function(x) nlevels(as.factor(x)))
cats

# apply MCA
mca1 = MCA(long_table_dt_sfct, graph = FALSE)
summary(mca1)
# explor(mca1)


# table of eigenvalues
mca1$eig

dim_1_perc <-  signif(mca1$eig[1,3], digits = 2)  
dim_2_perc <- signif(mca1$eig[2,3] - dim_1_perc, digits = 2)  


plot(mca1$eig[ , 2])

# column coordinates
head(mca1$var$coord)

# row coordinates
head(mca1$ind$coord)

# data frames for ggplot
mca1_vars_df = data.frame(mca1$var$coord, Variable = rep(names(cats), cats))
mca1_obs_df = data.frame(mca1$ind$coord)

     

# see https://rpkgs.datanovia.com/factoextra/reference/fviz_contrib.html
p_cntrb <- fviz_contrib(mca1, choice="var", axes = 1, top = 10, fill = "lightgray", color = "grey") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(title = element_blank()) +
  theme(plot.background = element_rect(colour = "black"))
ggsave("210312_998_r_summarize_results_mca_dim1.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1.3, width = 50, height = 25, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# MCA plot of observations and categories
p_mca <- ggplot(data = mca1_obs_df, aes(x = Dim.1, y = Dim.2)) + 
    geom_hline(yintercept = 0, colour = "gray70") + 
    geom_vline(xintercept = 0, colour = "gray70") + 
    geom_point(colour = "gray50", alpha = 0.7) + 
    geom_density2d(colour = "gray80") +
    scale_colour_discrete(name = "Variable") +
    geom_label_repel(data = mca1_vars_df, aes(x = Dim.1, y = Dim.2, label = rownames(mca1_vars_df), colour = Variable), max.overlaps = Inf, point.size = NA) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab(paste0("Dim. 1 (", dim_1_perc, "% Variance)")) + 
    ylab(paste0("Dim. 2 (", dim_2_perc, "% Variance)")) +
    annotation_custom(ggplotGrob(p_cntrb), xmin = -1.3, xmax = -0.4, ymin = 1.2, ymax = 2.3)
    
ggsave("210408_998_r_summarize_results_mca.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 250, height = 250, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# XI. Combine plots for manuscript 
# ================================

# 7-Apr-2021: needs revisions once imagery is received and completed above

# 1.) plot arrangement 1: Venn - NMDS - MCA
# -----------------------------------------
ggarrange(
  ggarrange(pPhl, pCls, pOrd, pFam, pGen, pSpc, ncol = 1, nrow = 6, labels = c("(a)")),
  ggarrange(
    ggarrange(p_nmds, p_cntrb,  ncol = 2, nrow = 1, labels = c("(b)", "(d)")),
    ggarrange(p_mca, ncol = 1, nrow = 1, labels = c("(c)")), 
    nrow = 2, heights = c(3, 7)
    ), ncol = 2, widths = c(2, 8)
  )    

ggsave("210312_998_r_summarize_results_fig2_draft_Venn.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 250, height = 300, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# 2.) plot arrangement 2: Venn - barplot - NMDS
# ---------------------------------------------
ggarrange(
  ggarrange(pPhl, pCls, pOrd, pFam, pGen, pSpc, ncol = 1, nrow = 6, labels = c("(a)")),
  ggarrange(p_barobs, ncol = 1, nrow = 1, labels = c("(b)")),
  ggarrange(p_nmds, ncol = 1, nrow = 1, labels = c("(c)")),
  ncol = 3, nrow = 1, widths = c(1, 2, 2))

ggsave("210407_998_r_summarize_results_fig2_draft_Venn.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 295, height = 150, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# 3.) plot arrangement 3: Venn - heatmap - NMDS - imagery
# -------------------------------------------------------

# read in imagery - now public domain (see REDME for credits) - later original artwork 
# img_bw <- readJPEG("/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210410_wikpedia_blue_whale.jpg")
# img_bd <- readJPEG("/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210410_wikpedia_bottlenose_dolphin.jpg")
img_bw <- readJPEG("/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210421_original_blue_whale.jpg")
img_bd <- readJPEG("/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210421_original_bottlenose_dolphin.jpg")



# get ggplot items and fix aspect rations
p_img_bd <- ggplot() + background_image(img_bd) + theme(plot.margin = margin(t=0.0, l=0.25, r=0.25, b=0.0, unit = "cm")) + coord_fixed(ratio=0.75)
p_img_bw <- ggplot() + background_image(img_bw) + theme(plot.margin = margin(t=0.0, l=0.25, r=0.25, b=0.0, unit = "cm")) + coord_fixed(ratio=0.75)

ggarrange(
  ggarrange(pPhl, pCls, pOrd, pFam, pGen, pSpc, ncol = 1, nrow = 6, labels = c("(a)")),
  ggarrange(p_heatobs, ncol = 1, nrow = 1, labels = c("(b)")),
  ggarrange(p_nmds,  p_img_bd, p_img_bw, ncol = 1, nrow = 3, labels = c("(c)", "(d)", "(e)")),
  ncol = 3, nrow = 1, widths = c(1, 2, 2))

ggsave("210407_998_r_summarize_results_fig2_draft_Venn_heat.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 200, height = 200, units = c("mm"),
         dpi = 500, limitsize = TRUE)




