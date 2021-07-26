#  **********************************
#  *                                *        
#  *  Get results from long tables  *
#  *                                *
#  **********************************

lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
rm(list = ls(all.names = TRUE))
gc()
options(tibble.print_max = Inf) 

# I. Load packages and define functions
# =====================================

library("tidyverse")   # because we can't stop using it anymore
library("magrittr")    # get the %<>% pipe

library("taxize")      # look up trivial names
Sys.setenv(ENTREZ_KEY="a634c6e9c96c3859bca27a2771f6d2872f08")
Sys.getenv("ENTREZ_KEY")


library("ggpubr") # combine plots -  http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/


library("sf")           # simple feature objects

library("flextable")   # format species lists as tables https://ardata-fr.github.io/flextable-book/
library("officer")     # format species lists as tables https://ardata-fr.github.io/flextable-book/
library("magick")      # convert flext table to ggplot grob
library("grid")        # convert flext table to ggplot grob
# library("cowplot")   # convert flext table to ggplot grob



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

# not in 
# -------
'%!in%' <- function(x,y)!('%in%'(x,y))


# get Euler objects for plotting
# ------------------------------
get_euler_object = function(level, tibl){
  require("eulerr")
  require("tidyverse")
  require("magrittr")
  
  # check if needed columns are in the input data
  stopifnot(c("BRUV.OBS.PRES", "EDNA.OBS.PRES", "OBIS.OBS.PRES", "PUBL.OBS.PRES") %in% names(tibl))
  stopifnot(level %in% c("SUPERKINGDOM", "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"))
  
  # isolate realvant columns for summary
  tibl %<>% select(BRUV.OBS.PRES, EDNA.OBS.PRES, OBIS.OBS.PRES, PUBL.OBS.PRES, SUPERKINGDOM,  PHYLUM,  CLASS,  ORDER,  FAMILY,  GENUS, SPECIES) %>% distinct()
  
  # debugging only
  # print(tibl)
  
  # sum up unique presences fr Euler plot
  tibl %<>% group_by(get(level)) %>% summarise(eDNA = as.logical(sum(EDNA.OBS.PRES)),
                                               BRUV = as.logical(sum(BRUV.OBS.PRES)),
                                               PUBL = as.logical(sum(PUBL.OBS.PRES)),
                                               OBIS = as.logical(sum(OBIS.OBS.PRES))
                                               )
  # debugging only
  print(tibl)
  
  return(euler(tibl[ , 2:5]))
}

# get Euler Ggplots 
# -----------------
get_euler_ggplot = function(level, euler_ob, plot_label = TRUE){
  require("tidyverse")
  require ("ggplotify")

  # sanitize input
  stopifnot( class(euler_ob)[1] == "euler")
  stopifnot(level %in% c("SUPERKINGDOM", "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"))
  
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

# get matrix - for distance calculations and 
get_matrix_or_table <- function(tibl, group_col = "RESERVE.GROUP.LOCATION", group_row = "SPECIES", obs_methods = NULL, tbl = FALSE){
  
  # debugging
  #   tibl = get("fish_biodiv")
  #   group_col = "RESERVE.GROUP.LOCATION"
  #   group_row = "SPECIES"
  #   obs_methods = c("eDNA")
  #   tbl = FALSE
  
  stopifnot(group_col %in% names(tibl))
  stopifnot(group_row %in% c("SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"))
  # stopifnot(length(obs_methods) == 1)
  stopifnot(obs_methods %in% c("BRUV", "eDNA", "OBIS", NULL))
  
  require(data.table) # re-use old code rather then finding out how to reimplement
  require(tidyverse)
  
  # subset for desired observation method if needed
  
  if(!is.null(obs_methods))  { 
    message("Keeping only observations of \"", obs_methods, "\".")
    
    # debugging checks
    #   filter(tibl, SAMPLE.TYPE %in% "BRUV") %>% select(SET.ID, REP.ID, SAMPLE.TYPE,RESERVE.GROUP.LOCATION, SPECIES, ANY.OBS.PRES, EDNA.OBS.PRES, BRUV.OBS.PRES, OBIS.OBS.PRES)
    #   filter(tibl, SAMPLE.TYPE %in% "eDNA") %>% select(SET.ID, REP.ID, SAMPLE.TYPE,RESERVE.GROUP.LOCATION, SPECIES, ANY.OBS.PRES, EDNA.OBS.PRES, BRUV.OBS.PRES, OBIS.OBS.PRES)
    #   filter(tibl, SAMPLE.TYPE %in% "OBIS") %>% select(SET.ID, REP.ID, SAMPLE.TYPE,RESERVE.GROUP.LOCATION, SPECIES, ANY.OBS.PRES, EDNA.OBS.PRES, BRUV.OBS.PRES, OBIS.OBS.PRES)
  
    tibl <- filter(tibl, SAMPLE.TYPE %in% obs_methods)
    message("Keeping:", tibl %>% ungroup %>% select(SAMPLE.TYPE) %>% distinct)
  }
  
  # filter for defined observations - NAs would get lumped together and give erroneous NMDS results
  message("Filtering out NA's of level ", group_row, " for distance calculations, taxon list will likely have shortened.")
  tibl <- filter(tibl, !is.na(get(group_row)))
  
  dtbl <-  as.data.table(tibl)   # re-use old code rather then finding out how to reimplement
  
  # aggregate table for provided grouping variables,
  #   get counts of occurrences per location, per taxonomic entity for subsequent matrix 
  #   not sure if all taxonomy level are needed, implementing just in case
  if (group_row == "SPECIES") {
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "NCBI.TAXID", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES", "TRIVIAL.SCPECIES"), .SDcols=c("ANY.OBS.PRES") ]
  } else if (group_row == "GENUS") {
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY",  "GENUS"), .SDcols=c("ANY.OBS.PRES") ]
  } else if (group_row == "FAMILY") {
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY"), .SDcols=c("ANY.OBS.PRES") ]
  } else if (group_row == "ORDER") {
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER"), .SDcols=c("ANY.OBS.PRES") ]
  } else if (group_row == "CLASS") {
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "SUPERKINGDOM", "PHYLUM", "CLASS"), .SDcols=c("ANY.OBS.PRES") ]
  } else {
    stop("\"group_row\" needs to be on of  \"CLASS\", \"ORDER\", \"FAMILY\", \"GENUS\", \"SPECIES\"")
  }
 
  
  # return matrix for NMDS as default - otherwise list
  if(tbl == TRUE)  {
    
    # print(as_tibble(dtbl_ag))
    return(as_tibble(dtbl_ag)) 
  
  } else if (tbl == FALSE) {
    
    dtbl_ag <- dcast(setDT(dtbl_ag), get(group_col)~get(group_row), value.var="ANY.OBS.PRES", sum, fill=0)
    mat <- as.matrix(dtbl_ag, rownames=TRUE)
    message("Matrix dimensions are", paste(dim(mat), collapse =  " "))
    return(mat)
  }
}

# get ANOSIM results for a given tibble ("tibl")
#   aggregate observations for "group_col" on taxonomic level "group_row"
#   aggregated observations ar replicates for ANOSIM as levels of "group_col_ano"
get_anosim <- function(tibl, group_col = NULL, group_row = NULL, group_col_ano = NULL, obs_methods = NULL, distance = NULL){
  
  # debugging
  #  tibl          = get(c("full_biodiv"))
  #  group_col     = c("SET.ID") 
  #  group_row     = c("SPECIES", "GENUS", "FAMILY", "ORDER", "CLASS")[1]
  #  group_col_ano = c("RESERVE.GROUP.LOCATION", "RESERVE.GROUP.INSIDE")[1] 
  #  obs_methods   = c("OBIS")
  #  distance      = "jaccard"
  
  # debugging 
  #  message("get_anosim() uses \"group_col\": ", group_col)
 
  # get matrix and group variables for Anosim
  ano_mat <- get_matrix_or_table(tibl, group_col, group_row, tbl = FALSE, obs_methods = obs_methods)
  ano_grp <- as_tibble(rownames(ano_mat)) %>% 
    rename(!!group_col := value)  %>% 
    merge(., unique(tibl[ ,c(as.character(group_col),group_col_ano ) ]))
  
  # debugging
  #  message(ano_grp)
    
  # prepare grouping in correct order 
  ano_grp <-  ano_grp[match(rownames(ano_mat), ano_grp[[group_col]]),]
   
  # isolate group column
  ano_grp <- ano_grp %>% pull(!! group_col_ano)
  
  # debugging
  #   message("get_anosim() is generated grouping var: ", paste(ano_grp, collapse = " "))
  #   message("                   for matrix rownames: ", paste(rownames(ano_mat), collapse = " "))
   
  # if not called from dataframes line needs to be: 
  # filter(get(group_col) %in% rownames(ano_mat)) %>% arrange(interp(group_col)) %>% pull(group_col_ano)

  
  # get vegan results
  anosim <- vegan::anosim(ano_mat, ano_grp, permutations = 9999, distance = distance)
  
  return(anosim)

}

# II. Read in data
# ================

# check input data of previous script
system("open -a \"Microsoft Excel\" \"/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/998_r_map_and_add_obis__full_data_raw.xlsx\"")

long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_map_and_add_obiss__full_data_raw.Rds")
long_table <- ungroup(long_table)

# III. Format data  
# =================
# - add trivial names 
# - mark non-NZ species  **(possibly needs to be re-worked)**
# - split "fish" and "full" data
# - filter for data completeness **(possibly needs to be re-worked)**

# clean taxonomy strings
# ----------------------

# erase all that may come after "sp." if there is an "sp." at all
long_table %<>% mutate(SPECIES = gsub("sp\\..*", "sp.", SPECIES))
long_table %<>% mutate(SPECIES = gsub("cf. ", "", SPECIES))
long_table %<>% mutate(SPECIES = str_replace(SPECIES, "\\s\\S*\\s\\S*(.*)", ""))

long_table %<>% mutate(ORDER = ifelse(ORDER == "Squalidae", "Squaliformes", ORDER))
long_table %<>% mutate(ORDER = ifelse(GENUS == "Callanthias", "Perciformes", ORDER))

# add trivial names 
# -----------------

# get trivial names
species_vector <- ungroup(long_table) |> select(SPECIES) |> distinct() |> pull(SPECIES)
trivial_obj <- sci2comm(species_vector, db = "ncbi", simplify = TRUE)

# format trivial names for left join
trivial_df <- data.frame(do.call(rbind , trivial_obj))
trivial_df <- rownames_to_column(trivial_df, var = "SPECIES")
trivial_df <- trivial_df |>  setNames( c("SPECIES", "TRIVIAL.SPECIES")) |> as_tibble()

# add trivial names to object 
long_table %<>% left_join(trivial_df)

#  mark non-NZ species  **(possibly needs to be re-worked)**
# ---------------------------------------------------------
#  16-Mar-2021 add asterisks ("*") to non-NZ species, and ("**") to non-fish (mammals and crustaceans)
#  after checking with list 
#  Roberts, C., Stewart, A., Struthers, C., Barker, J. & Kortet, S. 2019 Checklist of the Fishes of New Zealand. 

nonnz_fish <- c("Asterropteryx", "Banjos", "Benitochromis", "Bostrychus", "Bovichtus", "Caprodon", "Coptodon", "Engraulis", "Gobiesox", "Gymnoscopelus", "Helcogramma", "Microcanthus", "Opistognathus", "Phoxinus", "Sander", "Scobinichthys")
nonnz_othr <- c("Macroctopus", "Jasus", "Arctocephalus", "Balaenoptera", "Tursiops")

long_table %<>% mutate(SPECIES = 
                         case_when(GENUS %in% nonnz_fish ~ paste0(SPECIES, "*"),
                                   GENUS %in% nonnz_othr ~ paste0(SPECIES, "**"),
                                                  TRUE ~ SPECIES)
                                                   )
long_table %<>% mutate(GENUS = 
                         case_when(GENUS %in% nonnz_fish ~ paste0(GENUS, "*"),
                                   GENUS %in% nonnz_othr ~ paste0(GENUS, "**"),
                                                    TRUE ~ GENUS)
                                                    )
# save / load annotated object
# ------------------------------

save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__start_env.Rdata")
saveRDS(long_table, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_summarize_results__full_data_rev.Rds")
long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_summarize_results__full_data_rev.Rds")

# Filter for data completeness **(possibly needs to be re-worked)**
# ------------------------------------------------------------------

# - not done yet -

#  get equivalent of BOTH.PRES
# ------------------------------

# function possibly needs to 
#   get presence / absence on a {taxonomic level} (SPECIES)
#   per a {location} (SET.ID  RESERVE.GROUP.LOCATION )
#   ? check data completeness  {all "1" in BRUV.OBS.PRES	EDNA.OBS.PRES	OBIS.OBS.PRES}
#   implementing as follows

# for aggregation of identical taxonomic entities across multiple observations, copy observation typ to one
#   uniting variable (used to be "BOTH.PRES")
long_table %<>% mutate(ANY.OBS.PRES = case_when(BRUV.OBS.PRES == 1 ~ 1, 
                                                EDNA.OBS.PRES == 1 ~ 1, 
                                                OBIS.OBS.PRES == 1 ~ 1,
                                                PUBL.OBS.PRES == 1 ~ 1,
                                                TRUE ~ 0))

# Split "fish" and "full" data
# ----------------------------
full_biodiv <- long_table %>% distinct()
fish_biodiv <- long_table %>% distinct() %>% filter(CLASS %in% c("Actinopteri", "Chondrichthyes")) %>% filter(!(GENUS %in% c("Sardinops")))

save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__data_filtered.Rdata")


# IV. What data is available for Fiordland - table summary for supplemnet
# ========================================================================

# format data for flex table
obs_sums <- fish_biodiv %>% ungroup() %>% group_by(SPECIES) %>%  
  select(SPECIES, BRUV.OBS.PRES, EDNA.OBS.PRES, OBIS.OBS.PRES, PUBL.OBS.PRES) %>%   
  summarize(BRUV.OBS.PRES.SUM = ifelse(sum(BRUV.OBS.PRES), sum(BRUV.OBS.PRES), NA),
            EDNA.OBS.PRES.SUM = ifelse(sum(EDNA.OBS.PRES), sum(EDNA.OBS.PRES), NA),
            OBIS.OBS.PRES.SUM = ifelse(sum(OBIS.OBS.PRES), TRUE, NA), 
            PUBL.OBS.PRES.SUM = ifelse(sum(PUBL.OBS.PRES), TRUE, NA)
            )
 
spcies_obs_sums <- ungroup(fish_biodiv) %>% select(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, TRIVIAL.SPECIES) %>% distinct() %>%
  left_join(obs_sums) %>% arrange(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES) 

# generate flex table - for merging with ggplot object
ft_spcies_obs_sums <-  flextable(spcies_obs_sums) %>% 
   merge_v(j = "PHYLUM", target = "PHYLUM") %>%
   merge_v(j = "CLASS", target = "CLASS") %>%
   merge_v(j = "ORDER", target = "ORDER") %>%
   merge_v(j = "FAMILY", target = "FAMILY") %>%
   merge_v(j = "GENUS", target = "GENUS") %>% 
   valign(valign = "top") %>%
   italic(j = c(5,6)) %>%
   bold(j = c(7)) %>% 
   set_header_labels(values = list(PHYLUM = "Phylum",
     CLASS = "Class",
     ORDER = "Order",
     FAMILY = "Family",
     GENUS = "Genus",
     SPECIES = "Species",
     TRIVIAL.SPECIES = "Common name",
     BRUV.OBS.PRES.SUM = "BRUV observations",
     EDNA.OBS.PRES.SUM = "eDNA observations",
     OBIS.OBS.PRES.SUM = "OBIS records", 
     PUBL.OBS.PRES.SUM = "Literature records" )) %>%
     fontsize(part = "all", size = 9) %>%
     fit_to_width(max_width = 12, inc = 1L, max_iter = 20)

save_as_docx(ft_spcies_obs_sums, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210712_998_r_summarize_results__all_data.docx",
  pr_section = prop_section(
    page_size = page_size(orient = "portrait"), type = "continuous"
    ))

save_as_html(ft_spcies_obs_sums, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210712_998_r_summarize_results__all_data.html")


nrow(spcies_obs_sums)                                    # found 116 species across all data sets
nrow(spcies_obs_sums |> filter (CLASS == "Actinopteri")) #       106 Actinopteri
nrow(spcies_obs_sums |> filter (CLASS == "Chondrichthyes")) #     10 Chondrichthyes

nrow(spcies_obs_sums |> filter (!is.na(BRUV.OBS.PRES.SUM)))  # 25 BRUV (in study area)
nrow(spcies_obs_sums |> filter (!is.na(EDNA.OBS.PRES.SUM)))  # 44 EDNA (in study area)
nrow(spcies_obs_sums |> filter (!is.na(OBIS.OBS.PRES.SUM)))  # 25 OBIS (in circle)
nrow(spcies_obs_sums |> filter (!is.na(PUBL.OBS.PRES.SUM)))  # 59 PUBL (Fiordland)


# III. Get Euler plots
# ====================

# get euler analysis results for plotting / plot_label = TRUE shrinks plots a lot
# euler_obs_full_bio <- lapply(list("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"), get_euler_object, full_biodiv)
# euler_ggp_full_bio <- mapply(get_euler_ggplot, list("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"),  euler_obs_full_bio, plot_label = FALSE, SIMPLIFY = FALSE)

# plot euler analysis results
euler_obs_fish_bio <- lapply(list("SUPERKINGDOM", "PHYLUM", "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"), get_euler_object, fish_biodiv)
euler_ggp_fish_bio <- mapply(get_euler_ggplot, list("SUPERKINGDOM", "PHYLUM", "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES"),  euler_obs_fish_bio, plot_label = FALSE, SIMPLIFY = FALSE)

# create compound plot with better labels then with plot_label = TRUE above
ggarrange( plotlist = euler_ggp_fish_bio[4:7],
           labels = str_to_sentence(c("ORDER",  "FAMILY",  "GENUS", "SPECIES")),
           font.label = list(size = 12, color = "black", face = "bold.italic", family = NULL),
           ncol = 2, nrow = 2)
           )

# save compound plot with better labels then with plot_label = TRUE above
ggsave("210712_998_r_summarize_results__euler_edna_bruv_obis.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 2.0, width = 100, height = 100, units = c("mm"),
         dpi = 500, limitsize = TRUE)  


# continue here after 26-Jul-2021




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
      geom_sf(data = bbox_rgl_full_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) +
      geom_sf_label(data=bbox_rgl_full_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
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
      geom_sf(data = bbox_rgl_fish_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) + 
      geom_sf_label(data=bbox_rgl_fish_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
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
         scale = 2.5, width = 85, height = 85, units = c("mm"),
         dpi = 500, limitsize = TRUE)  
         
# see saveRDS(map_main, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_get_OBIS_and_map__mapggplot.Rds")

# Omitted - Get a tree 
# =====================

# - omitted - 

# use reticulate to call Python from within R (see https://cran.r-project.org/web/packages/reticulate/vignettes/calling_python.html)
# use Python's ETE toolkit (see http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html)
# ETE toolkit installed in conda environment "ete3" (http://etetoolkit.org/download/)

# Sys.setenv(RETICULATE_PYTHON = "/Users/paul/Applications/miniconda3/envs/ete3/bin/python")
# reticulate::use_condaenv("/Users/paul/Applications/miniconda3/envs/ete3",  required = TRUE) # ... in the righ environmnet...
# reticulate::py_config()
# source_python("/Users/paul/Documents/OU_eDNA/200901_scripts/get_tree_for_ncbi_taxid_vector.py") #  ...using the corrcet function
# 
# tax_ids <- c(9606, 9598, 10090, 7707, 8782)
# 
# # call python function (R multi-element vector should become list automatically)
# get_tree_for_ncbi_taxid_vector(tax_ids)
# 
# 
# # print obtained Newick tree using ggtree
# #   https://www.molecularecologist.com/2017/02/08/phylogenetic-trees-in-r-using-ggtree/


# V. Get biodiversity heat map and matching flex table
# ====================================================

# get plotting data sets
# ----------------------

htmp_tibl_fish <- bind_rows(
  get_matrix_or_table(fish_biodiv, obs_methods = "eDNA",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "eDNA"),
  get_matrix_or_table(fish_biodiv, obs_methods = "BRUV",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "BRUV"), 
  get_matrix_or_table(fish_biodiv, obs_methods = "OBIS",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "OBIS"), 
)

htmp_tibl_full <- bind_rows(
  get_matrix_or_table(full_biodiv, obs_methods = "eDNA",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "eDNA"),
  get_matrix_or_table(full_biodiv, obs_methods = "BRUV",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "BRUV"), 
  get_matrix_or_table(full_biodiv, obs_methods = "OBIS",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "OBIS"), 
)

# get margin sums 
# h_total <- long_table_dt_agg_gen %>% 
#   group_by(GENUS) %>% 
#   summarise(BOTH.PRES = sum(BOTH.PRES)) %>% 
#   mutate(RESERVE.GROUP.LOCATION = 'TOTAL')
# 
# v_total <- long_table_dt_agg_gen %>% 
#   group_by(RESERVE.GROUP.LOCATION) %>% 
#   summarise(BOTH.PRES = sum(BOTH.PRES)) %>% 
#   mutate(GENUS = 'TOTAL')
#   
# add margin sums column for merging
# v_total <- v_total %>% mutate(GENUS.SUM = BOTH.PRES) %>% select(RESERVE.GROUP.LOCATION, GENUS.SUM)
# h_total <- h_total %>% mutate(RESERVE.GROUP.LOCATION.SUM = BOTH.PRES) %>% select(GENUS, RESERVE.GROUP.LOCATION.SUM)

# aggregate discrete observation of wither method ("BOTH.PRES") per sampling area (RESERVE.GROUP.LOCATION) on GENUS level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
# long_table_dt_agg_gen <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]
# 
# RESERVE.GROUP.LOCATION.LABS <- list(
#   "LS MR"   = "LS MR\n (n = 4)",
#   "LS CTRL" = "LS CTRL\n (n = 4)",
#   "FF MR"   = "FF MR\n (n = 2)",
#   "FF CTRL" = "FF CTRL\n (n = 3)",
#   "WJ MR"   = "WJ MR\n (n = 4)",
#   "WJ CTRL" = "WJ CTRL\n (n = 4)"
# )

# get matching flex table
# -----------------------

# format data for flex table
tibl_plot <- htmp_tibl_full %>% select(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, TRIVIAL.SCPECIES) %>% 
 distinct() %>% arrange(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES) 
 
 
# generate flex table - for merging with ggplot object
ft <-  flextable(tibl_plot) %>% 
   merge_v(j = "PHYLUM", target = "PHYLUM") %>%
   merge_v(j = "CLASS", target = "CLASS") %>%
   merge_v(j = "ORDER", target = "ORDER") %>%
   merge_v(j = "FAMILY", target = "FAMILY") %>%
   merge_v(j = "GENUS", target = "GENUS") %>% 
   valign(valign = "top") %>%
   italic(j = c(5,6)) %>%
   bold(j = c(7)) %>% 
   set_header_labels(values = list(PHYLUM = "Phylum",
     CLASS = "Class",
     ORDER = "Order",
     FAMILY = "Family",
     GENUS = "Genus",
     SPECIES = "Species",
     TRIVIAL.SCPECIES = "Common name")) %>%
     fontsize(part = "all", size = 12) %>%
     autofit(add_w = 1, add_h = 0, part = c("all"))

# order factors in plotting object to match flex table  
y_axis_label_order <- fct_relevel(htmp_tibl_full$SPECIES, rev(c(tibl_plot$SPECIES)))  # y_axis_label_order <- reorder(htmp_tibl_full$SPECIES, desc(htmp_tibl_full$SPECIES))


# plot out data as tiles
# -----------------------
# - needs vectors ordered with tree
# - needs annotation of non NZ species
# - may need margin sums

plot_htmp <- ggplot(htmp_tibl_full, aes_string(x = "RESERVE.GROUP.LOCATION", y = y_axis_label_order)) +
    geom_tile(aes(fill = ANY.OBS.PRES) ) +
    scale_fill_gradient(low="black", high="red") +
    geom_text(size = 2, aes(label = ANY.OBS.PRES, color = "white")) +
    facet_grid(. ~ SAMPLE.TYPE, scales = "fixed") + 
    theme_bw() +
    theme(legend.position = "none", 
          strip.text.y = element_text(angle=0), 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7, face = "italic"), 
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()
          ) +
    xlab("Sampling Locations") # + ylab("Species Observations")

# Combine heatmap and flextable
# ------------------------------

# get flex table as ggplot object
ft_raster <- as_raster(ft) # webshot and magick.
ft_rgrob  <- rasterGrob(ft_raster)
plot_ft <- ggplot() + theme_void() + annotation_custom(ft_rgrob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

ggarrange(plot_ft, plot_htmp, ncol = 2, nrow = 1, labels = c("a","b"), widths = c(1,1))

ggsave("210712_998_r_summarize_results__biodiv_heat.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 2.5, width = 85, height = 85, units = c("mm"),
         dpi = 500, limitsize = TRUE)  


# VI. ANOSIM of observation types and variables
# ===============================================

# testing function with one data set
get_anosim(distance = "jaccard", tibl = full_biodiv, group_col = "SET.ID", group_row = c("SPECIES"), group_col_ano = "RESERVE.GROUP.LOCATION", obs_methods = "BRUV")


# Analysis for fish (as per eDNA markers) for BRUV and eDNA (data complete across all sets)
# ----------------------------------------------------------------------------------------------

# setting up parameter combinations for complete ANOSIM analysis
anosim_analysis_fish <- expand.grid(
  distance      = c("jaccard"),
  tibl          = c("fish_biodiv"),
  group_col     = c("SET.ID"), 
  group_row     = c("SPECIES", "GENUS", "FAMILY", "ORDER", "CLASS"),
  group_col_ano = c("RESERVE.GROUP.LOCATION", "RESERVE.GROUP.INSIDE"), 
  obs_methods   = c("eDNA", "BRUV")
  )

# run ANOSIM analysis
anosim_results_fish <- apply(anosim_analysis_fish, 1, FUN = function(x) try(get_anosim(distance = x[1], tibl = get(x[2]), group_col = x[3], group_row = x[4], group_col_ano = x[5], obs_methods = x[6])))

# inspect results
#  str(anosim_results[[1]])

# get results data frame
anosim_analysis_fish$statistic  <- unlist(lapply(anosim_results_fish, '[[', 5))
anosim_analysis_fish$significance <- unlist(lapply(anosim_results_fish, '[[', 2))
anosim_analysis_fish
anosim_analysis_fish %>% filter(significance <= 0.05 )


# Analysis only full biodiversity (for BRUV, eDNA, and OBIS) 
# -----------------------------------------------------------
#   with truntaed dataset (OBIS observations missing for some SET.ID's (22,27,28,29))
#   --------------------------------------------------------------------------------

# filter out data undefined for OBIS 
full_biodiv_obis <- full_biodiv %>% filter(SET.ID %!in% c(22,27,28,29))

# testing function with one data set
get_anosim(distance = "jaccard", tibl = full_biodiv_obis, group_col = "SET.ID", group_row = c("SPECIES"), group_col_ano = "RESERVE.GROUP.LOCATION", obs_methods = "OBIS")

anosim_analysis_full_obis <- expand.grid(
  distance      = c("jaccard"),
  tibl          = c("full_biodiv_obis"),
  group_col     = c("SET.ID"), 
  group_row     = c("SPECIES", "GENUS", "FAMILY", "ORDER", "CLASS"),
  group_col_ano = c("RESERVE.GROUP.LOCATION", "RESERVE.GROUP.INSIDE"), 
  obs_methods   = c("eDNA", "BRUV","OBIS")
  )

# run ANOSIM analysis
anosim_results_full_obis <- apply(anosim_analysis_full_obis, 1, FUN = function(x) try(get_anosim(distance = x[1], tibl = get(x[2]), group_col = x[3], group_row = x[4], group_col_ano = x[5], obs_methods = x[6])))

# get results data frame
anosim_analysis_full_obis$statistic    <- unlist(lapply(anosim_results_full_obis, '[[', 5))
anosim_analysis_full_obis$significance <- unlist(lapply(anosim_results_full_obis, '[[', 2))
anosim_analysis_full_obis
anosim_analysis_full_obis %>% filter(significance <= 0.05 )


# Analysis full biodiversity (for BRUV, eDNA) 
# -----------------------------------------------------------

# testing function with one data set
get_anosim(distance = "jaccard", tibl = full_biodiv, group_col = "SET.ID", group_row = c("GENUS"), group_col_ano = "RESERVE.GROUP.LOCATION")

anosim_analysis_full <- expand.grid(
  distance      = c("jaccard"),
  tibl          = c("full_biodiv"),
  group_col     = c("SET.ID"), 
  group_row     = c("SPECIES", "GENUS", "FAMILY", "ORDER", "CLASS"),
  group_col_ano = c("RESERVE.GROUP.LOCATION", "RESERVE.GROUP.INSIDE"), 
  obs_methods   = c("eDNA", "BRUV")
  )

# run ANOSIM analysis
anosim_results_full <- apply(anosim_analysis_full, 1, FUN = function(x) try(get_anosim(distance = x[1], tibl = get(x[2]), group_col = x[3], group_row = x[4], group_col_ano = x[5], obs_methods = x[6])))

# get results data frame
anosim_analysis_full$statistic    <- unlist(lapply(anosim_results_full, '[[', 5))
anosim_analysis_full$significance <- unlist(lapply(anosim_results_full, '[[', 2))
anosim_analysis_full
anosim_analysis_full %>% filter(significance <= 0.05 )

































# ###############################################################
# obis obs missing for some set.ids (22,27,28,29) - annotate in map 
# not done yet
full_biodiv %>% summarise(sum(OBIS.OBS.PRES)) %>% print(n = Inf)
# ###############################################################


















# IX. Indicator species analysis for significant ANOMSIM results
# ==============================================================


# - not done yet - 




pmap(list(x = 1, y = 1, z = 1 ), sum, na.rm = TRUE)

mapply(rep, 1:4, 4:1)
mapply(rep, times = 1:4, MoreArgs = list(x = 42))









# ******* OLD code below **********
# =================================









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




