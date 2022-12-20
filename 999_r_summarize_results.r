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

library("tidyverse")	# because we can't stop using it anymore
library("magrittr")		# get the %<>% pipe
library("taxize")		  # look up trivial names
library("ggpubr")		  # combine plots -  http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
library("sf")			    # simple feature objects
library("flextable")	# format species lists as tables https://ardata-fr.github.io/flextable-book/
library("officer")		# format species lists as tables https://ardata-fr.github.io/flextable-book/
library("magick")		  # convert flext table to ggplot grob
library("grid")			  # convert flext table to ggplot grob
library("sjPlot")		  # model plots

# API key for Taxcise - trivial names lookup
Sys.setenv(ENTREZ_KEY="YourNcbiApiKeyHere")
Sys.getenv("ENTREZ_KEY")

# API key for Fishbase - synonyms lookup
Sys.setenv(FISHBASE_API="sealifebase")
Sys.getenv("FISHBASE_API")

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
  bbox <- tibl %>% 
  group_by(across(all_of(location))) %>%
  summarise(xmin = min(MH.PPS.LONG) -0.01 ,ymin = min(MH.GPS.LAT) -0.01, xmax=max(MH.PPS.LONG) +0.01,  ymax = max(MH.GPS.LAT) +0.01) %>%
  gather(x,lon,c('xmin','xmax')) %>% gather(y,lat,c('ymin','ymax')) %>%
  st_as_sf(coords=c('lon','lat'),crs=4326,remove=F) %>%
  group_by(across(all_of(location))) %>% mutate(angle = calc_angle(lon,lat)) %>%
  arrange(angle) %>% 
  summarise(do_union=FALSE) %>% 
  st_cast('POLYGON')
  
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
  stopifnot(obs_methods %in% c("BRUV", "eDNA", "OBIS", "PUBL", NULL))
  
  require("data.table") # re-use old code rather then finding out how to reimplement
  require("tidyverse")
  
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
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "NCBI.TAXID", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES", "TRIVIAL.SPECIES"), .SDcols=c("ANY.OBS.PRES") ]
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
    message("Matrix dimensions are ", paste(dim(mat), collapse =  " x "))
    return(mat)
  }
}

# get ANOSIM results for a given tibble ("tibl")
#   aggregate observations for "group_col" on taxonomic level "group_row"
#   aggregated observations ar replicates for ANOSIM as levels of "group_col_ano"
get_vegan <- function(tibl, group_col = NULL, group_row = NULL, group_col_ano = NULL, obs_methods = NULL, distance = NULL, mp = FALSE){
  
  require("vegan")        # for NMDS 
  require("indicspecies") # indicator species  - see citation below

  
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

  # get analysis results - ANOSIM unless Indicator species analysis is requested at function start
  if (mp == FALSE) {
    
    anosim <- vegan::anosim(ano_mat, ano_grp, permutations = 9999, distance = distance)
    return(anosim)
    
  } else {
    
    mpl <- multipatt(ano_mat, ano_grp, func = "r.g", control = how(nperm=9999))
    return(mpl)
    
  }

}

# II. Read in data
# ================

# check input data of previous script
system("open -a \"Microsoft Excel\" \"/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/999_r_map_and_add_obis__full_data_raw.xlsx\"")

long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/999_r_map_and_add_obis__full_data_raw.Rds")
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

# 20.12.2022 - correcting for synonyms, reinspecting table for sanity
# -------------------------------------------------------------------

# isolate species names from analysis data
spc_read <- long_table %>% pull("SPECIES") %>% unique() %>% sort()

# fetch synonyms from Fishbase (20.12.2022)
spc_in_syn <- rfishbase::synonyms(species_list = spc_read)

# save state after fishbase lookup
save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__syn_lookup.Rdata")

# inspect synonyms
spc_in_syn %>% select(synonym, Species, Status)

# isolate fishbase dat that needs to be used to corrcet data in long table
fb_data_for_correction <- spc_in_syn %>% filter(Status == "synonym") %>% filter(synonym != "Opistognathus sp.") 

# inspect those data again
fb_data_for_correction %>% select(synonym, Species, Status)

# correct synonymous species and genus names in `long_table`, that can be named properly

# [continue here after 20.12.2022]


# [older code below]

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


# III. What data is available for Fiordland - table summary for supplement
# =======================================================================

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

# Summary: general species counts 
nrow(spcies_obs_sums)                                    # found 116 species across all data sets
nrow(spcies_obs_sums |> filter (CLASS == "Actinopteri")) #       106 Actinopteri
nrow(spcies_obs_sums |> filter (CLASS == "Chondrichthyes")) #     10 Chondrichthyes

nrow(spcies_obs_sums |> filter (!is.na(BRUV.OBS.PRES.SUM)))  # 25 BRUV (in study area)
nrow(spcies_obs_sums |> filter (!is.na(EDNA.OBS.PRES.SUM)))  # 44 EDNA (in study area)
nrow(spcies_obs_sums |> filter (!is.na(OBIS.OBS.PRES.SUM)))  # 25 OBIS (in circle)
nrow(spcies_obs_sums |> filter (!is.na(PUBL.OBS.PRES.SUM)))  # 59 PUBL (Fiordland)

# Summary: OBIS data in small circles 
fish_biodiv |> pull(SET.ID) |> unique() |> sort()
fish_biodiv |> 
  filter(SET.ID %!in% c(98, 99)) |>
  filter(SAMPLE.TYPE == "OBIS") |> 
  select(SPECIES, SAMPLE.TYPE, RESERVE.GROUP.LOCATION) |> 
  distinct() |> pull(RESERVE.GROUP.LOCATION) |> unique()

fish_biodiv |> 
  filter(SET.ID %!in% c(98, 99)) |>
  filter(SAMPLE.TYPE == "OBIS") |> 
  select(SET.ID, SPECIES, SAMPLE.TYPE, RESERVE.GROUP.LOCATION) |> 
  distinct() |> pull(SET.ID) |> unique()  # 7  9 11 12 17 18 19 21 26
  

fish_biodiv |> 
  filter(SET.ID %!in% c(98, 99)) |>
  filter(SAMPLE.TYPE == "OBIS") |> 
  select(SPECIES, SAMPLE.TYPE, RESERVE.GROUP.LOCATION) |> 
  distinct(SPECIES)

# 1 Acanthoclinus matti     
# 2 Lepidoperca tasmanica   
# 3 Pseudolabrus miles      
# 4 Parapercis colias       
# 5 Notolabrus cinctus      
# 6 Notolabrus fucicola     
# 7 Notolabrus celidotus    
# 8 Nemadactylus macropterus
# 9 Latris lineata          
# 10 Pseudophycis barbata    
# 11 Scorpaena papillosa     
# 12 Carcharodon carcharias

# Summary: Literture species counts
fish_biodiv |> filter(SET.ID %in% c(98, 99)) |> select(SPECIES) |> distinct()

# Summary: Fish in BRUV that ar not in literture
bruv_species <- fish_biodiv |> 
  select(SPECIES, SAMPLE.TYPE) |> filter(SAMPLE.TYPE == "BRUV") |> 
  distinct() |> pull(SPECIES)

publ_species <- fish_biodiv |> 
  select(SPECIES, SAMPLE.TYPE) |> filter(SAMPLE.TYPE %!in% c("BRUV", "eDNA")) |> 
  distinct() |> pull(SPECIES)

bruv_species[bruv_species %!in% publ_species] |> sort()

# IV. Get Euler plots
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
           
# save compound plot with better labels then with plot_label = TRUE above
ggsave("210712_998_r_summarize_results__euler_edna_bruv_obis.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 2.0, width = 100, height = 100, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# V. Get geographical maps with heat overlays
# =============================================

# compare script  ~/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r

# data preparation
# ----------------

# keep only data relevant for plotting 
fish_biodiv_gh <- fish_biodiv |> filter(SAMPLE.TYPE %in% c("eDNA","BRUV","OBIS") & SET.ID %!in% c(98,99))

# for mapping: get column-subset sf's with WGS 84 in degrees 
fish_biodiv_sf <- get_sf_biodiv(fish_biodiv_gh)

# for mapping: get map layers 
nzshp_hires_WGS84_sf <- read_sf("/Users/paul/GIS/NZ_coast/NZ_Coast_isl.shp") %>% st_transform(crs = 4326)
nzshp_lores_WGS84_sf <- rmapshaper::ms_simplify(input = as(nzshp_hires_WGS84_sf, 'Spatial')) %>% st_as_sf

# for mapping: define bounding boxes as in map in previous script 
#  field work area & sample groups
bbox_fwork <- st_as_sfc(st_bbox(c(xmin = (166.5-0.1), xmax = (167.0+0.1), ymax = (-46.04-0.1), ymin = (-45.52+0.1)), crs = st_crs(4326)))
#  boxes around default value RESERVE.GROUP.LOCATION  
bbox_rgl_fish_biodiv <- get_bbox_anyloc(fish_biodiv_gh) # must use original object, not sf

# for mapping and buffer calculations at correct scale: re-project all sf's to local km  
get_reprojection <- function(sf) st_transform(sf, crs = st_crs("+proj=utm +zone=58G +datum=WGS84 +units=km"))

fish_biodiv_sf_km <- get_reprojection(fish_biodiv_sf)

nzshp_hires_WGS84_sf_km <- get_reprojection(nzshp_hires_WGS84_sf)
nzshp_lores_WGS84_sf_km <- get_reprojection(nzshp_lores_WGS84_sf)

bbox_fwork_km <- get_reprojection(bbox_fwork)

bbox_rgl_fish_biodiv_km <- get_reprojection(bbox_rgl_fish_biodiv)

# calculate 2.5 km buffers
fish_biodiv_sf_km_sid_buff <- fish_biodiv_sf_km %>% select("SET.ID") %>% distinct %>% st_buffer(2.5)

# get dataframes suitable for plotting with below functions - write as function
fish_biodiv_df_edna <- get_plot_df(fish_biodiv_sf_km, "eDNA")
fish_biodiv_df_bruv <- get_plot_df(fish_biodiv_sf_km, "BRUV")
fish_biodiv_df_obis <- get_plot_df(fish_biodiv_sf_km, "OBIS")

# mapping
# --------

# map 1: sampling map from `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_get_OBIS_and_map.r`
map_a <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_get_OBIS_and_map__mapggplot.Rds")

ggsave("210809_998_r_summarize_results_map_main.pdf", plot = map_a, 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 152, height = 121, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# map 2: eDNA observations
map_b <- ggplot() +
      geom_density_2d_filled(data = get_plot_df(fish_biodiv_sf_km, "eDNA"), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
      # facet_grid(. ~ SAMPLE.TYPE) +
      geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
      # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") +
      geom_sf(data = bbox_rgl_fish_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) + 
      # geom_sf_label(data=bbox_rgl_fish_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
      stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "eDNA")}, aes(shape = RESERVE.GROUP), color = "grey20", size = 2) +
      stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "eDNA")}, aes(shape = RESERVE.GROUP), color = "white", size = 1) +
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

ggsave("210809_998_r_summarize_results_map_edna.pdf", plot = map_b, 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 152, height = 121, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# map 2: BRUV observations
map_c <- ggplot() +
      geom_density_2d_filled(data = get_plot_df(fish_biodiv_sf_km, "BRUV"), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
      # facet_grid(. ~ SAMPLE.TYPE) +
      geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
      # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") +
      geom_sf(data = bbox_rgl_fish_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) + 
      # geom_sf_label(data=bbox_rgl_fish_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
      stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "BRUV")}, aes(shape = RESERVE.GROUP), color = "grey20", size = 2) +
      stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "BRUV")}, aes(shape = RESERVE.GROUP), color = "white", size = 1) +
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

ggsave("210809_998_r_summarize_results_map_bruv.pdf", plot = map_c, 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 152, height = 121, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# map 3: local OBIS observations
map_d <- ggplot() +
      geom_density_2d_filled(data = get_plot_df(fish_biodiv_sf_km, "OBIS"), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
      # facet_grid(. ~ SAMPLE.TYPE) +
      geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
      # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") +
      geom_sf(data = bbox_rgl_fish_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) + 
      # geom_sf_label(data=bbox_rgl_fish_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
      stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "OBIS")}, aes(shape = RESERVE.GROUP), color = "grey20", size = 2) +
      stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "OBIS")}, color = "white", size = 1) +
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

ggsave("210809_998_r_summarize_results_map_obis.pdf", plot = map_d, 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 152, height = 121, units = c("mm"),
         dpi = 500, limitsize = TRUE)

save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__mapping.Rdata")

# slooooooooow 
ggarrange( 
  ggarrange(map_a,               ncol = 1, nrow = 1, labels = c("a")),
  ggarrange(map_b, map_c, map_d, ncol = 1, nrow = 3, labels = c("b","c", "d")),
  widths = c(2, 1), ncol = 2, nrow = 1  
  )
  
# save compound plot with better labels then with plot_label = TRUE above
ggsave("210712_998_r_summarize_results__geoheat_edna_bruv_obis.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 152, height = 121, units = c("mm"),
         dpi = 500, limitsize = TRUE)  

save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__mapping_done.Rdata")
         
# VI. Get biodiversity heat map and matching flex table
# ====================================================

# a.) get plotting data sets
# --------------------------

fish_biodiv_tbls <- fish_biodiv |> filter(!(SAMPLE.TYPE %in% c("OBIS") & SET.ID %in% c(1,3,4,5,7,8,9,10,11,12,17,18,19,21,22,23,24,26,27,28,29))) 

# export data fro analysis by MdL
saveRDS(fish_biodiv_tbls, "/Users/paul/Documents/OU_eDNA/201028_Robjects/210703_998_r_summarize_results__data_gtestimate_accumulation_curves.Rds")

  # |>
  # select(SET.ID, SAMPLE.TYPE, RESERVE.GROUP, RESERVE.GROUP.LOCATION) |>
  # distinct()

htmp_tibl_fish <- bind_rows(
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "eDNA",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "eDNA"),
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "BRUV",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "BRUV"), 
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "OBIS",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "OBIS"),
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "PUBL",  tbl = TRUE) %>% add_column(SAMPLE.TYPE = "PUBL"),
)

# b.) add and analyse eDNA BLAST results
# --------------------------------------------

# summary for ms - all species
fish_biodiv |> filter(SAMPLE.TYPE == "eDNA") |> select(ASV, SPECIES) |> distinct(SPECIES)

# summary for ms - all asvs
fish_biodiv |> filter(SAMPLE.TYPE == "eDNA") |> select(ASV, SPECIES) |> distinct(ASV)

# for regression analysis - isolate BLAST results from raw data 
fish_biodiv_blast <- fish_biodiv |> 
  filter(SAMPLE.TYPE == "eDNA") |> 
  select(ASV, RESERVE.GROUP.LOCATION, EDNA.OBS.PRES, FAMILY, GENUS, SPECIES, NCBI.LEVEL, NCBI.TAXDB.INC, NCBI.TAXID, NCBI.TAXID.INC, HSP.GAPS, HSP.IDENTITY.PERC) |>
  arrange(FAMILY, GENUS, SPECIES)

# for regression analysis - add  count locations per species and asv
fish_asv_at_locs <- fish_biodiv_blast |> group_by(SPECIES) |> summarize(LOC.PER.SPC = n_distinct(RESERVE.GROUP.LOCATION)) |> arrange(LOC.PER.SPC) 
fish_asv_at_locs <- fish_biodiv_blast |> select(ASV, EDNA.OBS.PRES, FAMILY, GENUS, SPECIES, RESERVE.GROUP.LOCATION, HSP.GAPS, HSP.IDENTITY.PERC) |>  left_join(fish_asv_at_locs)

# for regression analysis - get a column with non-nz species as per above
fish_asv_at_locs <- fish_asv_at_locs |> mutate(NOT.NZ = as.factor(ifelse( grepl("*", SPECIES, fixed = TRUE), TRUE, FALSE)))

# for regression analysis - re-scale percentages as per MdL - for CIs more easily understandable
fish_asv_at_locs <- fish_asv_at_locs |> mutate(HSP.IDENTITY.PERC = 100 * HSP.IDENTITY.PERC)

# for regression analysis -  save object for external inspection if desirable
# saveRDS({fish_asv_at_locs |> select(LOC.PER.SPC, HSP.GAPS, HSP.IDENTITY.PERC, NOT.NZ)}, "/Users/paul/Documents/OU_eDNA/201028_Robjects/210703_998_r_summarize_results__data_spc_distribution_vs_quality.Rds")
saveRDS(fish_asv_at_locs, "/Users/paul/Documents/OU_eDNA/201028_Robjects/210703_998_r_summarize_results__data_spc_distribution_vs_quality.Rds")

# for manuscript - count fraction of NOT.NZ species among unique eDNA assignments
fish_asv_at_locs |> select(ASV, NOT.NZ) |> distinct() |> group_by(NOT.NZ) |>
  arrange(NOT.NZ) # 92 ASV: 39 False (42.4% True) / 53 True (57% True)

fish_asv_at_locs |> select(SPECIES, NOT.NZ) |> distinct() |> group_by(NOT.NZ) |>
  arrange(NOT.NZ) # 44 SPECIES: 25 False (56.8% True) / 19 True (43.1% True)



# test relationship between  SPC.PER.LOC / HSP.GAPS / HSP.IDENTITY.PERC / NOT.NZ
# ------------------------------------------------------------------------------

# binomial regression
glm_mod <-  glm(NOT.NZ ~ HSP.GAPS + HSP.IDENTITY.PERC, weights = LOC.PER.SPC, family = binomial, data = fish_asv_at_locs)
summary(glm_mod)

# Call:
# glm(formula = NOT.NZ ~ HSP.GAPS + HSP.IDENTITY.PERC, family = binomial, 
#     data = fish_asv_at_locs, weights = LOC.PER.SPC)
# 
# Deviance Residuals: 
#    Min      1Q  Median      3Q     Max  
# -4.959  -1.983   1.008   1.839   2.472  
# 
# Coefficients:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -6.49027    1.91292  -3.393 0.000692 ***
# HSP.GAPS           0.33279    0.08460   3.934 8.36e-05 ***
# HSP.IDENTITY.PERC  0.07174    0.01997   3.592 0.000328 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 572.60  on 141  degrees of freedom
# Residual deviance: 552.31  on 139  degrees of freedom
# AIC: 558.31
# 
# Number of Fisher Scoring iterations: 5


confint(glm_mod, level = 0.95) # probabilities

# Calculating CIS on the scale of the link function and not the response scale
#   https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
#   getting inverse link function for calculating CIS

# unit changes of predictors
fam <- family(glm_mod)
ilink <- fam$linkinv
ilink(confint(glm_mod, level = 0.95))

# Waiting for profiling to be done...
#                          2.5 %     97.5 %
# (Intercept)       3.339469e-05 0.05834373
# HSP.GAPS          5.435871e-01 0.62425910
# HSP.IDENTITY.PERC 5.082516e-01 0.52786043

# check_observation among true/false for figure legend
fish_asv_at_locs |> group_by(NOT.NZ) |> summarise(across(c("SPECIES", "ASV"), list(n_distinct)))

# # A tibble: 2 × 3
#   NOT.NZ SPECIES_1 ASV_1
#   <fct>      <int> <int>
# 1 FALSE         25    39
# 2 TRUE          19    53

coeff_plot <- sjPlot::plot_model(glm_mod, vline.color = "red",  show.values = TRUE) +
  theme_bw() +
  scale_x_discrete(labels = c("Algn. Cov.", "Gaps")) +
  ggtitle("Influence of alignment quality on non-native status")
  
ggsave("210712_998_r_summarize_results__coeff_plot.pdf", plot = coeff_plot, 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 200, height = 135, units = c("mm"),
         dpi = 500, limitsize = TRUE)

model_plot <- plot_model(glm_mod, type = "pred", terms = c("HSP.IDENTITY.PERC", "HSP.GAPS"), show.data = TRUE, jitter = 0.0, ci.lvl = 0.95) +
                theme_bw() +       
                ylab("Probabilty of observing non-native species") +
                xlab("Alignment identity percentage") +
                ggtitle("Predicted probabilities of non-native status") +
                labs(col = "Gaps") +
                theme(legend.position = c(.9, .6),
                      legend.justification = c("right", "top"),
                      legend.box.just = "left",
                      legend.box.background = element_rect(color="grey30", size=0.5)
                      )

ggsave("210712_998_r_summarize_results__asv_bin_regression.pdf", plot = last_plot(), 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 1, width = 200, height = 135, units = c("mm"),
       dpi = 500, limitsize = TRUE)

tab_model(glm_mod)


# eDNA data with BLAST results - other reporting 
# ----------------------------------------------

# for reporting - summaries for gaps and query coverage
fish_biodiv_blast_unq <- fish_biodiv_blast |> distinct(across(c("ASV","FAMILY", "SPECIES","NCBI.LEVEL", "NCBI.TAXDB.INC", "NCBI.TAXID", "NCBI.TAXID.INC", "HSP.GAPS", "HSP.IDENTITY.PERC")))
nrow(fish_biodiv_blast_unq) # 92 ASV resolved to species

fish_biodiv_blast_unq |> filter(HSP.GAPS == 0) |> filter(HSP.IDENTITY.PERC == 1) 

fish_biodiv_blast_unq |> filter(HSP.IDENTITY.PERC != 1)
fish_biodiv_blast_unq |> filter(HSP.IDENTITY.PERC != 0) |> pull("FAMILY") |> unique()

fish_biodiv_blast_unq |> filter(HSP.GAPS != 0)
fish_biodiv_blast_unq |> filter(HSP.GAPS != 0) |> pull("FAMILY") |> unique()

summary(fish_biodiv_blast_unq$HSP.IDENTITY.PERC)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.7868  0.8933  0.9702  0.9321  0.9763  1.0000 
mean(fish_biodiv_blast_unq$HSP.IDENTITY.PERC)    # 0.9320525
sd(fish_biodiv_blast_unq$HSP.IDENTITY.PERC)      # 0.06863656

summary(fish_biodiv_blast_unq$HSP.GAPS)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       0       0       0       1       1      10
mean(fish_biodiv_blast_unq$HSP.GAPS)      # 1
sd(fish_biodiv_blast_unq$HSP.GAPS)      # 1.839732

fish_biodiv_blast_cov <- fish_biodiv_blast |> 
  group_by(SPECIES) |> 
  summarize(BLAST.COV.RNG =  ifelse( signif(100*min(HSP.IDENTITY.PERC), 3) != signif(100*max(HSP.IDENTITY.PERC),3),
                                     paste0( signif(100*min(HSP.IDENTITY.PERC), 3), "-", signif(100*max(HSP.IDENTITY.PERC),3), "%"),
                                     paste0( signif(100*min(HSP.IDENTITY.PERC),3), "%")
                                     ))
  
fish_biodiv_blast_gap <- fish_biodiv_blast |> 
  group_by(SPECIES) |> 
  summarize(BLAST.GAP.RNG =  ifelse( signif( min(HSP.GAPS), 3) != signif(max(HSP.GAPS),3),
                                     paste0( signif(min(HSP.GAPS), 3), "-", signif(max(HSP.GAPS),3)),
                                     paste0( signif(min(HSP.GAPS),3))
                                     ))

# extended data (possibly) for table plot 
htmp_tibl_fish_blrngs <-  htmp_tibl_fish |> left_join(fish_biodiv_blast_cov) |> left_join(fish_biodiv_blast_gap) 


# get matching flex table
# -----------------------

# format data for flex table
tibl_plot <- htmp_tibl_fish_blrngs %>% select(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, TRIVIAL.SPECIES, BLAST.COV.RNG, BLAST.GAP.RNG) %>% 
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
     TRIVIAL.SPECIES = "Common name",
     BLAST.COV.RNG = "Algn. covrg.",
     BLAST.GAP.RNG = "Algn. gaps"
     )) %>%
     fontsize(part = "all", size = 12) %>%
     autofit(add_w = 1, add_h = 0, part = c("all"))

save_as_html(ft, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210712_998_r_summarize_results__spcies_obs_matching_tiles.html")

# order factors in plotting object to match flex table  
y_axis_label_order <- fct_relevel(htmp_tibl_fish$SPECIES, rev(c(tibl_plot$SPECIES)))  # y_axis_label_order <- reorder(htmp_tibl_full$SPECIES, desc(htmp_tibl_full$SPECIES))

# plot out data as tiles
# -----------------------
# - needs vectors ordered with tree
# - needs annotation of non NZ species
# - may need margin sums
plot_htmp <- ggplot(htmp_tibl_fish, aes_string(x = "RESERVE.GROUP.LOCATION", y = y_axis_label_order)) +
    geom_tile(aes(fill = ANY.OBS.PRES) ) +
    scale_fill_gradient(low="black", high="red") +
    geom_text(size = 2, aes(label = ANY.OBS.PRES, color = "white")) +
    facet_grid(. ~ SAMPLE.TYPE, scales = "free") + 
    theme_bw() +
    theme(legend.position = "none", 
          strip.text.y = element_text(angle=0), 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7, face = "italic"), 
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()
          ) +
    xlab("Sampling Locations") # + ylab("Species Observations")

# adjust facet width - https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap/52422707
# convert ggplot object to grob object
grob_htmp <- ggplotGrob(plot_htmp)

# optional: take a look at the grob object's layout
gtable::gtable_show_layout(grob_htmp)

# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.columns <- grob_htmp$layout$l[grepl("panel", grob_htmp$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(plot_htmp)$layout$panel_scales_x,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
grob_htmp$widths[facet.columns] <- grob_htmp$widths[facet.columns] * x.var

# plot result
plot_htmp_adjusted <- as.ggplot(grob_htmp)

ggsave("210712_998_r_summarize_results__biodiv_tiles_only.pdf", plot = plot_htmp_adjusted, 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 155, height = 297, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# Combine heatmap and flextable
# ------------------------------

# get flex table as ggplot object
ft_raster <- as_raster(ft) # webshot and magick.
ft_rgrob  <- rasterGrob(ft_raster)
plot_ft <- ggplot() + theme_void() + annotation_custom(ft_rgrob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

ggarrange(plot_ft, plot_htmp_adjusted, ncol = 2, nrow = 1, labels = c("a","b"), widths = c(1,1))

ggsave("210712_998_r_summarize_results__biodiv_heat.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1.5, width = 210, height = 297, units = c("mm"),
         dpi = 500, limitsize = TRUE)

save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__dis_done.Rdata")

# VI. ANOSIM of observation types and variables
# ===============================================

# subset to local observation without OBIS, OBIS data is too sparse to allow meaningful conclusions
fish_biodiv_local <- fish_biodiv |> filter(SET.ID %!in% c(98, 99)) |> filter(SAMPLE.TYPE %in% c("eDNA", "BRUV", "OBIS"))


# testing function with one data set
get_vegan(distance = "jaccard", tibl = fish_biodiv_local, group_col = "SET.ID", group_row = c("SPECIES"), group_col_ano = "RESERVE.GROUP.LOCATION", obs_methods = "BRUV")


# Analysis for fish (as per eDNA markers) for BRUV and eDNA (data complete across all sets)
# ----------------------------------------------------------------------------------------------

# setting up parameter combinations for complete ANOSIM analysis
anosim_analysis_fish <- expand.grid(
  distance      = c("jaccard"),
  tibl          = c("fish_biodiv_local"),
  group_col     = c("SET.ID"), 
  group_row     = c("SPECIES", "GENUS", "FAMILY", "ORDER"),
  group_col_ano = c("RESERVE.GROUP.LOCATION", "RESERVE.GROUP.INSIDE"), 
  obs_methods   = c("eDNA", "BRUV")
  )

# run ANOSIM analysis
anosim_results_fish <- apply(anosim_analysis_fish, 1, FUN = function(x) try(get_vegan(distance = x[1], tibl = get(x[2]), group_col = x[3], group_row = x[4], group_col_ano = x[5], obs_methods = x[6])))

# inspect results
#  str(anosim_results_fish[[1]])

# get results data frame
anosim_analysis_fish$statistic  <- unlist(lapply(anosim_results_fish, '[[', 5))
anosim_analysis_fish$significance <- unlist(lapply(anosim_results_fish, '[[', 2))
anosim_analysis_fish
anosim_analysis_fish %>% filter(significance <= 0.05 )


# generate flex table - for supplement 
ft_anosim <-  flextable(anosim_analysis_fish) %>% 
   merge_v(j = "distance", target = "distance") %>%
   merge_v(j = "tibl", target = "tibl") %>%
   merge_v(j = "group_col", target = "group_col") %>%
   merge_v(j = "group_row", target = "group_row") %>%
   merge_v(j = "group_col_ano", target = "group_col_ano") %>% 
   merge_v(j = "obs_methods", target = "obs_methods") %>%
   merge_v(j = "significance", target = "significance") %>%
   valign(valign = "top") %>%
   set_header_labels(values = list(distance = "Distance",
     tibl = "Data set",
     group_col = "Replication over",
     group_row = "Tax. level",
     group_col_ano = "Location grouping",
     obs_methods = "Obs. method",
     statistic = "ANOSIM R",
     significance = "Significance")) %>%
     highlight(j = ~ significance, color = function(x) {ifelse(x < 0.05, "lightgreen", NA)} ) %>%
     fontsize(part = "all", size = 9) %>%
     fit_to_width(max_width = 12, inc = 1L, max_iter = 20)
save_as_html(ft_anosim, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210712_998_r_summarize_results__ANOSIM.html")

save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__anaosim_done.Rdata")


# VII. Indicator species analysis for significant ANOMSIM results
# ==============================================================

# testing indicator species analysis - note flag "map" set to TRUE
get_vegan(distance = "jaccard", tibl = fish_biodiv_local, group_col = "SET.ID", group_row = c("SPECIES"), group_col_ano = "RESERVE.GROUP.LOCATION", obs_methods = "BRUV", mp = TRUE)

# Analysis for fish (as per eDNA markers) for BRUV and eDNA (data complete across all sets)
# ----------------------------------------------------------------------------------------------

# setting up parameter combinations for complete ANOSIM analysis
mpatt_analysis_fish <- expand.grid(
  distance      = c("jaccard"),
  tibl          = c("fish_biodiv_local"),
  group_col     = c("SET.ID"), 
  group_row     = c("SPECIES", "GENUS", "FAMILY", "ORDER"),
  group_col_ano = c("RESERVE.GROUP.LOCATION"), 
  obs_methods   = c("BRUV")
  )

# run Indicator species analysis
mpatt_results_fish <- apply(mpatt_analysis_fish, 1, FUN = function(x) try(get_vegan(distance = x[1], tibl = get(x[2]), group_col = x[3], group_row = x[4], group_col_ano = x[5], obs_methods = x[6], mp = TRUE)))

# inspect results
str(mpatt_results_fish)
str(mpatt_results_fish[[1]])
mpatt_results_fish[[1]] # species level - order follows rows in mpatt_analysis_fish
mpatt_results_fish[[2]]
mpatt_results_fish[[3]]
mpatt_results_fish[[4]]

summary(mpatt_results_fish[[1]])

#  Multilevel pattern analysis
#  ---------------------------
# 
#  Association function: r.g
#  Significance level (alpha): 0.05
# 
#  Total number of species: 25
#  Selected number of species: 1 
#  Number of species associated to 1 group: 0 
#  Number of species associated to 2 groups: 1 
#  Number of species associated to 3 groups: 0 
#  Number of species associated to 4 groups: 0 
#  Number of species associated to 5 groups: 0 
# 
#  List of species associated to each combination: 
# 
#  Group WJ CTRL+WJ MR  #sps.  1 
#                        stat p.value  
# Bodianus unimaculatus 0.725  0.0282 *
# ---
# Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1 

summary(mpatt_results_fish[[2]])

# Multilevel pattern analysis
#  ---------------------------
# 
#  Association function: r.g
#  Significance level (alpha): 0.05
# 
#  Total number of species: 23
#  Selected number of species: 1 
#  Number of species associated to 1 group: 0 
#  Number of species associated to 2 groups: 1 
#  Number of species associated to 3 groups: 0 
#  Number of species associated to 4 groups: 0 
#  Number of species associated to 5 groups: 0 
# 
#  List of species associated to each combination: 
# 
#  Group WJ CTRL+WJ MR  #sps.  1 
#           stat p.value  
# Bodianus 0.725    0.03 *
# ---
# Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1 


summary(mpatt_results_fish[[3]])

# Multilevel pattern analysis
#  ---------------------------
# 
#  Association function: r.g
#  Significance level (alpha): 0.05
# 
#  Total number of species: 18
#  Selected number of species: 1 
#  Number of species associated to 1 group: 0 
#  Number of species associated to 2 groups: 0 
#  Number of species associated to 3 groups: 0 
#  Number of species associated to 4 groups: 1 
#  Number of species associated to 5 groups: 0 
# 
#  List of species associated to each combination: 
# 
#  Group FF CTRL+FF MR+WJ CTRL+WJ MR  #sps.  1 
#          stat p.value  
# Labridae 0.67  0.0385 *
# ---
# Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1 

summary(mpatt_results_fish[[4]])

# Multilevel pattern analysis
#  ---------------------------
# 
#  Association function: r.g
#  Significance level (alpha): 0.05
# 
#  Total number of species: 11
#  Selected number of species: 1 
#  Number of species associated to 1 group: 0 
#  Number of species associated to 2 groups: 0 
#  Number of species associated to 3 groups: 1 
#  Number of species associated to 4 groups: 0 
#  Number of species associated to 5 groups: 0 
# 
#  List of species associated to each combination: 
# 
#  Group FF MR+WJ CTRL+WJ MR  #sps.  1 
#              stat p.value  
# Perciformes 0.699  0.0301 *
# ---
# Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1 

save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__multipatt_done.Rdata")

