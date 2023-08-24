#  **********************************
#  *                                *        
#  *  Get results from long tables  *
#  *                                *
#  **********************************

lapply(paste('package:', names(sessionInfo()$otherPkgs),sep=""),detach, character.only=TRUE, unload=TRUE)
rm(list = ls(all.names = TRUE))
gc()
options(tibble.print_max = Inf) 

# Environment setup ----

# _1.) Packages ----

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
library("jtools")		  # model plots

# _2.) API keys ----

# Entrez API key - NCBI taxonomy lookup 
Sys.setenv(ENTREZ_KEY="ecc505b227f772d346fb57816cac0bfda408")
Sys.getenv("ENTREZ_KEY")

# Fishbase API key  - synonyms lookup 
Sys.setenv(FISHBASE_API="sealifebase")
Sys.getenv("FISHBASE_API")

# _3.) Functions  ----

# "not in"

'%!in%' <- function(x,y)!('%in%'(x,y))

# Euler objects

get_euler_object = function(level, tibl){
  require("eulerr")
  require("tidyverse")
  require("magrittr")
  
  # check if needed columns are in the input data
  stopifnot(c("BRUV.OBS.PRES", "EDNA.OBS.PRES", "OBIS.OBS.PRES", "PUBL.OBS.PRES") %in% names(tibl))
  stopifnot(level %in% c("SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"))
  
  # isolate realvant columns for summary
  tibl %<>% select(BRUV.OBS.PRES, EDNA.OBS.PRES, OBIS.OBS.PRES, PUBL.OBS.PRES, SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES) %>% distinct()
  
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

# Euler plots

get_euler_ggplot = function(level, euler_ob, plot_label = TRUE){
  require("tidyverse")
  require ("ggplotify")

  # sanitize input
  stopifnot( class(euler_ob)[1] == "euler")
  stopifnot(level %in% c("SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"))
  
  euler_ggplot <- as.ggplot(
    plot(euler_ob, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels = list(font=1, cex=0.8))
    ) + {if(plot_label == TRUE) labs(subtitle = str_to_sentence(level))}
  
  return(euler_ggplot)
}

# Get spatial points data frames for map generation

#   (from full_biodiv or fish_biodiv )

get_sf_biodiv =  function(tibl){
  require("tidyverse")
  require("magrittr")
  require("sf")
  
  # define columns for mapping add input verification
  cols <- c("SET.ID", "MH.GPS.LAT", "MH.PPS.LONG", "RESERVE.GROUP", "RESERVE.GROUP.INSIDE",
            "RESERVE.GROUP.LOCATION", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS",
            "SPECIES", "ASV", "ABUNDANCE", "SAMPLE.TYPE", "BRUV.OBS.PRES", "EDNA.OBS.PRES", "OBIS.OBS.PRES")
  stopifnot(cols %in% names(tibl))
  
  # select relavant data fro mapping
  tibl %<>% ungroup %>% select(all_of(cols)) %>% arrange(SET.ID) 
  
  # get simple feature df for mapping, define coordinates as WGS84 (degrees)  
  tibl %<>% st_as_sf(coords=c("MH.PPS.LONG","MH.GPS.LAT")) %>% st_set_crs(4326) 

}

# Get bounding box around an area defined by a variable (here default: RESERVE.GROUP.LOCATION)
# - from https://stackoverflow.com/questions/54696440/create-polygons-representing-bounding-boxes-for-subgroups-using-sf

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
  summarise(xmin = min(MH.PPS.LONG) -0.01 ,ymin = min(MH.GPS.LAT) -0.01, xmax=max(MH.PPS.LONG) +0.01, ymax = max(MH.GPS.LAT) +0.01) %>%
  gather(x,lon,c('xmin','xmax')) %>% gather(y,lat,c('ymin','ymax')) %>%
  st_as_sf(coords=c('lon','lat'),crs=4326,remove=F) %>%
  group_by(across(all_of(location))) %>% mutate(angle = calc_angle(lon,lat)) %>%
  arrange(angle) %>% 
  summarise(do_union=FALSE) %>% 
  st_cast('POLYGON')
  
  return(bbox)
}

# Get data frame for heat-map plotting
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

# Get matrix for distance calculations needed for ANOSIM

get_matrix_or_table <- function(tibl, group_col = "RESERVE.GROUP.LOCATION", group_row = "SPECIES", obs_methods = NULL, tbl = FALSE){
  
  # debugging
  #   tibl = get("fish_biodiv")
  #   group_col = "RESERVE.GROUP.LOCATION"
  #   group_row = "SPECIES"
  #   obs_methods = c("eDNA")
  #   tbl = FALSE
  
  stopifnot(group_col %in% names(tibl))
  stopifnot(group_row %in% c("SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"))
  # stopifnot(length(obs_methods) == 1)
  stopifnot(obs_methods %in% c("BRUV", "eDNA", "OBIS", "PUBL", NULL))
  
  require("data.table") # re-use old code rather then finding out how to re implement
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
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "NCBI.TAXID", "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", "TRIVIAL.SPECIES"), .SDcols=c("ANY.OBS.PRES") ]
  } else if (group_row == "GENUS") {
    dtbl_ag <- dtbl[, lapply(.SD, sum, na.rm=TRUE), by=c(group_col, "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS"), .SDcols=c("ANY.OBS.PRES") ]
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

# Get ANOSIM results for a given tibble ("tibl")
# - aggregate observations for "group_col" on taxonomic level "group_row"
# - aggregated observations ar replicates for ANOSIM as levels of "group_col_ano"

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

# Strip asterisk annotations from species strings where needed

get_clean_strings = function(x) { 
  
  x <- stringr::str_replace(x, "\\*+", "")
  x <- stringr::str_replace(x, "\\*+", "")
  x <- stringr::str_replace(x, " $", "")
  
  return(x)
  }

# Logit to probabilities

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# Data read in ----

# check input data of previous script
system("open -a \"Microsoft Excel\" \"/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/999_r_map_and_add_obis__full_data_raw.xlsx\"")

long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/999_r_map_and_add_obis__full_data_raw.Rds")
long_table <- ungroup(long_table)

# Data shaping and formatting ----

# - add trivial names 
# - mark non-NZ species  **(possibly needs to be re-worked)**
# - split "fish" and "full" data
# - filter for data completeness **(possibly needs to be re-worked)**

# _1.) Tidy scientific taxonomy strings ----

# Needed for manuscript and lookup of trivial names 

# erase all that may come after "sp." if there is an "sp." at all
long_table %<>% mutate(SPECIES = gsub("sp\\..*", "sp.", SPECIES))
long_table %<>% mutate(SPECIES = gsub("cf. ", "", SPECIES))
long_table %<>% mutate(SPECIES = str_replace(SPECIES, "\\s\\S*\\s\\S*(.*)", ""))

long_table %<>% mutate(ORDER = ifelse(ORDER == "Squalidae", "Squaliformes", ORDER))
long_table %<>% mutate(ORDER = ifelse(GENUS == "Callanthias", "Perciformes", ORDER))

# _2.) Add trivial names ----

# __a) Fetch  trivial name synonyms from fishbase ----

# isolate species names from analysis data
spc_read <- long_table %>% pull("SPECIES") %>% unique() %>% sort()

# fetch synonyms from Fishbase (20.12.2022)
spc_in_syn <- rfishbase::synonyms(species_list = spc_read)


# __b)  Save state ----

# save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__syn_lookup.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__syn_lookup.Rdata")

# __c) Inspect trivial name synonyms from fishbase ----

# inspect synonyms
spc_in_syn %>% select(synonym, Species, Status)

# isolate fishbase data that needs to be used to corrected in long table
fb_data_for_correction <- spc_in_syn %>% filter(Status == "synonym") %>% filter(synonym != "Opistognathus sp.") 

# inspect those data again
fb_data_for_correction %>% select(synonym, Species, Status)

# correct synonymous species and genus names in `long_table`, that can be named properly
# Monocentris japonicus     Monocentris japonica synonym
# Helicolenus hilgendorfi Helicolenus hilgendorfii synonym
# Cephaloscyllium isabellum Cephaloscyllium isabella synonym
# Cheilodactylus zonatus       Goniistius zonatus synonym

# confirm that synonyms are in long table
long_table %>% filter(SPECIES %in% c("Monocentris japonicus", "Helicolenus hilgendorfi", 
                                     "Cephaloscyllium isabellum", "Cheilodactylus zonatus")) %>% dplyr::select(FAMILY, GENUS, SPECIES)


# writing correct scientific SPECIES names into table overwriting synonymous names
long_table %<>%  mutate(SPECIES = case_when(SPECIES == "Monocentris japonicus"     ~ "Monocentris japonica",
                                           SPECIES == "Helicolenus hilgendorfi"   ~ "Helicolenus hilgendorfii",
                                           SPECIES == "Cephaloscyllium isabellum" ~ "Cephaloscyllium isabella",
                                           SPECIES == "Cheilodactylus zonatus"    ~ "Goniistius zonatus",
                                                                             TRUE ~ SPECIES))
# writing correct scientific GENUS names into table overwriting synonymous names
long_table %<>%  mutate(GENUS = case_when(SPECIES == "Monocentris japonica" ~ "Monocentris",
                                          SPECIES == "Helicolenus hilgendorfii" ~ "Helicolenus",
                                          SPECIES == "Cephaloscyllium isabella" ~ "Cephaloscyllium",
                                          SPECIES == "Goniistius zonatus" ~ "Goniistius",
                                            TRUE ~ GENUS))

# confirm that synonyms have been replaced - leaving one old species in query  - shoul only return one value
long_table %>% filter(SPECIES %in% c("Monocentris japonica", "Helicolenus hilgendorfii"))

# __d) Add trivial names to main object  ---- 

# get trivial names
species_vector <- ungroup(long_table) |> select(SPECIES) |> distinct() |> pull(SPECIES)
trivial_obj <- sci2comm(species_vector, db = "ncbi", simplify = TRUE) # <- running also without API key

# format trivial names for left join
trivial_df <- data.frame(do.call(rbind , trivial_obj))
trivial_df <- rownames_to_column(trivial_df, var = "SPECIES")
trivial_df <- trivial_df |>  setNames( c("SPECIES", "TRIVIAL.SPECIES")) |> as_tibble()

# add trivial names to object 
long_table %<>% left_join(trivial_df)

# __e) Save state, but currently loading from here ----

# save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__got_trivial-names.Rdata")
load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__got_trivial-names.Rdata")

# __f) Add more trivial names by manual look up ---- 

# check initial state with of trivial species names:
long_table %>% dplyr::select(SPECIES, TRIVIAL.SPECIES) %>% distinct() %>% arrange(SPECIES) %>% print(n = Inf)

# reworked on 23-Jan-2023 and 6 Mar 2023 - looked up trivial names that could not be found automatically
long_table %<>% mutate(TRIVIAL.SPECIES = 
                        case_when(SPECIES == "Acanthoclinus fuscus" ~ "olive rockfish",
                                  SPECIES == "Acanthoclinus littoreus" ~ "New Zealand rockfish",                             
                                  SPECIES == "Acanthoclinus marilynae" ~ "Stout rockfish",
                                  SPECIES == "Acanthoclinus matti" ~ "New Zealand longfin",                             
                                  SPECIES == "Acanthoclinus rua" ~ "little rockfish",
                                  SPECIES == "Aplidium adamsi" ~ "[squirt]",
                                  SPECIES == "Aplidium coronum" ~ "[squirt]",
                                  SPECIES == "Aplidium phortax" ~ "[squirt]",
                                  SPECIES == "Aplidium powelli" ~ "[squirt]",
                                  SPECIES == "Aplodactylus arctidens" ~ "marblefish",
                                  SPECIES == "Arctocephalus australis" ~ "South American fur seal", 
                                  SPECIES == "Botrylloides leachii" ~ "[tunicate]",
                                  SPECIES == "Botryllus stewartensis" ~ "[tunicate]",                         
                                  SPECIES == "Bovichtus angustifrons" ~ " horny thornfish",
                                  SPECIES == "Caesioperca lepidoptera" ~ "butterfly perch",
                                  SPECIES == "Callanthias allporti" ~ "splendid sea perch",
                                  SPECIES == "Callanthias japonicus" ~ "yellowsail red bass",
                                  SPECIES == "Cephaloscyllium isabella" ~ "draughtsboard shark",
                                  SPECIES == "Chaetodon zanzibarensis" ~ "Zanzibar butterflyfish",                             
                                  SPECIES == "Cheilodactylus variegatus" ~ "Peruvian morwong",
                                  SPECIES == "Cnemidocarpa bicornuta" ~ "[tunicate]",                             
                                  SPECIES == "Cnemidocarpa nisiotis" ~ "[tunicate]",                             
                                  SPECIES == "Cominella sp." ~ "[snail]",                             
                                  SPECIES == "Cryptichthys jojettae" ~ "cryptic triplefin",                             
                                  SPECIES == "Didemnum inveteratum" ~ "[tunicate]",                             
                                  SPECIES == "Diplosoma listerianum" ~ "[tunicate]",                     
                                  SPECIES == "Diplosoma velatum" ~ "[tunicate]",                            
                                  SPECIES == "Eptatretus cirrahtus" ~ "broadgilled hagfish",                             
                                  SPECIES == "Eudistoma circumvallatum" ~ "[tunicate]",                            
                                  SPECIES == "Fiordichthys slartibartfasti" ~ "Fiordland brotula",                             
                                  SPECIES == "Forsterygion malcolmi" ~ "mottled triplefin",                          
                                  SPECIES == "Forsterygion maryannae" ~ "oblique-swimming triplefin",                             
                                  SPECIES == "Gaidropsarus novaezelandi" ~ "New Zealand rockling",                             
                                  SPECIES == "Galaxias argenteus" ~ "giant kōkopu",                            
                                  SPECIES == "Galaxias eldoni" ~ "Eldon\'s galaxias",                             
                                  SPECIES == "Gobiopsis atrata" ~ "New Zealand black goby",                             
                                  SPECIES == "Gymnoscopelus nicholsi" ~ "Nichol's lanternfish",                             
                                  SPECIES == "Helcogramma striata" ~ "tropical striped triplefin",                           
                                  SPECIES == "Helicolenus hilgendorfii" ~ "Hilgendorf's saucord",                             
                                  SPECIES == "Helicolenus percoides" ~ "red gurnard perch",                             
                                  SPECIES == "Hemerocoetes monopterygius" ~ "opalfish",                             
                                  SPECIES == "Hygophum hygomii" ~ "Bermuda lantern fish",                             
                                  SPECIES == "Hypoplectrodes huntii" ~ "redbanded perch",                             
                                  SPECIES == "Karalepis stewarti" ~ "scaly-headed triplefin",                             
                                  SPECIES == "Lepidoperca tasmanica" ~ "Tasmanian perch",                             
                                  SPECIES == "Lissocampus filum" ~ "shortsnout pipefish",                             
                                  SPECIES == "Lissoclinum notti" ~ "[tunicate]",                          
                                  SPECIES == "Lotella phycis" ~ "Beardie",
                                  SPECIES == "Maurolicus muelleri" ~ "pennant pearlside",
                                  SPECIES == "Mendosoma lineatum" ~ "telescope fish",                             
                                  SPECIES == "Modicus minimus" ~ "small clingfish",                             
                                  SPECIES == "Modicus tangaroa" ~ "eyespot clingfish",                             
                                  SPECIES == "Monocentris japonica" ~ "Japanese pineapplefish",                             
                                  SPECIES == "Morus serrator" ~ "Australasian gannet - [bird]",                             
                                  SPECIES == "Mustelus asterias" ~ "starry smooth-hound",                             
                                  SPECIES == "Notoclinops caerulepunctus" ~ "blue dot triplefin",                             
                                  SPECIES == "Notoclinops segmentatus" ~ "blue-eyed triplefin",                             
                                  SPECIES == "Notoclinus compressus" ~ "Brown topknot",                             
                                  SPECIES == "Notoclinus fenestratus" ~ "New Zealand topknot",                             
                                  SPECIES == "Notolabrus cinctus" ~ "girdled wrasse",                             
                                  SPECIES == "Opistognathus iyonis" ~ "well-building jawfish",                             
                                  SPECIES == "Opistognathus sp." ~ "jawfish",                             
                                  SPECIES == "Parapercis decemfasciata" ~ NA,                            
                                  SPECIES == "Patiriella regularis" ~ "New Zealand common cushion star [sea star]",                             
                                  SPECIES == "Polyprion oxygeneios" ~ "hāpuku",                             
                                  SPECIES == "Pseudolabrus miles" ~ "Scarlet wrasse",                             
                                  SPECIES == "Ritterella sigillinoides" ~ "[tunicate]",                             
                                  SPECIES == "Ruanoho decemdigitatus" ~ "longfinned triplefin",                             
                                  SPECIES == "Scorpaena papillosa" ~ "red scorpionfish",                             
                                  SPECIES == "Synoicum kuranui" ~ "[tunicate]",                             
                                  SPECIES == "Synoicum occidentalis" ~ "[tunicate]",                             
                                  SPECIES == "Synoicum stewartense" ~ "[tunicate]",                              
                                  SPECIES == "Thalasseleotris iota" ~ "New Zealand pygmy sleeper",                             
                                  SPECIES == "Trididemnum shawi" ~ "[tunicate]",
                                  TRUE ~ TRIVIAL.SPECIES))

# check altered state with of trivial species names:
long_table %>% dplyr::select(SPECIES, TRIVIAL.SPECIES) %>% distinct() %>% arrange(SPECIES) %>% print(n = Inf)

# __g) Correcting spelling mistakes  ----

long_table %<>% mutate(SPECIES = 
                         case_when(SPECIES %in% "Eptatretus cirrahtus"  ~ "Eptatretus cirrhatus",
                                   TRUE ~ SPECIES)) #  %>% select(SPECIES) %>% distinct() %>% arrange(SPECIES)


# __h) Save state ----

# save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__got_more_trivial-names.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__got_more_trivial-names.Rdata")

# _3.) Mark non-fish and non-New Zealand taxa  ----

# __a) Show current state of taxonomy strings ----

long_table %>% dplyr::select(GENUS, SPECIES, TRIVIAL.SPECIES) %>% distinct() %>% arrange(SPECIES) %>% print(n = Inf) 

# __b) Export data for manual inspection ----

# export taxonomy strings for manual lookup  
long_table %>% 
  dplyr::select(CLASS, ORDER, FAMILY, GENUS, SPECIES, TRIVIAL.SPECIES) %>% distinct() %>% 
  arrange(SPECIES) %>% 
  print(n = Inf) %>% 
  writexl::write_xlsx(., "/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part.xlsx")

# add here from `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part_annotated.xlsx`
# marked non-NZ species  
# - originally using list: Roberts, C., Stewart, A., Struthers, C., Barker, J. & Kortet, S. 2019 Checklist of the Fishes of New Zealand. 
# - using list: Checklist of the Fishes of New Zealand: version 1.2 July 2020 CD Roberts, AL Stewart, CD Struthers, JJ Barker and S Kortet Museum of New Zealand Te Papa Tongarewa
#         at : /Users/paul/Documents/OU_eDNA/200224_references/210908_MA_DOC001887_TePapa_Checklist-of-Fishes-of_full.pdf


# __c) Define non-NZ fish taxa after manual inspection ----

# fish genera not known from NZ waters (found manually by literature search)
# - see `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part_annotated.xlsx`

nonnz_fish <- c(
  "Alburnus alburnus", "Aplocheilus lineatus",  "Asterropteryx semipunctata",  "Atherinomorus lacunosus",  "Bostrychus zonatus",  "Bovichtus angustifrons",  "Callanthias japonicus",  "Caprodon schlegelii",  "Chaetodon zanzibarensis",  "Cheilodactylus variegatus",  "Chelidonichthys spinosus",  "Conodon nobilis",  "Engraulis japonicus",  "Erythrocles schlegelii",  "Gaidropsarus argentatus",  "Gobiesox maeandricus",  "Goniistius zonatus",  "Gymnoscopelus nicholsi",  "Helcogramma striata",  "Helicolenus hilgendorfii",  "Helicolenus percoides",  "Opistognathus iyonis",  "Opistognathus liturus",  "Opistognathus punctatus",  "Opistognathus sp.",  "Parapercis decemfasciata",  "Pseudophycis barbata",  "Scobinichthys granulatus",  "Scomber japonicus",  "Trachurus japonicus",  "Mustelus asterias",  "Squalus suckleyi"
  ) 

# __d) Mark non-fish and non-NZ status among taxonomy strings ----

#   16-Mar-2021 add asterisks ("*") to non-NZ species, and ("**") to non-fish (mammals and crustaceans)
#   see `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part_annotated.xlsx` for an annoated list

long_table %<>% mutate(SPECIES = 
                         case_when(SPECIES %in% nonnz_fish                                  ~ paste0(SPECIES, "*"),
                                   CLASS %!in% c("Actinopteri", "Chondrichthyes", "Myxini") ~ paste0(SPECIES, "**"),
                                                                                              TRUE ~ SPECIES)
                       )

long_table %<>% mutate(GENUS = 
                         case_when(SPECIES %in% nonnz_fish                                   ~ paste0(GENUS, "*"),
                                   CLASS %!in% c("Actinopteri", "Chondrichthyes", "Myxini")  ~ paste0(GENUS, "**"),
                                                                                               TRUE ~ GENUS)
                       )


# __h) Save state ----

# save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__taxa_marked.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__taxa_marked.Rdata")

# _4.) Encode and Evaluate MEGAN species assignments ----

# __a) Get a list of spcies detected by MEGAN ----

# get the following vector from edited file /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_new_analysis_outputs/751_12S_single_end_ee3-seq_blast-noenv-ex__taxon_to_count_edited.txt

# manually created
megan_species <- c("Notorynchus cepedianus", "Anguilla australis",  "Trachurus japonicus", "Peltorhamphus novaezeelandiae",
                   "Oplegnathus fasciatus", "Callanthias japonicus", "Asterropteryx semipunctata", "Entomacrodus stellifer",
                   "Poecilia reticulata", "Macruronus novaezelandiae", "Oncorhynchus mykiss", "Lampanyctodes hectoris", 
                   "Caprodon schlegelii", "Engraulis japonicus", "Cervus elaphus", "Tursiops truncatus", 
                   "Pontoscolex corethrurus", "Lamellibrachia barhami", "Lamellibrachia satsuma", "Canis lupus", 
                   "Arctocephalus forsteri", "Katsuwonus pelamis", "Mus musculus", "Homo sapiens", "Pan troglodytes",
                   "Bos taurus", "Ovis aries", "Sus scrofa", "Porcellio scaber", "Brontispa longissima", "Odax pullus")

# manually created
megan_fish_species <- c("Anguilla australis", "Asterropteryx semipunctata", "Callanthias japonicus", "Caprodon schlegelii", "Engraulis japonicus", "Katsuwonus pelamis", 
                        "Macruronus novaezelandiae", "Notorynchus cepedianus", "Oncorhynchus mykiss", "Trachurus japonicus" )

# __b) Show eDNA species detected by MEGAN ----

long_table %>% filter(SAMPLE.TYPE == "eDNA") %>% filter(get_clean_strings(SPECIES) %in% megan_species) %>% select(SPECIES) %>% distinct()

# __c) Show eDNA species not detected by MEGAN ----

long_table %>% filter(SAMPLE.TYPE == "eDNA") %>% filter(get_clean_strings(SPECIES) %!in% megan_species) %>% select(SPECIES) %>% distinct()
                        
# __d) Show eDNA fish species detected by MEGAN data filtered for fish  ----

long_table %>% filter(SAMPLE.TYPE == "eDNA") %>% filter(CLASS %in% c("Actinopteri", "Chondrichthyes", "Myxini")) %>%  filter(get_clean_strings(SPECIES) %in% megan_fish_species) %>% select(SPECIES) %>% distinct()

# __e) Show eDNA fish species not detected by MEGAN data filtered for fish  ----

long_table %>% filter(SAMPLE.TYPE == "eDNA") %>% filter(CLASS %in% c("Actinopteri", "Chondrichthyes", "Myxini")) %>%  filter(get_clean_strings(SPECIES) %!in% megan_fish_species) %>% select(SPECIES) %>% distinct()

# __f) Mark MEGAN-detected species in long table ----

long_table %<>% mutate(SPECIES = case_when(get_clean_strings(SPECIES) %in% megan_species  ~ paste0(SPECIES, " ***"), TRUE ~ SPECIES))

# __g) Save state ----

# save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__finished_megan_integration.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__finished_megan_integration.Rdata")

# _6.) Get equivalent of `BOTH.PRES` ----

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

# _7.) Split data into sets for "fish", "other", and "full" data ----

# __a) Split data ----

long_table %>% distinct()
fish_biodiv <- long_table %>% distinct() %>% filter(CLASS %in% c("Actinopteri", "Chondrichthyes", "Myxini")) %>% filter(!(GENUS %in% c("Sardinops")))

# __b) Save state ----

# save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__data_filtered.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__data_filtered.Rdata")

# __c) Export data fot analysis by MdL ----

fish_biodiv_tbls <- fish_biodiv |> filter(!(SAMPLE.TYPE %in% c("OBIS") & SET.ID %in% c(1,3,4,5,7,8,9,10,11,12,17,18,19,21,22,23,24,26,27,28,29))) 
saveRDS(fish_biodiv_tbls, "/Users/paul/Documents/OU_eDNA/201028_Robjects/230515_999_r_summarize_results__data_gtestimate_accumulation_curves.Rds")

# Get and Export Plain Flextable ----

# _1.) Format data for flex table ----

obs_sums <- fish_biodiv %>% ungroup() %>% group_by(SPECIES) %>%  
  select(SPECIES, BRUV.OBS.PRES, EDNA.OBS.PRES, OBIS.OBS.PRES, PUBL.OBS.PRES) %>%   
  summarize(BRUV.OBS.PRES.SUM = ifelse(sum(BRUV.OBS.PRES), sum(BRUV.OBS.PRES), NA),
            EDNA.OBS.PRES.SUM = ifelse(sum(EDNA.OBS.PRES), sum(EDNA.OBS.PRES), NA),
            OBIS.OBS.PRES.SUM = ifelse(sum(OBIS.OBS.PRES), TRUE, NA), 
            PUBL.OBS.PRES.SUM = ifelse(sum(PUBL.OBS.PRES), TRUE, NA)
            )
 
spcies_obs_sums <- ungroup(fish_biodiv) %>% select(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, TRIVIAL.SPECIES) %>% distinct() %>%
  left_join(obs_sums) %>% arrange(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES) 

# _2.) Generate Flextable ----

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

# _3.) Save Flextable as .docx file ----

save_as_docx(ft_spcies_obs_sums, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230512_999_r_summarize_results__all_data.docx",
             pr_section = prop_section(
               page_size = page_size(orient = "portrait"), type = "continuous"
             ))

# _4.) Save Flextable as .html file ----

save_as_html(ft_spcies_obs_sums, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230512_998_r_summarize_results__all_data.html")

# Get Figure 1 - geographical maps with heat overlays (Figure 1)----

# compare script  ~/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r

# _1.) Data preparation ----

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
get_reprojection <- function(sf) st_transform(sf, crs = st_crs(2193))

fish_biodiv_sf_km <- get_reprojection(fish_biodiv_sf)

nzshp_hires_WGS84_sf_km <- get_reprojection(nzshp_hires_WGS84_sf)
nzshp_lores_WGS84_sf_km <- get_reprojection(nzshp_lores_WGS84_sf)

bbox_fwork_km <- get_reprojection(bbox_fwork)

bbox_rgl_fish_biodiv_km <- get_reprojection(bbox_rgl_fish_biodiv)

# calculate 2.5 km buffers
fish_biodiv_sf_km_sid_buff <- fish_biodiv_sf_km %>% select("SET.ID") %>% distinct %>% st_buffer(2.5)

# get data frames suitable for plotting with below functions - write as function
fish_biodiv_df_edna <- get_plot_df(fish_biodiv_sf_km, "eDNA")
fish_biodiv_df_bruv <- get_plot_df(fish_biodiv_sf_km, "BRUV")
fish_biodiv_df_obis <- get_plot_df(fish_biodiv_sf_km, "OBIS")

# mapping

# _2.) Main map ----

# map 1: sampling map from `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_get_OBIS_and_map.r`
# map_a <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_get_OBIS_and_map__mapggplot.Rds")

map_a <- readRDS(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/999_r_get_OBIS_and_map__mapggplot.Rds")

ggsave("230515_999_r_summarize_results_map_main.pdf", plot = map_a, 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 1, width = 152, height = 121, units = c("mm"),
       dpi = 500, limitsize = TRUE)

# _3.) Mini maps ----

# map 2: eDNA observations
map_b <- ggplot() +
  geom_density_2d_filled(data = get_plot_df(fish_biodiv_sf_km, "eDNA"), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
  # facet_grid(. ~ SAMPLE.TYPE) +
  geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
  # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") +
  geom_sf(data = bbox_rgl_fish_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) + 
  # geom_sf_label(data=bbox_rgl_fish_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
  stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "eDNA")}, aes(shape = RESERVE.GROUP), color = "grey20", size = 3) +
  stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "eDNA")}, aes(shape = RESERVE.GROUP), color = "white", size = 2) +
  coord_sf(xlim = c((1125643-40000), (1125643+4000)), ylim = c((4909254-30000),(4909254+31600)), expand = FALSE) +
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

ggsave("230515_999_r_summarize_results_map_edna.pdf", plot = map_b, 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 0.5, width = 152, height = 121, units = c("mm"),
       dpi = 500, limitsize = TRUE)

# map 2: BRUV observations
map_c <- ggplot() +
  geom_density_2d_filled(data = get_plot_df(fish_biodiv_sf_km, "BRUV"), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
  # facet_grid(. ~ SAMPLE.TYPE) +
  geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
  # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") +
  geom_sf(data = bbox_rgl_fish_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) + 
  # geom_sf_label(data=bbox_rgl_fish_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
  stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "BRUV")}, aes(shape = RESERVE.GROUP), color = "grey20", size = 3) +
  stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "BRUV")}, aes(shape = RESERVE.GROUP), color = "white", size = 2) +
  coord_sf(xlim = c((1125643-40000), (1125643+4000)), ylim = c((4909254-30000),(4909254+31600)), expand = FALSE) +
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

ggsave("230515_999_r_summarize_results_map_bruv.pdf", plot = map_c, 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 0.5, width = 152, height = 121, units = c("mm"),
       dpi = 500, limitsize = TRUE)

# map 3: local OBIS observations
map_d <- ggplot() +
  geom_density_2d_filled(data = get_plot_df(fish_biodiv_sf_km, "OBIS"), aes(x= lon , y = lat), contour_var = "count", alpha = 0.5) +
  # facet_grid(. ~ SAMPLE.TYPE) +
  geom_sf(data = nzshp_lores_WGS84_sf_km, color=alpha("grey20",1), alpha = 0.8) +
  # geom_sf(data = fish_biodiv_sf_km_sid_buff, fill = NA, colour = "darkgrey") +
  geom_sf(data = bbox_rgl_fish_biodiv_km, fill = NA, colour = "grey20", linetype = "dotted", size = 0.5) + 
  # geom_sf_label(data=bbox_rgl_fish_biodiv_km, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 7, nudge_y = 6.5) +
  stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "OBIS")}, aes(shape = RESERVE.GROUP), color = "grey20", size = 3) +
  stat_sf_coordinates(data = {fish_biodiv_sf_km |> filter(SAMPLE.TYPE == "OBIS")}, aes(shape = RESERVE.GROUP), color = "white", size = 2) +
  coord_sf(xlim = c((1125643-40000), (1125643-800)), ylim = c((4909254-30000),(4909254+30750)), expand = FALSE) +
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

ggsave("230515_999_r_summarize_results_map_obis.pdf", plot = map_d, 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = .5, width = 152, height = 121, units = c("mm"),
       dpi = 500, limitsize = TRUE)

# _4.) Save state ----

# save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping.Rdata")

# _5.) Arrange and save compound map ----

# slooooooooow
ggarrange(
  ggarrange(
    map_a,
    ncol = 1,
    nrow = 1,
    labels = c("a")
  ),
  ggarrange(
    map_b,
    map_c,
    map_d,
    ncol = 1,
    nrow = 3,
    labels = c("b", "c", "d"),
    align = "hv"
  ),
  widths = c(2, 1),
  ncol = 2,
  nrow = 1
)

# save compound plot with better labels then with plot_label = TRUE above
ggsave("230515_999_r_summarize_results__geoheat_edna_bruv_obis.pdf", plot = last_plot(), 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 1, width = 152, height = 121, units = c("mm"),
       dpi = 500, limitsize = TRUE)  

# _6.)  Save state ----

# save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping_done.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping_done.Rdata")

# Results for main text (1 of 2)  ----

# _0.) Results - Data collation ----

# __a) Count Sites that yielded BRUV or eDNA data ----
fish_biodiv |> 
  filter(SET.ID %!in% c(98, 99)) |>
  filter(SAMPLE.TYPE %in% c("BRUV", "EDNA")) |> 
  select(SET.ID, SPECIES, SAMPLE.TYPE, RESERVE.GROUP.LOCATION) |> 
  distinct() |> pull(SET.ID) |> unique()  |> length() 

# __b) Show areas that yielded OBIS data ---- 
fish_biodiv |> 
  filter(SET.ID %!in% c(98, 99)) |>
  filter(SAMPLE.TYPE == "OBIS") |> 
  select(SPECIES, SAMPLE.TYPE, RESERVE.GROUP.LOCATION) |> 
  distinct() |> pull(RESERVE.GROUP.LOCATION) |> unique()
# "LS CTRL" "FF MR"   "FF CTRL" "WJ MR"   "WJ CTRL"

# __c) Count Sites that yielded OBIS data ----
fish_biodiv |> 
  filter(SET.ID %!in% c(98, 99)) |>
  filter(SAMPLE.TYPE == "OBIS") |> 
  select(SET.ID, SPECIES, SAMPLE.TYPE, RESERVE.GROUP.LOCATION) |> 
  distinct() |> pull(SET.ID) |> unique()  |> length()  
# 9

# _1.) Results - Obtaining a fish species list of southern Fjordland ----

# __a) What does the list look like? ----

# ____ General species counts ----

nrow(spcies_obs_sums)                                    # found 116 species across all data sets
nrow(spcies_obs_sums |> filter (CLASS == "Actinopteri")) #       105 Actinopteri
nrow(spcies_obs_sums |> filter (CLASS == "Chondrichthyes")) #     10 Chondrichthyes
nrow(spcies_obs_sums |> filter (CLASS == "Myxini")) #              1 Myxini
 
# ____ General species counts for each data type ----

nrow(spcies_obs_sums |> filter (!is.na(PUBL.OBS.PRES.SUM)))  # 59 -> 61 PUBL (Fiordland)
nrow(spcies_obs_sums |> filter (!is.na(OBIS.OBS.PRES.SUM)))  # 25 -> 28 OBIS (in circle)
nrow(spcies_obs_sums |> filter (!is.na(BRUV.OBS.PRES.SUM)))  # 25 -> 26 BRUV (in study area)
nrow(spcies_obs_sums |> filter (!is.na(EDNA.OBS.PRES.SUM)))  # 44 -> 43 EDNA (in study area)

# show the 43 eDNA species - for LCA see way above in script at section `Encode and Evaluate MEGAN species assignment`
spcies_obs_sums |> filter (!is.na(EDNA.OBS.PRES.SUM)) %>% select(-c(PHYLUM, CLASS, ORDER, FAMILY, GENUS))

# ____ Augmentation of  OBIS and literature by eDNA and BRUV ----


# ____ - Get fish species in OBIS and literature ----

pbob_species <- fish_biodiv |> select(SPECIES, SAMPLE.TYPE) |> filter(SAMPLE.TYPE %in% c("OBIS", "PUBL")) |> distinct()
unique(pbob_species[["SPECIES"]]) # Species in OBIS and literature - 70 

# ____ - Get fish species in eDNA  ----

edna_species <- fish_biodiv |> select(SPECIES, SAMPLE.TYPE) |> filter(SAMPLE.TYPE == "eDNA") |> distinct()
edna_species |> print(n = Inf)

# ____ - Get BRUV fish species ----

bruv_species <- fish_biodiv |> select(SPECIES, SAMPLE.TYPE) |> filter(SAMPLE.TYPE == "BRUV") |> distinct()
bruv_species |> print(n = Inf)

# ____ - Intersect BRUV and literture ----

intersect(bruv_species$SPECIES, pbob_species$SPECIES)

# ____ - Intersect eDNA and literature ----

intersect(edna_species$SPECIES, pbob_species$SPECIES)


# __b) Interesting fish hits ----

# ____ - Notable fish species detected only using eDNA - but see main text and table 

setdiff(
  # eDNA with perfect alignments
  long_table %>% filter(HSP.GAPS == 0 & HSP.IDENTITY.PERC == 1) %>% filter(CLASS %in% c("Actinopteri", "Chondrichthyes", "Myxini")) %>% pull(SPECIES) %>% unique(),
  
  # all others
  union(bruv_species$SPECIES, pbob_species$SPECIES)
)

# ____ - Notable fish species detected only using BRUV - but see main text and table 

setdiff(
  # eDNA with perfect alignments
  bruv_species$SPECIES,
  
  # all others
  union(edna_species$SPECIES, pbob_species$SPECIES)
)

# ____ - Notable fish species detected using BRUV and eDNA - but see main text and table 

intersect(bruv_species$SPECIES, edna_species$SPECIES)

# __c)  Interesting non-fish hits ----

# For eDNA and LCA see MEGAN file `/Users/paul/Documents/OU_eDNA/201126_preprocessing/megan/751_12S_single_end_ee3-seq_blast-noenv.rma6`

# For BRUV previous script `/Users/paul/Documents/OU_eDNA/200901_scripts/996_r_get_BRUV_long_table.r`


# _2.) Results - Evaluating the usefulness of eDNA survey without local reference data ----

# __a) Comparison to LCA ----

# See section above `Encode and Evaluate MEGAN species assignments`

# All MEAGN assignmnets not New Zealand
setdiff(
  
edna_species$SPECIES[grep("\\* \\*\\*\\*",  edna_species$SPECIES)], 

# all others
union(bruv_species$SPECIES, pbob_species$SPECIES)
)


# All eDNA assignmnets not New Zealand

setdiff(
  # eDNA with perfect alignments
  long_table %>% filter(CLASS %in% c("Actinopteri", "Chondrichthyes", "Myxini")) %>% pull(SPECIES) %>% unique(),
  
  # all others
  union(bruv_species$SPECIES, pbob_species$SPECIES)
)

# __b) Summary of alignment qualities  (Supplement)----

# Among eDNA data there are 43 distinct species observations
fish_biodiv |> select(ASV, SAMPLE.TYPE, SPECIES) |> filter(SAMPLE.TYPE == "eDNA") |> distinct(SPECIES)

# Among eDNA data there are 156 non-distinct species observations
fish_biodiv |> select(ASV, SAMPLE.TYPE, SPECIES) |> filter(SAMPLE.TYPE == "eDNA") 

# There are 96 eDNA distcict observations
fish_biodiv_blast_unq <- fish_biodiv_blast |> distinct(across(c("ASV","FAMILY", "SPECIES","NCBI.LEVEL", "NCBI.TAXDB.INC", "NCBI.TAXID", "NCBI.TAXID.INC", contains("HSP"))))
nrow(fish_biodiv_blast_unq) 

# ASV with 0 gaps and full coverage - "eight yielded flawless alignments" 
fish_biodiv_blast_unq |> filter(HSP.GAPS == 0) |> filter(HSP.IDENTITY.PERC == 1) 

# Eighty-eight ASVs had variable query coverage, and families 
fish_biodiv_blast_unq |> filter(HSP.IDENTITY.PERC != 1)
fish_biodiv_blast_unq |> filter(HSP.IDENTITY.PERC != 0) |> pull("FAMILY") |> unique()

# 33 ASVs had variable gap counts 
fish_biodiv_blast_unq |> filter(HSP.GAPS != 0)
fish_biodiv_blast_unq |> filter(HSP.GAPS != 0) |> pull("FAMILY") |> unique()

# find worst matching families - not in manuscript
fish_biodiv_blast_unq |> group_by(FAMILY) |> slice(which.max(HSP.GAPS)) |> slice(which.min(HSP.IDENTITY.PERC)) |> arrange(desc(HSP.GAPS))

# Query coverage statistsics
summary(fish_biodiv_blast_unq$HSP.IDENTITY.PERC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7868  0.9037  0.9702  0.9379  0.9827  1.0000
sd(fish_biodiv_blast_unq$HSP.IDENTITY.PERC)    # 0.06043447

# Gap count statistix
summary(fish_biodiv_blast_unq$HSP.GAPS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.9792  1.0000 10.0000
sd(fish_biodiv_blast_unq$HSP.GAPS)    # 1.800463

# __c) Regression analysis ----

# ____ Format data for regression analysis ----

# summary for ms - all asvs
ms_all_asvs <- fish_biodiv |> filter(SAMPLE.TYPE == "eDNA") |> select(ASV, SPECIES) |> distinct(ASV)

# for regression analysis - isolate BLAST results from raw data 
fish_biodiv_blast <- fish_biodiv |> 
  filter(SAMPLE.TYPE == "eDNA") |> 
  select(ASV, RESERVE.GROUP.LOCATION, EDNA.OBS.PRES, FAMILY, GENUS, SPECIES, NCBI.LEVEL, NCBI.TAXDB.INC, NCBI.TAXID, NCBI.TAXID.INC, HSP.GAPS, HSP.IDENTITY.PERC, contains("HSP")) |>
  arrange(FAMILY, GENUS, SPECIES)

# for regression analysis - add  count locations per species and asv
fish_asv_at_locs <- fish_biodiv_blast |> group_by(SPECIES) |> summarize(LOC.PER.SPC = n_distinct(RESERVE.GROUP.LOCATION)) |> arrange(LOC.PER.SPC) 
fish_asv_at_locs <- fish_biodiv_blast |> select(ASV, EDNA.OBS.PRES, FAMILY, GENUS, SPECIES, RESERVE.GROUP.LOCATION, HSP.GAPS, HSP.IDENTITY.PERC, contains("HSP")) |>  left_join(fish_asv_at_locs)

# for regression analysis - get a column with non-nz species as per above
fish_asv_at_locs <- fish_asv_at_locs |> mutate(NOT.NZ = as.factor(ifelse( grepl("*", SPECIES, fixed = TRUE), TRUE, FALSE)))

# for regression analysis - re-scale percentages as per MdL - for CIs more easily understandable
fish_asv_at_locs <- fish_asv_at_locs |> mutate(HSP.IDENTITY.PERC = 100 * HSP.IDENTITY.PERC)

# for regression analysis -  save object for external inspection if desirable
saveRDS(fish_asv_at_locs, "/Users/paul/Documents/OU_eDNA/201028_Robjects/230514_999_r_summarize_results__data_spc_distribution_vs_quality.Rds")

# for manuscript - count fraction of NOT.NZ species among unique eDNA assignments
fish_asv_at_locs |> select(ASV, SPECIES, NOT.NZ)

fish_asv_at_locs |> select(ASV, NOT.NZ) |> distinct() |> group_by(NOT.NZ) |>
  arrange(NOT.NZ) # was: 92 ASV: 39 False (42.4% False) / 53 True (57% True)
# now: 96 ASV: 20 False (21%  False) / 76 True (79% True) 

fish_asv_at_locs |> select(SPECIES, NOT.NZ) |> distinct() |> group_by(NOT.NZ) |>
  arrange(NOT.NZ) # was: 44 SPECIES: 25 False (56.8% True) / 19 True (43.1% True)
# now: 43 SPECIES:  9 False (32.5% True) / 34 True (67.4% True)

# ____ Regression analysis of alignment parameters ---- 

# ____ - Inspect modelling data

glimpse(fish_asv_at_locs)

# ____ - Build logistic regressions

# Null model

glm_mod_0 <-  glm(NOT.NZ ~ 1, weights = LOC.PER.SPC, family = binomial, data = fish_asv_at_locs)
summary(glm_mod_0) # AIC: 385.11

# Gaps only

glm_mod_1 <-  glm(NOT.NZ ~ HSP.GAPS, weights = LOC.PER.SPC, family = binomial, data = fish_asv_at_locs)
summary(glm_mod_1) # AIC: 386.34

# Identity percentage only 

glm_mod_2 <-  glm(NOT.NZ ~ HSP.IDENTITY.PERC, weights = LOC.PER.SPC, family = binomial, data = fish_asv_at_locs)
summary(glm_mod_2) # AIC:  374.28

# Both gaps and coverage 

glm_mod_3 <-  glm(NOT.NZ ~ HSP.IDENTITY.PERC + HSP.GAPS, weights = LOC.PER.SPC, family = binomial, data = fish_asv_at_locs)
summary(glm_mod_3) # AIC: 351.97

coefficients(glm_mod_3)
logit2prob(coefficients(glm_mod_3))
confint(glm_mod_3)
logit2prob(confint(glm_mod_3))

exp(coefficients(glm_mod_3))/(1+exp(coefficients(glm_mod_3)))
exp(confint(glm_mod_3))/(1+ exp(confint(glm_mod_3)))

# ____ - Testing model significance

anova(glm_mod_0, glm_mod_1, test="Chisq") # Gaps only **not significant**
anova(glm_mod_0, glm_mod_2, test="Chisq") # Match percentage **significant**
anova(glm_mod_2, glm_mod_3, test="Chisq") # both together **significant**
anova(glm_mod_1, glm_mod_3, test="Chisq") # both together **significant**
anova(glm_mod_0, glm_mod_3, test="Chisq") # both together **significant**

# ____ - Ranking models

AIC(glm_mod_0)
AIC(glm_mod_1)
AIC(glm_mod_2)
AIC(glm_mod_3) # best, AIC 351.9707

# ____ - Inspect best model

jtools::summ(glm_mod_3)

plot <- plot_summs(glm_mod_3)
plot + theme_bw() + 
  labs(title = "Coefficient estimates", 
       subtitle = paste("R model formula: ", as.character(paste(deparse(formula(glm_mod_3), width.cutoff = 500), collapse=""))),
       x="estimate, incl. CI [log odds]", y = "model coefficients")

ggsave("/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_new_analysis_outputs/230515_999_logistic_rgression_alignmnet_parameters.pdf", scale = 1.65, width = 4, height = 3, units = "in", dpi = 300)
ggsave("/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_si_di_development/7_model_coefficients_v1.pdf", scale = 1.65, width = 4, height = 3, units = "in", dpi = 300)


# ____ - Other model reporting

# check_observation among true/false for figure legend
fish_asv_at_locs |> group_by(NOT.NZ) |> summarise(across(c("SPECIES", "ASV"), list(n_distinct)))

coeff_plot <- sjPlot::plot_model(glm_mod_3, vline.color = "red", show.values = TRUE, type = "est") +
  theme_bw() +
  scale_x_discrete(labels = c("Algn. Cov.", "Gaps")) +
  ggtitle("Relationship between non-native status and alignment quality")

ggsave("230515_999_r_summarize_results__coeff_plot.pdf", plot = coeff_plot, 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 0.75, width = 200, height = 135, units = c("mm"),
       dpi = 500, limitsize = TRUE)

ggsave("7_model_coefficients_v2.pdf", plot = coeff_plot, 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_si_di_development",
       scale = 0.75, width = 200, height = 135, units = c("mm"),
       dpi = 500, limitsize = TRUE)

model_plot <- plot_model(glm_mod_3, type = "pred", terms = c("HSP.IDENTITY.PERC", "HSP.GAPS"), show.data = TRUE, jitter = 0.0, ci.lvl = 0.95) +
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

ggsave("230515_999_r_summarize_results__asv_bin_regression.pdf", plot = last_plot(), 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 1, width = 200, height = 135, units = c("mm"),
       dpi = 500, limitsize = TRUE)

tab_model(glm_mod_3)

# # ____ - Save state

# save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__moidelling_done.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__moidelling_done.Rdata")

# __d) Inspected species accumulation curves ----

# Done by MDL

# __d) Calculated Good-Turing estimators ----

# Done by MDL

# __f) Analysed divergence of species observations by data source ----

# ____ Get Figure 3 - Euler plots ----

# "Concordance of taxonomic information between the four data sources diverged with lowering taxonomic levels, most pronounced among eDNA data (Fig. 3"

# ____ - Create Euler plots ----
euler_obs_fish_bio <- lapply(list("SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"), get_euler_object, fish_biodiv)
euler_ggp_fish_bio <- mapply(get_euler_ggplot, list("SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"), euler_obs_fish_bio, plot_label = FALSE, SIMPLIFY = FALSE)

# ____ - Create compound plot with better labels then with plot_label = TRUE above ----
ggarrange( plotlist = euler_ggp_fish_bio[4:7],
           labels = str_to_sentence(c("ORDER", "FAMILY", "GENUS", "SPECIES")),
           font.label = list(size = 12, color = "black", face = "bold.italic", family = NULL),
           ncol = 2, nrow = 2)

# ____ - Save compound plot with better labels then with plot_label = TRUE above ----
ggsave("230515_999_r_summarize_results__euler_edna_bruv_obis.pdf", plot = last_plot(), 
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
       scale = 2.0, width = 100, height = 100, units = c("mm"),
       dpi = 500, limitsize = TRUE)


# Unused code bits possibly needed for some later revisions ----

# _1.) List non-unique and unique species observations in Literature and OBIS ----

# Species in OBIS and literature 
pbob_species <- fish_biodiv |> select(SPECIES, SAMPLE.TYPE) |> filter(SAMPLE.TYPE %in% c("OBIS", "PUBL")) |> distinct()
unique(pbob_species[["SPECIES"]]) # Species in OBIS and literature - 70 

# _2.) eDNA species that are also in OBIS or Literature ----

sum(pbob_species[["SPECIES"]] %in% edna_species[["SPECIES"]]) # 2

# these are those two species - in literure - not in OBIS
intersect(get_clean_strings(pbob_species[["SPECIES"]]), get_clean_strings(edna_species[["SPECIES"]])) # "Aldrichetta forsteri" "Thyrsites atun"
intersect(get_clean_strings(obis_species[["SPECIES"]]), get_clean_strings(edna_species[["SPECIES"]])) # "Aldrichetta forsteri" "Thyrsites atun"
intersect(get_clean_strings(publ_species[["SPECIES"]]), get_clean_strings(edna_species[["SPECIES"]])) # "Aldrichetta forsteri" "Thyrsites atun"

# _3.) OBIS or Literature species not detected by eDNA ----

setdiff(get_clean_strings(pbob_species[["SPECIES"]]), get_clean_strings(edna_species[["SPECIES"]]))

# _4.) eDNA species that are not in OBIS or Literature ----

setdiff(get_clean_strings(edna_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]]))

# _5.) eDNA species that are not in OBIS or Literature, not listed in Roberts 2020----

intersect(setdiff(get_clean_strings(edna_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish)

# _6.) eDNA species that are not in OBIS or Literature, listed in Roberts 2020 ----

setdiff(setdiff(get_clean_strings(edna_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish)

# _7.) eDNA species that are not in OBIS or Literature, listed in Roberts 2020, detected by MEGAN ----

intersect(setdiff(setdiff(get_clean_strings(edna_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish), megan_fish_species)

# _8.) eDNA species that are not in OBIS or Literature, listed in Roberts 2020, not detected by MEGAN ----

setdiff(setdiff(setdiff(get_clean_strings(edna_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish), megan_fish_species)

# _9.) OBIS or Literature species not detected by BRUV ----

setdiff(get_clean_strings(pbob_species[["SPECIES"]]), get_clean_strings(bruv_species[["SPECIES"]]))

# _10.) BRUV species that are not in OBIS or Literature ----

setdiff(get_clean_strings(bruv_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]]))

# _11.) BRUV species that are not in OBIS or Literature, not listed in Roberts 2020----

intersect(setdiff(get_clean_strings(bruv_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish)

# _12.) BRUV species that are not in OBIS or Literature, listed in Roberts 2020 ----

setdiff(setdiff(get_clean_strings(bruv_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish)

# _13.) BRUV species that are not in OBIS or Literature, listed in Roberts 2020, detected by MEGAN ----

intersect(setdiff(setdiff(get_clean_strings(bruv_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish), megan_fish_species)

# _14.) BRUV species that are not in OBIS or Literature, listed in Roberts 2020, not detected by MEGAN ----

setdiff(setdiff(setdiff(get_clean_strings(bruv_species[["SPECIES"]]), get_clean_strings(pbob_species[["SPECIES"]])), nonnz_fish), megan_fish_species)

# _16.) Save state  ----

# save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__summaries_done.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__summaries_done.Rdata")

# Results for main text (2 of 2)   ----

# _1.) Get and Export Expanded Flextable with BLAST results (for Table 1 and Figure 2)----

# __a) Creating table columns for manuscript: HSP.IDENTITY.PERC ----

fish_biodiv_blast_cov <- fish_biodiv_blast |> 
  group_by(SPECIES) |> 
  summarize(BLAST.COV.RNG =  ifelse( signif(100*min(HSP.IDENTITY.PERC), 3) != signif(100*max(HSP.IDENTITY.PERC),3),
                                     paste0( signif(100*min(HSP.IDENTITY.PERC), 3), "-", signif(100*max(HSP.IDENTITY.PERC),3), "%"),
                                     paste0( signif(100*min(HSP.IDENTITY.PERC),3), "%")
                                     ))

# __b) Creating table columns for manuscript: HSP.GAPS ----

fish_biodiv_blast_gap <- fish_biodiv_blast |> 
  group_by(SPECIES) |> 
  summarize(BLAST.GAP.RNG =  ifelse( signif( min(HSP.GAPS), 3) != signif(max(HSP.GAPS),3),
                                     paste0( signif(min(HSP.GAPS), 3), "-", signif(max(HSP.GAPS),3)),
                                     paste0( signif(min(HSP.GAPS), 3))
                                     ))

# __c) Creating table columns for manuscript: AVG.HSP.BIT.SCORE  ----

# 15-05-2023 - adding  Bit score parmeters
glimpse(fish_biodiv_blast)

fish_biodiv_blast_avgbitsc <- fish_biodiv_blast |> 
  group_by(SPECIES) |> 
  summarize(BLAST.AVG.RNG = case_when(
    signif(min(AVG.HSP.BIT.SCORE), 2) == signif(max(AVG.HSP.BIT.SCORE), 2) ~ as.character(signif(max(AVG.HSP.BIT.SCORE),2)), 
    signif(min(AVG.HSP.BIT.SCORE), 2) != signif(max(AVG.HSP.BIT.SCORE), 2) ~ paste0(signif(min(AVG.HSP.BIT.SCORE),2), "-" , signif(max(AVG.HSP.BIT.SCORE),2))
    )
    )
              
# __d) Creating table columns for manuscript: MAX.HSP.BIT.SCORE  ----

fish_biodiv_blast_maxbitsc <- fish_biodiv_blast |> 
  group_by(SPECIES) |> 
  summarize(BLAST.MAX.RNG = case_when(
    signif(min(MAX.HSP.BIT.SCORE), 2) == signif(max(MAX.HSP.BIT.SCORE), 2) ~ paste0("(", as.character(signif(max(MAX.HSP.BIT.SCORE),2)), ")"), 
    signif(min(MAX.HSP.BIT.SCORE), 2) != signif(max(MAX.HSP.BIT.SCORE), 2) ~ paste0("(", signif(min(MAX.HSP.BIT.SCORE),2), "-" , signif(max(MAX.HSP.BIT.SCORE),2),")" )
  ))

# __e) Get useful bit score columns ----

# use this to use all values
fish_biodiv_blast_bitsc <- left_join(fish_biodiv_blast_avgbitsc, fish_biodiv_blast_maxbitsc)
fish_biodiv_blast_bitsc <-  fish_biodiv_blast_bitsc |> mutate(BLAST.BSC.RNG = paste0(as.character(BLAST.AVG.RNG), " ", as.character(BLAST.MAX.RNG))) |> select(-c("BLAST.MAX.RNG", "BLAST.BSC.RNG"))

# use this to use average values only
fish_biodiv_blast_bitsc <- left_join(fish_biodiv_blast_avgbitsc, fish_biodiv_blast_maxbitsc)
fish_biodiv_blast_bitsc <- fish_biodiv_blast_bitsc |> mutate(BLAST.BSC.RNG = BLAST.AVG.RNG) |> select(-c("BLAST.MAX.RNG", "BLAST.AVG.RNG"))

# __f) Extend table plot ----

htmp_tibl_fish <- bind_rows(
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "eDNA", tbl = TRUE) %>% add_column(SAMPLE.TYPE = "eDNA"),
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "BRUV", tbl = TRUE) %>% add_column(SAMPLE.TYPE = "BRUV"), 
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "OBIS", tbl = TRUE) %>% add_column(SAMPLE.TYPE = "OBIS"),
  get_matrix_or_table(fish_biodiv_tbls, obs_methods = "PUBL", tbl = TRUE) %>% add_column(SAMPLE.TYPE = "PUBL"),
)


htmp_tibl_fish_blrngs <- htmp_tibl_fish |> left_join(fish_biodiv_blast_cov) |> left_join(fish_biodiv_blast_gap) |> left_join(fish_biodiv_blast_bitsc)

# __g) Get Flextable ----

# format data for flex table
tibl_plot <- htmp_tibl_fish_blrngs %>% select(PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, TRIVIAL.SPECIES, BLAST.COV.RNG, BLAST.GAP.RNG, BLAST.BSC.RNG) %>% 
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
     BLAST.GAP.RNG = "Algn. gaps",
     BLAST.BSC.RNG = "Algn. bitsc."
     )) %>%
     fontsize(part = "all", size = 12) %>%
     autofit(add_w = 1, add_h = 0, part = c("all"))

save_as_html(ft, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__spcies_obs_matching_tiles.html")
save_as_docx(ft, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__spcies_obs_matching_tiles.docx")

# _2.) Get Figure 2 - Tile display item of all biodiversity data ----

#  __a) get the heat map (Fig 2) ----

# order factors in plotting object to match flex table  
y_axis_label_order <- fct_relevel(htmp_tibl_fish$SPECIES, rev(c(tibl_plot$SPECIES)))  # y_axis_label_order <- reorder(htmp_tibl_full$SPECIES, desc(htmp_tibl_full$SPECIES))

# plot out data as tiles
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
          axis.text.y = element_text(angle = 0, hjust = 1, size = 7, face = "italic"), 
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

ggsave("230515_999_r_summarize_results__biodiv_tiles_only.pdf", plot = plot_htmp_adjusted, 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 155, height = 297, units = c("mm"),
         dpi = 500, limitsize = TRUE)


#  __b) Combine heat map and flextable ----

# get flex table as ggplot object
ft_raster <- as_raster(ft) # webshot and magick.
ft_rgrob  <- rasterGrob(ft_raster)
plot_ft <- ggplot() + theme_void() + annotation_custom(ft_rgrob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

ggarrange(plot_ft, plot_htmp_adjusted, ncol = 2, nrow = 1, labels = c("a","b"), widths = c(1,1))

ggsave("230515_999_r_summarize_results__biodiv_heat.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1.5, width = 210, height = 297, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# __c) Save state ----

# save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__dis_done.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__dis_done.Rdata")

# _3.) Results – Show that local reference data is needed to address ecological questions in the Fiordland region ----

# __a.) ANOSIM of observation types and variables ----

# ____ Set up analysis ----

# subset to local observation without OBIS, OBIS data is too sparse to allow meaningful conclusions
fish_biodiv_local <- fish_biodiv |> filter(SET.ID %!in% c(98, 99)) |> filter(SAMPLE.TYPE %in% c("eDNA", "BRUV", "OBIS"))

# testing function with one data set
get_vegan(distance = "jaccard", tibl = fish_biodiv_local, group_col = "SET.ID", group_row = c("SPECIES"), group_col_ano = "RESERVE.GROUP.LOCATION", obs_methods = "BRUV")


# Analysis for fish (as per eDNA markers) for BRUV and eDNA (data complete across all sets)

# setting up parameter combinations for complete ANOSIM analysis
anosim_analysis_fish <- expand.grid(
  distance      = c("jaccard"),
  tibl          = c("fish_biodiv_local"),
  group_col     = c("SET.ID"), 
  group_row     = c("SPECIES", "GENUS", "FAMILY", "ORDER"),
  group_col_ano = c("RESERVE.GROUP.LOCATION", "RESERVE.GROUP.INSIDE"), 
  obs_methods   = c("eDNA", "BRUV")
  )

# ____ Run ANOSIM analysis ----

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
     highlight(j = ~ significance, color = function(x) {ifelse(x < 0.05, "lightgreen", "white")} ) %>%
     fontsize(part = "all", size = 9) %>%
     fit_to_width(max_width = 12, inc = 1L, max_iter = 20)

# ____  Export ANOSIM results ----  

# save_as_html(ft_anosim, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210712_998_r_summarize_results__ANOSIM.html")
  save_as_html(ft_anosim, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__ANOSIM.html")
  save_as_docx(ft_anosim, path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__ANOSIM.docx")

# ____ Save state ----  

# save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__anaosim_done.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__anaosim_done.Rdata")
  
# __b.) Indicator species analysis for significant ANOMSIM results  ----

# ____  Testing indicator species analysis ----
  
# - note flag "map" set to TRUE
get_vegan(distance = "jaccard", tibl = fish_biodiv_local, group_col = "SET.ID", group_row = c("SPECIES"), group_col_ano = "RESERVE.GROUP.LOCATION", obs_methods = "BRUV", mp = TRUE)

# ____  Running analysis ----

# Analysis for fish (as per eDNA markers) for BRUV and eDNA (data complete across all sets)

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

# ____  Inspect analysis results ----

# inspect results
str(mpatt_results_fish)
str(mpatt_results_fish[[1]])
mpatt_results_fish[[1]] # species level - order follows rows in mpatt_analysis_fish
mpatt_results_fish[[2]]
mpatt_results_fish[[3]]
mpatt_results_fish[[4]]

summary(mpatt_results_fish[[1]])

# Multilevel pattern analysis
#   
# Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 26
# Selected number of species: 1 
# Number of species associated to 1 group: 0 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# Number of species associated to 5 groups: 0 
# 
# List of species associated to each combination: 
#   
# Group WJ CTRL+WJ MR  #sps.  1 
# stat p.value  
# Bodianus unimaculatus 0.725  0.0304 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(mpatt_results_fish[[2]])

# Multilevel pattern analysis
#   
# Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 24
# Selected number of species: 1 
# Number of species associated to 1 group: 0 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# Number of species associated to 5 groups: 0 
# 
# List of species associated to each combination: 
#   
# Group WJ CTRL+WJ MR  #sps.  1 
# stat p.value  
# Bodianus 0.725  0.0285 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

summary(mpatt_results_fish[[3]])

# Multilevel pattern analysis
#   
# Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 19
# Selected number of species: 1 
# Number of species associated to 1 group: 0 
# Number of species associated to 2 groups: 0 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 1 
# Number of species associated to 5 groups: 0 
# 
# List of species associated to each combination: 
#   
# Group FF CTRL+FF MR+WJ CTRL+WJ MR  #sps.  1 
# stat p.value  
# Labridae 0.67  0.0418 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

summary(mpatt_results_fish[[4]])

# Multilevel pattern analysis
#   
# Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 12
# Selected number of species: 1 
# Number of species associated to 1 group: 0 
# Number of species associated to 2 groups: 0 
# Number of species associated to 3 groups: 1 
# Number of species associated to 4 groups: 0 
# Number of species associated to 5 groups: 0 
# 
# List of species associated to each combination: 
#   
# Group FF MR+WJ CTRL+WJ MR  #sps.  1 
# stat p.value  
# Perciformes 0.699  0.0299 *
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# _4.) Save state ----

# save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__multipatt_done.Rdata")
# load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__multipatt_done.Rdata")

