#  *******************************************************
#  *                                                     *   
#  *  Including Fiordland species lists from literature  *
#  *                                                     *
#  *******************************************************
#
#  Paul Czechowski - paul.czechowski@gmail.com

# I. Environment setup 
# ====================
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
rm(list = ls(all.names = TRUE))
gc()

library("taxize")
library("dplyr")
library("tidyr")

Sys.setenv(ENTREZ_KEY="a634c6e9c96c3859bca27a2771f6d2872f08")
Sys.getenv("ENTREZ_KEY")

options(tibble.print_max = Inf) 
`%!in%` <- Negate(`%in%`)

# II. Data read-in  
# ================

# Define species from literature
# ------------------------------
#
# "fish" (Actinopterygii, Elasmobranchii) species records from 
# 
#  Milford Sound First baseline survey for non-indigenous marine 
#  species (Research Project ZBS2005/19) 
#  MAF Biosecurity New Zealand Technical Paper No: 2008/01 
#  https://niwa.co.nz/static/marine-biosecurity/Inglis%20et%20al%202008%20milford%20resurvey%20report.pdf
#
#  and **many references** therein!
# 
#  Wing, S. R., and L. Jack. 2013. Marine reserve networks conserve biodiversity 
#  by stabilizing communities and maintaining food web structure. Ecosphere 4(11):135. 
#  http://dx.doi.org/10.1890/ES13-00257.1 

spc_in <- c(
  "Acanthoclinus fuscus", "Acanthoclinus littoreus", "Acanthoclinus marilynae", "Acanthoclinus matti", 
  "Acanthoclinus rua", "Aldrichetta forsteri", "Aplidium adamsi", "Aplodactylus arctidens", "Aplodactylus arctidens",
  "Atherinomorus lacunosus", "Bellapiscis lesleyae", "Bellapiscis medius", "Bovichtus variegatus",
  "Caesioperca lepidoptera", "Callanthias allporti", "Cephaloscyllium isabellum", "Cominella sp.", 
  "Conger verreauxi", "Cryptichthys jojettae", "Eptatretus cirrahtus", "Fiordichthys slartibartfasti", 
  "Forsterygion flavonigrum", "Forsterygion lapillum", "Forsterygion malcolmi", "Forsterygion varium", 
  "Gaidropsarus novaezelandi", "Modicus minimus", "Gobiopsis atrata", "Grahamina capito", "Helicolenus percoides", 
  "Hypoplectrodes huntii", "Karalepis stewarti", "Latridopsis ciliaris", "Latridopsis forsteri", "Latris lineata", 
  "Lepidoperca tasmanica", "Lissocampus filum", "Lotella rhacina", "Mendosoma lineatum", "Mendosoma lineatum", 
  "Modicus tangaroa", "Nemadactylus macropterus", "Notoclinops caerulepunctus", "Notoclinops segmentatus", 
  "Notoclinus fenestratus", "Notolabrus celidotus", "Notolabrus cinctus", "Notolabrus fucicola", "Forsterygion maryannae", 
  "Odax pullus", "Parapercis colias", "Parapercis gilliesii", "Paratrachichthys trailli", "Meuschenia scaber", 
  "Patiriella regularis", "Polyprion oxygeneios", "Pseudolabrus miles", "Pseudophycis barbata",
  "Retropinna retropinna", "Rhombosolea plebeia", "Ruanoho decemdigitatus",
  "Ruanoho whero", "Scorpaena papillosa", "Scorpis lineolata", "Squalus acanthias", "Thyrsites atun"
  ) 
  
# Download NCBI annotations (to match taxonomy records of other surveys)
# ----------------------------------------------------------------------

# sort strings
spc <- spc_in  |> sort() |> unique()   # for assembly of object 
gen <- gsub( " .*$", "", spc)          # for assembly of object 
gen_uniq <- gen |> sort() |> unique()  # for reporting only (below)

# download annotations - list of data frames
gen_ls <- classification(gen, db = "ncbi")
spc_ls <- classification(spc, db = "ncbi")
gen_uniq_list <- classification(gen_uniq, db = "ncbi")

# 58 of 64 species found in NCBI (and 22 not found in NCBI excluded from eDNA)
# ---------------------------------------------------------------------------
spc_in_nf <- tibble( SPECIES = c(
  "Acanthoclinus littoreus",
  "Acanthoclinus marilynae",
  "Acanthoclinus matti",
  "Acanthoclinus rua",
  "Aplidium adamsi",
  "Callanthias allporti",
  "Cephaloscyllium isabellum",
  "Cominella sp.",
  "Cryptichthys jojettae",
  "Gaidropsarus novaezelandi",
  "Hypoplectrodes huntii",
  "Eptatretus cirrahtus",
  "Fiordichthys slartibartfasti",
  "Forsterygion malcolmi",
  "Forsterygion maryannae",
  "Gaidropsarus novaezelandi",
  "Modicus minimus",
  "Modicus tangaroa",
  "Notoclinops caerulepunctus",
  "Notoclinops segmentatus",
  "Notoclinus fenestratus",
  "Ruanoho decemdigitatus" 
   ))

# 44 of 48 genera found in NCBI (and 4 not found in NCBI excluded from eDNA)
# ---------------------------------------------------------------------------

# already among species names - not needed
# gen_in_nf <- tibble( GENUS = c(
#   "Cryptichthys",
#   "Fiordichthys",
#   "Modicus",
#   "Notoclinops"
#   ))

# III. Data formatting
# ====================

# get basic species lists
# -----------------------

# get list of tibbles and strip class attributes for row binding
spc_tibl_ls <- sapply(spc_ls, as_tibble)
gen_tibl_ls <- sapply(gen_uniq_list, as_tibble)

length(gen_tibl_ls) # 48 genera
length(spc_tibl_ls) # 64 spcies

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
unique(spc$GENUS); unique(spc$SPECIES) # 34 Genera, 42 species fully resolved
gen_not_in_spc <- gen |> filter(GENUS %!in% unique(spc$GENUS))
# add species not found and later fill higher taxonomy columns up 
spc <- bind_rows(spc, spc_in_nf, gen ) 

#  add missing species manually for which genus information is available
# ----------------------------------------------------------------------

# sort stuff
col_order <- c("SUPERKINGDOM", "PHYLUM",  "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")                     
spc <- spc |> mutate(GENUS = case_when( is.na(GENUS) == TRUE ~  gsub( " .*$", "", SPECIES), TRUE ~ as.character(GENUS)))
spc <- spc |> relocate(col_order) |> arrange(across( rev(col_order[1:(length(col_order)-1)]) ))

# fill NA's that can be filled
spc <- spc |> group_by(GENUS)  |> fill(FAMILY, .direction = c("updown"))
spc <- spc |> group_by(FAMILY) |> fill(ORDER, .direction = c("updown"))
spc <- spc |> group_by(ORDER)  |> fill(CLASS, .direction = c("updown")) 
spc <- spc |> group_by(CLASS) |> fill(PHYLUM, .direction = c("updown"))
spc <- spc |> group_by(PHYLUM) |> fill(SUPERKINGDOM, .direction = c("updown"))


# plug holes manually 
spc <- spc |> mutate(ORDER = ifelse(FAMILY == "Plesiopidae", "Ovalentaria", ORDER)) 

spc <- spc |> mutate(FAMILY = ifelse(GENUS == "Cryptichthys", "Tripterygiidae", FAMILY)) 
spc <- spc |> mutate(ORDER  = ifelse(GENUS == "Cryptichthys", "Blenniiformes" , ORDER)) 

spc <- spc |> mutate(FAMILY = ifelse(GENUS == "Fiordichthys", "Bythitidae", FAMILY)) 
spc <- spc |> mutate(ORDER  = ifelse(GENUS == "Fiordichthys", "Ophidiiformes" , ORDER)) 

spc <- spc |> mutate(FAMILY = ifelse(GENUS == "Modicus", "Gobiesocidae", FAMILY)) 
spc <- spc |> mutate(ORDER  = ifelse(GENUS == "Modicus", "Gobiesociformes" , ORDER)) 

spc <- spc |> mutate(FAMILY = ifelse(GENUS == "Notoclinops", "Tripterygiidae", FAMILY)) 
spc <- spc |> mutate(ORDER  = ifelse(GENUS == "Notoclinops", "Blenniiformes" , ORDER)) 


spc <- spc |> mutate(FAMILY = ifelse(GENUS == "Callanthias", "Callanthiidae", FAMILY)) 
spc <- spc |> mutate(ORDER = ifelse(GENUS == "Callanthias", "Eupercaria", ORDER)) 


# fill NA's that can be filled
spc <- spc |> relocate(col_order) |> arrange(across( rev(col_order[1:(length(col_order)-1)]) ))
spc <- spc |> group_by(GENUS)  |> fill(FAMILY, .direction = c("updown"))
spc <- spc |> group_by(FAMILY) |> fill(ORDER, .direction = c("updown"))
spc <- spc |> group_by(ORDER)  |> fill(CLASS, .direction = c("updown")) 
spc <- spc |> group_by(CLASS)  |> fill(PHYLUM, .direction = c("updown"))
spc <- spc |> group_by(PHYLUM) |> fill(SUPERKINGDOM, .direction = c("updown"))

# remove holes that can be
spc <- spc |> filter(!is.na(SPECIES)) |> distinct()

# sort again
spc <- spc |> relocate(col_order) |> arrange(across(col_order))


# add other variables for downstream compatibility
# -----------------------------------------------

spc <- spc |> mutate(NCBI.LEVEL = ifelse(GENUS != "Cominella", "species", "genus"))
spc <- spc |> mutate(NCBI.TAXID = ifelse(!is.na(NCBI.TAXID), NCBI.TAXID, as.character("0"))) 
spc <- spc |> mutate(NCBI.TAXID = as.numeric(NCBI.TAXID))
spc <- spc |> mutate(NCBI.TAXID.INC = ifelse(NCBI.TAXID == 0, TRUE, FALSE)) 
spc <- spc |> mutate(SAMPLE.TYPE = "PUBL") |> mutate(ABUNDANCE = 1) |> mutate(SET.ID = 98)
spc <- spc |> mutate(PUBL.OBS.PRES = 1)
spc <- spc |> mutate(LOC.NAME = "Fiordland")
spc <- spc |> mutate(RESERVE.GROUP = "FI")
spc <- spc |> mutate(RESERVE.GROUP.INSIDE = FALSE)
spc <- spc |> mutate(RESERVE.GROUP.LOCATION = "FI CTRL")


# IV. Data export
# ===============
publ_obs <- ungroup(spc)

# save or load molten state 
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210723_995_r_get_PUBL_long_table__image.Rdata")
saveRDS(publ_obs, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210723_995_r_get_PUBL_long_table__publ_obs.Rds")
saveRDS(publ_obs, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210723_995_r_get_PUBL_long_table__publ_obs.Rds")

# Appendix. Unused species
# =========================

# MSNIS - NIWA report - not fish
#  
# Champia  affinis
# Crella incrustans
# Didemnum sp.
# Alexandrium affine
# Alexandrium minutum
# Alexandrium ostenfeldii
# Alexandrium tamarense
# Alexandrium catenella
# Gymnodinium catenatum
# Diplosoma velatum
# Esperiopsis edwardii 
# Haliclona clathrata
# Heterosigma akashiwo
# Leucosolenia challengeri
# Leucosolenia discoveryi
# Orthopyxis integra
# Polysiphonia brodiei
# Polysiphonia constricta
# Polysiphonia sertularioides
# Polysiphonia subtilissima
# Raspailia agminata 
# Sargassum verruculosum
# Scruparia ambigua
# Scruparia ambigua
# Tethya bergquistae
# 
# Aglaophamus macroura
# Owenia petersenae
# Scoloplos simplex
# Prionospio australiensis
# Timarete anchylochaetus
# Waitangi rakiura
# Myadora striata
# Scalpomactra scalpellum
# Notocallista multistriata
# Austrofusus glans
# Amalda australis
# Amalda novaezelandiae
# Pervicacia tristis
# Zeacumantus subcarinatus
# Antisolarium eoenum

