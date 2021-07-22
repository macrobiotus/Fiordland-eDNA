# 
#
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
  "Hypoplectrodes huntii", "Karalepis stewarti", "Latridopsis ciliaris", "Latridopsis forsteri ", "Latris lineata", 
  "Lepidoperca tasmanica", "Lissocampus filum", "Lotella rhacina", "Mendosoma lineatum", "Mendosoma lineatum", 
  "Modicus tangaroa", "Nemadactylus macropterus", "Notoclinops caerulepunctus", "Notoclinops segmentatus", 
  "Notoclinus fenestratus", "Notolabrus celidotus", "Notolabrus cinctus", "Notolabrus fucicola", "Forsterygion maryannae", 
  "Odax pullus", "Parapercis colias", "Parapercis gilliesii", "Paratrachichthys trailli", "Meuschenia scaber", 
  "Patiriella regularis", "Polyprion oxygeneios", "Pseudolabrus miles ", "Pseudophycis barbata",
  "Retropinna retropinna", "Rhombosolea plebeia", "Ruanoho decemdigitatus",
  "Ruanoho whero", "Scorpaena papillosa", "Scorpis lineolata", "Squalus acanthias", "Thyrsites atun"
  ) 
  
# Download NCBI annotations (to match taxonomy records of other surveys)
# ----------------------------------------------------------------------

# sort strings
spc <- species_in  |> sort() |> unique() # for assembly of object 
gen <- gsub( " .*$", "", species)             # for assembly of object 
gen_uniq <- gen |> sort() |> unique()       # for reporting only (below)

# download annotations - list of data frames
gen_ls   <- classification(genus, db = "ncbi")
spc_ls <- classification(species, db = "ncbi")
gen_uniq_list <- classification(genus_uniq, db = "ncbi")

# 58 of 64 species found in NCBI (and 22 not found in NCBI excluded from eDNA)
# ---------------------------------------------------------------------------

# Acanthoclinus littoreus
# Acanthoclinus marilynae
# Acanthoclinus matti
# Acanthoclinus rua
# Aplidium adamsi
# Callanthias allporti
# Cephaloscyllium isabellum
# Cominella sp.
# Cryptichthys jojettae
# Gaidropsarus novaezelandi
# Hypoplectrodes huntii
# Eptatretus cirrahtus
# Fiordichthys slartibartfasti
# Forsterygion malcolmi
# Forsterygion maryannae
# Gaidropsarus novaezelandi
# Modicus minimus
# Modicus tangaroa
# Notoclinops caerulepunctus
# Notoclinops segmentatus
# Notoclinus fenestratus
# Ruanoho decemdigitatus

# 44 of 48 genera found in NCBI (and 4 not found in NCBI excluded from eDNA)
# ---------------------------------------------------------------------------
# Cryptichthys
# Fiordichthys
# Modicus
# Notoclinops

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

# select genera that haven't been found on species level already
unique(spc$GENUS); unique(spc$SPECIES) # 34 Genera, 42 species fully resolved
gen_not_in_spc <- gen |> filter(GENUS %!in% unique(spc$GENUS))
spc <- bind_rows(spc, gen_not_in_spc ) 

#  add missing species manually for which genus information is available
# ----------------------------------------------------------------------

# continue here after 22-Jul-2021


# add other variables for downstream compatibility
# -----------------------------------------------

# SET.ID will be PUBL



# get counts and verify table
# ---------------------------


# IV. Data export
# ===============


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

