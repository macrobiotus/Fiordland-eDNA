# **********************************************
# * Create, filter, and write Physloseq object *
# **********************************************
# 12-Dec-2020, 08-Jan-2020, 11-Jan-2020, 20-Jan-2020

# load packages
# =============
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse") # work using tibbles
library("phyloseq") # handle Qiime 2 data in R - best used to generate long dataframe

library("data.table") # faster handling of large tables - such as molten Phyloseq objects
library("future.apply") # faster handling of large tables - such as molten Phyloseq objects

library("decontam") # decontamination - check `https://benjjneb.github.io/decontam/vignettes/decontam_intro.html`
library("openxlsx") # write Excel tables

# functions
# =========

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)

# Create phyloseq object.
# -----------------------

# Creates `phyloseq` objects from Qiime` compatible data.
get_phsq_ob <- function(biom_path, sequ_path, tree_path) {

	# read data into R 
	btab <- phyloseq::import_biom(biom_path)
	# tree <- ape::read.tree(tree_path)
	sequ <- Biostrings::readDNAStringSet(sequ_path)

	# construct object  
	phsq_ob <- phyloseq::merge_phyloseq(btab, tree_path, sequ)

	# return object
	return(phsq_ob)
}

# Remove empty data
# -----------------

# This function removes "0" count phylotypes from samples and samples with "0"
# phylotypes.
remove_empty <- function(phsq_ob){

  # filter Phylotypes
  phsq_ob <- phyloseq::prune_taxa (taxa_sums (phsq_ob) > 0, phsq_ob)
  
  # filter samples
  phsq_ob <- phyloseq::prune_samples (sample_sums (phsq_ob) > 0, phsq_ob)
  
  # return object
  return (phsq_ob)
}


# Clean taxonomy strings
# ----------------------

clean_molten_tax_strings <- function(molten_phyloseq_object){

  require(dplyr)
  
  molten_phyloseq_object <- molten_phyloseq_object %>% 
    mutate(across(c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), .funs = gsub, pattern = "__", replacement = " "))

  molten_phyloseq_object <- molten_phyloseq_object %>% 
    mutate(across(c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), .funs = gsub, pattern = "_", replacement = " "))

  molten_phyloseq_object <- molten_phyloseq_object %>% 
    mutate(across(c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), .funs = gsub, pattern = "  ", replacement = " "))

  molten_phyloseq_object <- molten_phyloseq_object %>% arrange(superkingdom, phylum, class, order, family, genus, species)
  
  return(molten_phyloseq_object)

}

# Melt and tidy Phyloseq object
# -----------------------------

get_tidy_molten_ps = function(psob){

  require("tidyverse")
  library("phyloseq")
  
  psob_molten <- psmelt(psob)
  
  # get a tibble
  psob_molten <- psob_molten %>% as_tibble(.)

  # remove empty data 
  psob_molten <- psob_molten %>% filter(Abundance > 0) 

  # clean taxonomy strings
  psob_molten <- clean_molten_tax_strings(psob_molten)

  # re-label available data categories for plotting
  psob_molten$phylum[which(psob_molten$phylum %in% c("nomatch"))]  <- "no Blast hit"
  psob_molten$phylum[which(psob_molten$phylum %in% c("undefined"))] <-  "missing taxonomy"
  
  # create uniform column names
  psob_molten <- psob_molten %>% rename(ASV = OTU) %>% rename_with(toupper)

  return(psob_molten)
  
}


# Get an easily digestible plot of a molten Phyloseq object (e.g. before, during, after filtering)
# ----------------------------------------------------------------------------------------------------

get_default_plot = function (psob_molten,taxlev, ptitl, pxlab, pylab ){

  require(ggplot2)
  
  ggplot(psob_molten, aes_string(x = taxlev, y = "ABUNDANCE", fill = taxlev)) +
    geom_bar(stat = "identity", position = "stack", colour = NA, size=0) +
    facet_grid(LOC.NAME ~ ., shrink = TRUE, scales = "fixed") +
    theme_bw() +
    theme(legend.position = "none") +
    theme(strip.text.y = element_text(angle=0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7), 
          axis.ticks.y = element_blank()) +
    labs( title = ptitl) + xlab(pxlab) + ylab(pylab)

}


# Get an easily digestible summary of a molten Phyloseq object (e.g. before, during, after filtering)
# ----------------------------------------------------------------------------------------------------

get_molten_ps_description <- function(ps, rank_level, rank_name){
   
  require("tidyverse")

  # total sequencing effort
  f_count <- sum(ps$ABUNDANCE)

  # total sample count 
  s_count <- length(unique(ps$SAMPLE))

  # total ASV count
  a_count <- length(unique(ps$ASV))

  # ASV count - matching search criterium 
  ss_count <- ps %>% filter( get(rank_level) %in% c(rank_name)) %>% distinct(ASV) %>% nrow()
  ss_string <- ps %>% filter( get(rank_level) %in% c(rank_name)) %>% distinct(get(rank_level)) %>% pull() %>% paste0(collapse = " ")

  # ASV count - not matching search criterium 
  nss_count <- ps %>% filter( get(rank_level) %!in% c(rank_name)) %>% distinct(ASV) %>% nrow()
  nss_string  <- ps %>% filter( get(rank_level) %!in% c(rank_name)) %>% distinct(get(rank_level)) %>% pull() %>% paste0(collapse = ",  ")

  # get sample coverages
  coverage_per_sample <- aggregate(ps$ABUNDANCE, by=list(Sample=ps$SAMPLE), FUN=sum)
  summary(coverage_per_sample)
  cov_s_min  <- summary(coverage_per_sample)[c(1) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_s_med  <- summary(coverage_per_sample)[c(3) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_s_mean <- summary(coverage_per_sample)[c(4) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_s_max <- summary(coverage_per_sample)[c(6) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)

  # get ASV coverages
  coverage_per_asv <- aggregate(ps$ABUNDANCE, by=list(ASV=ps$ASV), FUN=sum)
  # print(head(coverage_per_asv))
  summary(coverage_per_asv)
  cov_a_min  <- summary(coverage_per_asv)[c(1) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_a_med  <- summary(coverage_per_asv)[c(3) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_a_mean <- summary(coverage_per_asv)[c(4) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_a_max <- summary(coverage_per_asv)[c(6) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)

  # print summary text
  message(paste0("\nThe current data set contains ", f_count, " sequences across ", s_count, " samples and ",   
                 a_count, " ASV's (", ss_count," ", ss_string,
                 ", as well as ", nss_count, " non-", ss_string, ", i.e. ", nss_string,
                 "). Sample mean (min., med., max.) coverage is ", cov_s_mean, 
                 " reads (", cov_s_min, ", ", cov_s_med, ", ", cov_s_max, 
                 "), and ASV mean (min, median, max) coverage is ", cov_a_mean, 
                 " reads (", cov_a_min, ", ", cov_a_med, ", ", cov_a_max, ").\n"
                 ))

  # summarize distinct values across the long data frame
  message("\nIn the following summary variable names shown are hard-coded in function \"get_molten_ps_stats\": ") 
  show_vars <- c("ASV", "ABUNDANCE", "SAMPLE", "LOC.NAME", "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")
  ps %>% select(any_of(show_vars)) %>% summarize_all(n_distinct, na.rm = TRUE) %>% print()
  
  # get most ASV's (not species) ordered by frequency
  message("\nGetting and returning highest-covered ASV's (not: species) ordered by frequency:")
  ps_asv_list <- left_join(coverage_per_asv , ps, by = c("ASV" = "ASV")) %>%
    distinct_at(vars("ASV", "x", "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"))
    
  ps_asv_list %>% arrange(desc(x)) %>% head(., n = 10) %>% print()


  # get most spcies's (not ASV's) ordered by frequency
  message("\nGetting and returning highest-covered species (not: ASV's) ordered by frequency:")
  
  ps_species_list <- aggregate(ps$ABUNDANCE, by=list(PHYLUM=ps$PHYLUM, CLASS=ps$CLASS, ORDER=ps$ORDER, FAMILY=ps$FAMILY, GENUS=ps$GENUS, SPECIES=ps$SPECIES), FUN=sum) 
  ps_species_list %>% arrange(desc(x)) %>% head(., n = 10) %>% print()

}


# I. Import Phyloseq object
# =========================

# Read in
# -------

# https://rdrr.io/github/jbisanz/qiime2R/f/vignettes/vignette.Rmd
#  error message stems from "#"  in taxonomy table
psob <- get_phsq_ob(biom_path = "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/980_12S_single_end_ee3-tab_q2_export/features-tax-meta.biom", 
                    sequ_path = "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/980_12S_single_end_ee3-tab_q2_export/dna-sequences.fasta", 
	                  tree_path = NULL
	                  )

# Rewrite rank names, as they seem to have been lost
# --------------------------------------------------
#   rank names from `/Users/paul/Documents/OU_pcm_eukaryotes/Github/150_r_get_q2_tax-tab.r`
colnames(tax_table(psob)) <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
head(tax_table(psob))

# Check and correct sample metadata
# ---------------------------------

# check REDAME.md 08-Jan-2021 and commit `64f61205e54334e9282616b0fa221c5362aa338c`
# Spelling and order should be should
#   be `sampleid,key,ng-ul,lat-dd,long-dd,pool-content,primer-direction,primer-label,adapter-fwd,flupad,index-fwd,template-primer-fwd,complete-primer-fwd,loc-name,sample-time,sample-type,inside-reserve,sample-name,vol-l,depth-m,notes,xtr-date,row,plate,col,primer-direction_rev,primer-label-rev,adapter-rev,index-rev,template-primer-rev,complete-primer-rev`
#   as in `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/850_prep_q2_predictor-tab__metadata.tsv`
#   as written by `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh`

names(sample_data(psob)) # names are spelled correctly but ordered alphabetically
head(sample_data(psob))  # names are changed from "-" to "." and ordered alphabetically in head command - PhyloSeq sloppy coding

# correct Sample metadata, as seen in plot later on
# - some "loc.names" are "NA": positive controls, negative controls, unassigned welle
# - "Five Fingers MR" are should be "Five Fingers"
sample_data(psob)[which(sample_data(psob)$`loc-name` == "Five Fingers MR"), ]$loc.name <- "Five Fingers"
sample_data(psob)[which(sample_data(psob)$`loc-name` == "NA"), ]$loc.name <- "Otago University"
sample_data(psob)[which(sample_data(psob)$`loc-name` == "Undaria Zone"), ]$loc.name <- "\"Undaria Zone\""

# Save or load initial state after import
# ---------------------------------------
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210108_990_r_get_eDNA_phyloseq__import_workspace.Rdata")
save(psob, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210108_990_r_get_eDNA_phyloseq__imported-psob.Rdata")
save(psob, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/zenodo/210108_990_r_get_eDNA_phyloseq__imported-psob.Rds")


# II. Inspect and summarize raw data
# ==================================

# melt and tidy Phyloseq object 
(psob_molten <- get_tidy_molten_ps(psob))

# create and save plot


get_default_plot (psob_molten, taxlev = "PHYLUM", ptitl = "Phyla across all locations before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")

ggsave("210108_990_r_get_eDNA_phyloseq__psob-unfiltered-phylum-overview.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# get summary 
get_molten_ps_description(psob_molten, rank_level = "SUPERKINGDOM", rank_name = "Eukaryota")

# remove temporary object
rm(psob_molten)


# III. Clean data using blank information - and summarise  
# =======================================================

# melt and tidy Phyloseq object 
(psob_molten <- get_tidy_molten_ps(psob))

# Describe controls
# -----------------

# split controls - what sample type are available? 
unique(psob_molten$SAMPLE.TYPE) # "eDNA": eDNA sample
                                # "blank": eDNA collection blank - bottel with mol grade water
                                # "pcntrol-blnd" positive controls (fish blend)
                                # "pcntrol-zebra": positive controls (Danio rerio)
                                # "ncntrl-pcr": amplification blank
                                # "ncntrol-xtr"
                                # "empty": unassigned wells  

# split controls -  how many of each sample type are available? 
psob_molten %>% group_by(KEY) %>% summarise(KEY = KEY, CNTROL = SAMPLE.TYPE ) %>% distinct() %>% pull("CNTROL") %>% table()

# split controls
all_cntrl <- psob_molten %>% filter(SAMPLE.TYPE %in% c( "pcntrol-zebra", "pcntrol-blnd", "ncntrl-pcr", "ncntrol-xtr", "blank"))
shp_cntrl <- psob_molten %>% filter(SAMPLE.TYPE %in% c( "blank"))
lab_cntrl <- psob_molten %>% filter(SAMPLE.TYPE %in% c( "ncntrl-pcr", "ncntrol-xtr"))
pos_cntrl <- psob_molten %>% filter(SAMPLE.TYPE %in% c( "pcntrol-zebra", "pcntrol-blnd")) 

# plot all controls 
get_default_plot (all_cntrl, taxlev = "PHYLUM", ptitl = "All controls across all locations before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")
get_default_plot (shp_cntrl, taxlev = "PHYLUM", ptitl = "eDNA controls across all locations before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")
get_default_plot (lab_cntrl, taxlev = "PHYLUM", ptitl = "PCR and Extraction controls controls all locations before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")
get_default_plot (pos_cntrl, taxlev = "PHYLUM", ptitl = "Fish positve controls all locations before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")

# describe all controls 
get_molten_ps_description(all_cntrl, rank_level = "SUPERKINGDOM", rank_name = "Eukaryota")
get_molten_ps_description(shp_cntrl, rank_level = "SUPERKINGDOM", rank_name = "Eukaryota")
get_molten_ps_description(lab_cntrl, rank_level = "SUPERKINGDOM", rank_name = "Eukaryota")
get_molten_ps_description(pos_cntrl, rank_level = "SUPERKINGDOM", rank_name = "Eukaryota")


# Removal of contamination with package `decontam`
# ------------------------------------------------
# following:
#   "https://benjjneb.github.io/decontam/vignettes/decontam_intro.html"

# plot initial library sizes
psob_df <- as.data.frame(sample_data(psob)) # Put sample_data into a ggplot-friendly data.frame
psob_df$LibrarySize <- sample_sums(psob)
psob_df <- psob_df[order(psob_df$LibrarySize),]
psob_df$Index <- seq(nrow(psob_df))

ggplot(data = psob_df, aes(x=Index, y=LibrarySize, color=sample.type)) + 
  geom_point() +
  theme_bw()

ggsave("210108_990_r_get_eDNA_phyloseq__decontam_library_coverages.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 50, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# set type for subsequent code -define new variable for sanity reasons 
sample_data(psob)$libconc <- as.numeric(sample_data(psob)$`ng-ul`)


# -- continue here after 20-01-2021 - untested code below-- 


# get combined factor identifying negative controls
sample_data(psob)$is.neg <- sample_data(psob)$`sample-type` %in% c("ncntrl-pcr", "ncntrol-xtr")

# accomodating package decontam - negative controls can't have 0 concentration
as_tibble(sample_data(psob)) %>% filter(is.neg == TRUE) %>% select(is.neg, libconc) %>% print(n = Inf)

sum(sample_data(psob)$is.neg)
sample_data(psob)$is.neg [which(sample_data(psob)$is.neg == TRUE & sample_data(psob)$libconc == 0)] <- FALSE
sum(sample_data(psob)$is.neg)

as_tibble(sample_data(psob)) %>% filter(is.neg == TRUE) %>% select(is.neg, libconc) %>% print(n = Inf)

# frequency and prevalence based contamination identification
contamdf.freq <- isContaminant(psob, method="combined", conc="libconc", neg="is.neg", threshold=0.6)
table(contamdf.freq$contaminant) # FALSE 8917  TRUE 1635 
head(which(contamdf.freq$contaminant))

# randomwly inspecting 81 ASV
set.seed(100)
plot_frequency(psob, taxa_names(psob)[sample(which(contamdf.freq$contaminant),81)], conc="RibLibConcAvg") +
    xlab("DNA Concentration (PicoGreen fluorescent intensity)")

ggsave("200814_decontam_81_likely-contaminats.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_pcm_eukaryotes/Manuscript/200622_display_item_development",
         scale = 3, width = 200, height = 200, units = c("mm"),
         dpi = 500, limitsize = TRUE)

ggsave("200814_decontam_81_likely-contaminats.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/ProcessingPlots",
         scale = 3, width = 200, height = 200, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# "In this plot the dashed black line shows the model of a noncontaminant
# sequence feature for which frequency is expected to be independent of the input
# DNA concentration. The red line shows the model of a contaminant sequence
# feature, for which frequency is expected to be inversely proportional to input
# DNA concentration, as contaminating DNA will make up a larger fraction of the
# total DNA in samples with very little total DNA. Clearly Seq3 fits the red
# contaminant model very well, while Seq1 does not."

# Make phyloseq object of presence-absence in negative controls and true samples
psob.pa <- transform_sample_counts(psob, function(abund) 1*(abund>0))
psob.pa.neg <- prune_samples(sample_data(psob.pa)$Description %in% c("ACNTRL", "XCNTRL"), psob.pa)
psob.pa.pos <- prune_samples(sample_data(psob.pa)$Description == "PCMNT", psob.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos = taxa_sums(psob.pa.pos), pa.neg = taxa_sums(psob.pa.neg),
                      contaminant=contamdf.freq$contaminant)

# look at the incidences of taxa observations in negative controls and positive samples.
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("prevalence (PCR and Extraction Controls)") +
  ylab("prevalence (PCM Samples)") +
  theme_bw()

ggsave("200814_decontam_contaminant_detection_check.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_pcm_eukaryotes/Manuscript/200622_display_item_development",
         scale = 3, width = 50, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

ggsave("200814_decontam_contaminant_detection_check.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/ProcessingPlots",
         scale = 3, width = 50, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# contamination filtering
psob_decontam <- prune_taxa(!contamdf.freq$contaminant, psob)


# -- continue here after 20-01-2021 - untested code above -- 







# remove temporary object
rm(psob_molten)




# IV. Retain Chordates only - and summarise 
# =========================================



# V. Retain fish only - and summarise 
# =====================================





# --- after 12-01-2021 --- 

# Numerical summaries
# -------------------
length(unique(psob_raw_molten$Sample)) # 154 samples from study sites
unique(psob_raw_molten$Sample) %>% grepl(".MM", . , fixed = TRUE) %>% sum # 26 sample Mount Menzies
unique(psob_raw_molten$Sample) %>% grepl(".ME", . , fixed = TRUE) %>% sum # 70 samples Mawson Escarpment
unique(psob_raw_molten$Sample) %>% grepl(".LT", . , fixed = TRUE) %>% sum # 58 samples Lake Terrasovoje
sum(26 + 70 + 58) # 154 - ok

# get total seqencing effort
sum(psob_raw_molten$Abundance) #  16 524 031 sequences total - after import filtering


# A tibble: 1 x 13
#     OTU Abundance Sample BarcodeSequence Location Description superkingdom phylum class order family genus species
#   <int>     <int>  <int>           <int>    <int>       <int>        <int>  <int> <int> <int>  <int> <int>   <int>
# 1  8097      2847    154             154        3           1            5     47   131   342    621  1161    1904

# ASV counts unfiltered data
length(unique(psob_raw_molten$OTU)) # 8097
# count eukaryotes and non-eukaryotes 
psob_raw_molten %>% filter(SUPERKINGDOM %in% c("Eukaryota")) %>% distinct(OTU)  # 2,656 Eukaryota ASV
psob_raw_molten %>% filter(SUPERKINGDOM %!in% c("Eukaryota")) %>% distinct(OTU) # 5,441 non-Eukaryota ASV

# Analyze coverages per samples
coverage_per_sample <- aggregate(psob_raw_molten$Abundance, by=list(Sample=psob_raw_molten$Sample), FUN=sum)
summary(coverage_per_sample)
coverage_per_sample %>% filter(., grepl(".MM", Sample , fixed = TRUE)) %>% summary(x) 
coverage_per_sample %>% filter(., grepl(".ME", Sample , fixed = TRUE)) %>% summary(x) 
coverage_per_sample %>% filter(., grepl(".LT", Sample , fixed = TRUE)) %>% summary(x) 

# get coverages per ASV and analyze
coverage_per_asv <- aggregate(psob_raw_molten$Abundance, by=list(ASV=psob_raw_molten$OTU), FUN=sum)
coverage_per_asv <- coverage_per_asv %>% arrange(desc(x))
summary(coverage_per_asv)


# get most species ordered by frequency
psob_asv_list <- left_join(coverage_per_asv , psob_raw_molten, by = c("ASV" = "OTU")) %>%
    distinct_at(vars("ASV", "X", "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")) 

psob_asv_list %>% head(., n = 12)

psob_asv_list %>% 
  arrange(superkingdom, phylum, class, order, family, genus, species) %>%
  write.xlsx(.,"/Users/paul/Documents/OU_pcm_eukaryotes/Manuscript/200622_display_item_development/200814_all_unfiltered_phyla_at_all_locations.xlsx", overwrite = TRUE)

psob_asv_list %>% 
  arrange(superkingdom, phylum, class, order, family, genus, species) %>%
  write.xlsx(.,"/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/SpeciesLists/200814_all_unfiltered_phyla_at_all_locations.xlsx", overwrite = TRUE)

# save or load molten state 
save.image(file = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/R/200_all_data_long_export_raw-image.Rdata")
save(psob_raw_molten, file = "/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/R/200_all_data_long_export_raw.Rdata")
load("/Users/paul/Documents/OU_pcm_eukaryotes/Zenodo/R/200_all_data_long_export_raw.Rdata")

