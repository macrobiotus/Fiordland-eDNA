# **********************************************
# * Create, filter, and write Physloseq object *
# **********************************************
# 12-Dec-2020, 08-Jan-2020, 11-Jan-2020, 20-Jan-2020, 21-Jan-2020,  10-Feb-2020 


# Packages loading 
# ==================
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse") # work using tibbles
library("magrittr") # more pipes
library("phyloseq") # handle Qiime 2 data in R - best used to generate long dataframe

library("data.table") # faster handling of large tables - such as molten Phyloseq objects
library("future.apply") # faster handling of large tables - such as molten Phyloseq objects

library("decontam") # decontamination - check `https://benjjneb.github.io/decontam/vignettes/decontam_intro.html`
library("openxlsx") # write Excel tables



# Functions
# =========

# Define new operator "not in"
"%!in%" <- function(x, y) !(x %in% y)

# Create Phyloseq object.
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

# This function removes "0" count phylotypes from samples and samples with "0" phylotypes.
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
  
  molten_phyloseq_object <- molten_phyloseq_object %>% mutate(across(c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), .funs = gsub, pattern = "__", replacement = " "))

  molten_phyloseq_object <- molten_phyloseq_object %>% mutate(across(c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), .funs = gsub, pattern = "_", replacement = " "))

  molten_phyloseq_object <- molten_phyloseq_object %>% mutate(across(c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), .funs = gsub, pattern = "  ", replacement = " "))

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


# Get an easily digestible coverage plot of a molten Phyloseq object (e.g. before, during, after filtering)
# ----------------------------------------------------------------------------------------------------

get_default_coverage_plot = function (psob_molten, taxlev, taxlev_fill = taxlev, ptitl, pxlab, pylab, facet_var = "LOC.NAME" ){

  require(ggplot2)
  require(data.table)

  # melting Phyloseq object to data table for merging and speed
  psob_molten <- data.table(psob_molten)
 
  # set sorting key properly
  setkey(psob_molten,ASV)
 
  
  # aggregate on chosen level
  #   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
  psob_molten <- psob_molten[, lapply(.SD, sum, na.rm=TRUE), by=c(facet_var, "ASV", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS",  "SPECIES"), .SDcols=c("ABUNDANCE") ]

  #  resort for clarity
  keycol <-c("ASV",facet_var)
  setorderv(psob_molten, keycol)
  
  # debugging
  print(head(psob_molten))

  # add presence-absence abundance column
  psob_molten <- psob_molten[ , ASVPRESENT :=  fifelse(ABUNDANCE == 0 , 0, 1, na=NA)]
  
   # debugging
  print(head(psob_molten))
    
  ggplot(psob_molten, aes_string(x = taxlev, y = "ABUNDANCE", fill = taxlev_fill)) +
    geom_bar(stat = "identity", position = "stack", colour = NA, size=0) +
    facet_grid(get(facet_var) ~ ., shrink = TRUE, scales = "fixed") +
    theme_bw() +
    theme(legend.position = "none") +
    theme(strip.text.y = element_text(angle=0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7), 
          axis.ticks.y = element_blank()) +
    labs( title = ptitl) + xlab(pxlab) + ylab(pylab)

}

# Get an easily digestible occurence plot of a molten Phyloseq object (e.g. before, during, after filtering)
# ----------------------------------------------------------------------------------------------------

get_default_ocurrence_plot = function (psob_molten, taxlev, taxlev_fill = taxlev, ptitl, pxlab, pylab, facet_var = "LOC.NAME" ){
  
  require(ggplot2)
  require(data.table)

  # melting Phyloseq object to data table for merging and speed
  psob_molten <- data.table(psob_molten)
 
  # set sorting key properly
  setkey(psob_molten,ASV)
 
  
  # aggregate on chosen level
  #   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
  psob_molten <- psob_molten[, lapply(.SD, sum, na.rm=TRUE), by=c(facet_var, "ASV", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS",  "SPECIES"), .SDcols=c("ABUNDANCE") ]

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
    facet_grid(get(facet_var) ~ ., shrink = TRUE, scales = "fixed") +
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

get_molten_ps_description = function(ps, rank_level = "SUPERKINGDOM", rank_name = "Eukaryota", prnt_n = Inf){

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
  cov_s_sd <- sd(coverage_per_sample$x)
  
  # get ASV coverages
  coverage_per_asv <- aggregate(ps$ABUNDANCE, by=list(ASV=ps$ASV), FUN=sum)
  
  # print(head(coverage_per_asv))
  summary(coverage_per_asv)
  cov_a_min  <- summary(coverage_per_asv)[c(1) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_a_med  <- summary(coverage_per_asv)[c(3) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_a_mean <- summary(coverage_per_asv)[c(4) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_a_max <- summary(coverage_per_asv)[c(6) ,2] %>% str_squish() %>% gsub("[^0-9]", "", .)
  cov_a_sd <- sd(coverage_per_asv$x)
  
  # print summary text
  print(paste0("The current data set contains ", f_count, " sequences across ", s_count, " samples and ",   
                 a_count, " ASV's (", ss_count," ", ss_string,
                 ", as well as ", nss_count, " non-", ss_string, ", i.e. ", nss_string,
                 "). Sample mean (min., med., max.) coverage is ", cov_s_mean, 
                 " reads (", cov_s_min, ", ", cov_s_med, ", ", cov_s_max, 
                 "), and ASV mean (min, median, max) coverage is ", cov_a_mean, 
                 " reads (", cov_a_min, ", ", cov_a_med, ", ", cov_a_max, ")."
                 ))
  print(paste0("Sample coverage sd is: ", cov_s_sd))
  print(paste0("ASV coverage sd is:", cov_a_sd))

  # summarize distinct values across the long data frame
  print("In the following summary variable names shown are hard-coded in function \"get_molten_ps_stats\": ") 
  show_vars <- c("ASV", "ABUNDANCE", "SAMPLE", "LOC.NAME", "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")
  ps %>% dplyr::select(any_of(show_vars)) %>% summarize_all(n_distinct, na.rm = TRUE) %>% print()
  
  # get most ASV's (not species) ordered by frequency
  print("Getting and returning highest-covered ASV's (not: species) ordered by frequency:")
  ps_asv_list <- left_join(coverage_per_asv , ps, by = c("ASV" = "ASV")) %>%
    distinct_at(vars("ASV", "x", "SUPERKINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"))
    
  ps_asv_list %>% arrange(desc(x)) %>% head(., n = prnt_n) %>% print()

  # get most spcies's (not ASV's) ordered by frequency
  print("Getting and returning highest-covered species (not: ASV's) ordered by frequency:")
  
  ps_species_list <- aggregate(ps$ABUNDANCE, by=list(PHYLUM=ps$PHYLUM, CLASS=ps$CLASS, ORDER=ps$ORDER, FAMILY=ps$FAMILY, GENUS=ps$GENUS, SPECIES=ps$SPECIES), FUN=sum) 
  ps_species_list %>% arrange(desc(x)) %>% head(., n = prnt_n) %>% print()

}

# get any variable in plate format 
# ---------------------------------

show_plate_loading = function(psob_molten, level = "PHYLUM", tax_vector = c("Chordata"), ggplx = "PHYLUM", ggply = "ABUNDANCE", ggfill = "PHYLUM", ggpltxt = "SAMPLE.TYPE"){
  
  require("data.table")
  
  message(paste("\"level\"=", level))
  message(paste("\"tax_vector\"=", tax_vector))
  
  # melting Phyloseq object to data table for merging and speed
  psob_dt <- data.table(psob_molten)
 
  # set sorting key properly
  setkey(psob_dt,KEY)
 
  # aggregate selected variables
  #   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
  #   for plate position (KEY) get
  #   ASV, each summed ABUNDANCE
  #   ASV, TAXONOMY
  #   SAMPLE.NAME
  #   SAMPLE.TYPE
  psob_dt <- psob_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("KEY", "PLATE", "ROW", "COL", "SAMPLE.NAME", "SAMPLE.TYPE", "ASV", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS",  "SPECIES"), .SDcols=c("ABUNDANCE") ]

  #  resort for clarity
  keycol <-c("KEY")
  setorderv(psob_dt, keycol)
  
  # add presence-absence abundance column
  psob_dt <- psob_dt[ , ASVPRESENT :=  fifelse(ABUNDANCE == 0 , 0, 1, na=NA)]

  # keep only chordates
  psob_dt <- psob_dt[get(level) %in% tax_vector ]
  
  # work around ggplot2 idiosyncrasies
  setnames(psob_dt, "COL", "CL")
  setnames(psob_dt, "ROW", "RW")

  psob_dt$CL <- factor(psob_dt$CL, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  psob_dt$PLATE <- factor(psob_dt$PLATE, levels = c("U.1", "U.2", "E.1", "E.2"))
  
  # show available vars
  message("The following variables can be plotted:")
  print(head(psob_dt))

  # plot out 
  ggplot(psob_dt, aes_string(x = ggplx, y = ggply, fill = ggfill)) +
    geom_bar(stat = "identity", position = "stack", colour = NA, size=0) +
    facet_grid(RW ~ PLATE + CL , shrink = TRUE, scales = "fixed") +
    geom_text(x = 1, y = 3, aes(label = get(ggpltxt)), data = psob_dt, size = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(strip.text.y = element_text(angle=0)) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7), 
          axis.ticks.y = element_blank()) +
    labs( title = "Plate overview") + xlab("plate column") + ylab("plate row") %>% quartz()
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

sample_data(psob)$`sample-type`

# more filtering determined as ok previously
#   remove zebra fish control - appears to have been forgotten to be pipetted in
#   Package decontam works only with defined data, remove undefined data from PHYLOSEQ object
psob <- psob %>% subset_samples(., sample.type != "empty" ) %>% subset_samples(., ng.ul != "0" ) %>% subset_samples(., sample.type != "pcntrol-zebra")

# Save or load initial state after import
# ---------------------------------------
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210203_990_r_get_eDNA_phyloseq__import_workspace.Rdata")
save(psob, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210203_990_r_get_eDNA_phyloseq__imported-psob.Rdata")
save(psob, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210210_990_r_get_eDNA_phyloseq__imported-psob.Rds")


# II. Inspect and summarize raw data
# ==================================

# melt and tidy Phyloseq object 
(psob_molten <- get_tidy_molten_ps(psob))

# Check composition and ASV count of entire data prior to cleanup
# ---------------------------------------------------------------

# plot composition per sample type
get_default_coverage_plot (psob_molten, facet_var = "SAMPLE.TYPE", taxlev = "PHYLUM", ptitl = "Phyla across all sample types before filtering", pxlab = "phyla (NCBI taxonomy) ", pylab =  "read counts in sample category (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered-sample-type-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_ocurrence_plot (psob_molten, facet_var = "SAMPLE.TYPE", taxlev = "PHYLUM", ptitl = "Phyla across all sample types before filtering", pxlab = "phyla (NCBI taxonomy) ", pylab =  "ASV counts in sample category (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered-sample-type-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# plot composition per location
get_default_coverage_plot (psob_molten, taxlev = "PHYLUM", ptitl = "Phyla across all locations before filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "read counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered-location-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_ocurrence_plot (psob_molten, taxlev = "PHYLUM", ptitl = "Phyla across all locations before filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "ASV counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered-location-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# written summaries
capture.output(get_molten_ps_description(psob_molten) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-unfiltered_summary.txt")

# plate loadings - chordates only otherwise plotter crashes
# - asvs

show_plate_loading(psob_molten,  ggply = "ASVPRESENT")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered-plate-loading-asv-counts_chordates_only.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# - read counts
show_plate_loading(psob_molten,  ggply = "ABUNDANCE")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered-plate-loading-read-counts_chordates_only.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# III. Inspect and summarize control data
# =======================================

# melt and tidy Phyloseq object 
(psob_molten <- get_tidy_molten_ps(psob))

# check sample types - what sample type are available? 
unique(psob_molten$SAMPLE.TYPE) # "eDNA": eDNA sample
                                # "blank": eDNA collection blank - bottel with mol grade water
                                # "pcntrol-blnd" positive controls (fish blend)
                                # "ncntrl-pcr": amplification blank
                                # "ncntrol-xtr"
                                # "empty": unassigned wells  

# split controls -  how many of each sample type are available? 
psob_molten %>% group_by(KEY) %>% summarise(KEY = KEY, CNTROL = SAMPLE.TYPE ) %>% distinct() %>% pull("CNTROL") %>% table()
#  blank         eDNA        empty   ncntrl-pcr  ncntrol-xtr pcntrol-blnd 
#    33          101            3           12            9            2 

# split data types and print text summaries
# ------------------------------------------
dna_cntrl <- psob_molten %>% filter(SAMPLE.TYPE %in% c("blank"))
pos_cntrl <- psob_molten %>% filter(SAMPLE.TYPE %in% c("pcntrol-blnd"))
neg_cntrl <- psob_molten %>% filter(SAMPLE.TYPE %in% c("ncntrl-pcr", "ncntrol-xtr"))
dna_smpls <- psob_molten %>% filter(SAMPLE.TYPE %in% c("eDNA"))

# some cross contamination - but really low, should affect results if removed check summaries
capture.output(get_molten_ps_description(dna_cntrl) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-dna_cntrl_summary.txt")
capture.output(get_molten_ps_description(pos_cntrl) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-pos_cntrl_summary.txt")
capture.output(get_molten_ps_description(neg_cntrl) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-neg_cntrl_summary.txt")
capture.output(get_molten_ps_description(dna_smpls) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-dna_smpls_summary.txt")

# Check cross contamination - visualize positive control data
# ----------------------------------------------------------

# keep read counts from  positive controls
# Species likley from controls
#   in controls: Corydoras sterbai, Chromobotia macracanthus, Poecilia reticulata, Poecilia sphenops, and Danio rerio
#   in data :    Crydoras aeneus,   Chromobotia macracanthus, Poecilia reticulata, Poecilia formosa
psob_molten_fish_controls <- get_tidy_molten_ps(psob) %>% mutate(ABUNDANCE = ifelse(GENUS %in% c("Crydoras", "Chromobotia","Poecilia"), ABUNDANCE, 0))

# show distribution of positive controls across plates
# reads only visible in the right places at this scale
show_plate_loading(psob_molten_fish_controls,  ggply = "ABUNDANCE")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-pcontrol-plate-loading-read-counts_chordates_only.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

show_plate_loading(psob_molten_fish_controls,  ggply = "ASVPRESENT")
# some cross contamination in 3 samples
ggsave("210211_990_r_get_eDNA_phyloseq__psob-pcontrol-plate-loading-asv-counts_chordates_only.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# check read counts in U.2.A.2 and U.2.A.3 and (U.2.D.3)
sort(unique(psob_molten_fish_controls$SAMPLE)) #i.e. U-2-A-2-PC191218-AZ, U-2-A-3-PC191218-AYAZ-blk, U-2-D-3-ncntrl-pcr
# they are quite low 
cross_cont <- psob_molten_fish_controls %>% filter(SAMPLE %in% c("U-2-A-2-PC191218-AZ", "U-2-A-3-PC191218-AYAZ-blk", "U-2-D-3-ncntrl-pcr")) %>%
  filter(ABUNDANCE != 0) %>% unique() %>% select(SAMPLE, ABUNDANCE)

# from the observed cross contamination get a simple model of how mmany reads to remove
# assume Poisson distribution of cross contamination counts
lmbda <- MASS::fitdistr(cross_cont %>% pull(ABUNDANCE), "Poisson")

# find upper read count limit of 95% of typical cross-contamination as observed
treshhold <- qpois(0.95, lmbda$estimate) # will remove 18S with 18 reads or less as possibly cross-contaminated

rm(psob_molten_fish_controls)


# IV. Check  cross contamination with package `decontam`
# ======================================================
# following:
#   "https://benjjneb.github.io/decontam/vignettes/decontam_intro.html"

# Summary for decontam package input data is not necessary anymore as done above already 
# copy input data for sanity reasons
psob_input <- psob 

# initial plotting and data preparation
# -------------------------------------

# plot initial library sizes
psob_df <- as.data.frame(sample_data(psob_input)) # Put sample_data into a ggplot-friendly data.frame
psob_df$LibrarySize <- sample_sums(psob_input)
psob_df <- psob_df[order(psob_df$LibrarySize),]
psob_df$Index <- seq(nrow(psob_df))

ggplot(data = psob_df, aes(x=Index, y=LibrarySize, color=sample.type)) + geom_point() + theme_bw() + 
    theme(strip.text.y = element_text(angle=0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1,  size = 7), 
          axis.ticks.y = element_blank()) +
          theme(legend.title=element_blank()) +
          labs(title = "Read Counts per library, and library types") +
          xlab("libraries, sorted by size") + 
          ylab("read counts")

ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered_library_size_per_sample_type.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)



# Identify  contamination  based on read frequency
# -------------------------------------------------

# set type for subsequent code 
sample_data(psob_input)$ng.ul <- as.numeric(sample_data(psob_input)$ng.ul)

# get combined factor identifying negative controls
sample_data(psob_input)$is.neg <- sample_data(psob_input)$sample.type %in% c("ncntrl-pcr", "ncntrol-xtr", "blank")

# frequency and prevalence based contamination identification
contamdf.freq <- isContaminant(psob_input, method="combined", conc="ng.ul", neg="is.neg", threshold=0.6)
table(contamdf.freq$contaminant) # FALSE 2057 TRUE 114 
which(contamdf.freq$contaminant) # indices for contaminant ASVs

# inspecting and plotting contaminant ASV
plot_frequency(psob_input, taxa_names(psob)[which(contamdf.freq$contaminant)], conc="ng.ul") +
    xlab("DNA Concentration (PicoGreen fluorescent intensity)")

ggsave("210211_990_r_get_eDNA_phyloseq__psob-unfiltered_likely_contaminants.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 250, height = 125, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# "In this plot the dashed black line shows the model of a noncontaminant
# sequence feature for which frequency is expected to be independent of the input
# DNA concentration. The red line shows the model of a contaminant sequence
# feature, for which frequency is expected to be inversely proportional to input
# DNA concentration, as contaminating DNA will make up a larger fraction of the
# total DNA in samples with very little total DNA. Clearly Seq3 fits the red
# contaminant model very well, while Seq1 does not."


# further inspect likely contaminants
# -----------------------------------

# separate assumed contamination from assumed true signals as phyloseq object 
psob_cand_cont <- prune_taxa(contamdf.freq$contaminant, psob_input)
psob_cand_eDNA <- prune_taxa(!contamdf.freq$contaminant, psob_input)

# melt both 
psob_cont_molten <- get_tidy_molten_ps(psob_cand_cont)
psob_eDNA_molten <- get_tidy_molten_ps(psob_cand_eDNA)

# just keep the important taxa for inspection  
psob_cont_molten_chrd <- psob_cont_molten %>% mutate(ABUNDANCE = ifelse(PHYLUM %in% c("Chordata"), ABUNDANCE, 0))
psob_eDNA_molten_chrd <- psob_eDNA_molten %>% mutate(ABUNDANCE = ifelse(PHYLUM %in% c("Chordata"), ABUNDANCE, 0))

# get text summaries
capture.output(get_molten_ps_description(psob_cont_molten_chrd) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-dna_cont_summary.txt")
capture.output(get_molten_ps_description(psob_eDNA_molten_chrd) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-dna_eDNA_summary.txt")

# check plate loadings of both categories
show_plate_loading(psob_cont_molten_chrd,  ggply = "ABUNDANCE")
show_plate_loading(psob_cont_molten_chrd,  ggply = "ASVPRESENT")

show_plate_loading(psob_eDNA_molten_chrd,  ggply = "ASVPRESENT")
show_plate_loading(psob_eDNA_molten_chrd,  ggply = "ABUNDANCE")

# decontam doesn't work well with the data at hand as contaminants are still in the data and mixed with obviously real samples
# resorting to simple subtraction

# V. Remove contamination with subtractive approach 
# ===================================================

# melt and tidy Phyloseq object 
(psob_molten <- get_tidy_molten_ps(psob))

# check sample types - what sample type are available? 
unique(psob_molten$SAMPLE.TYPE) # "eDNA": eDNA sample
                                # "blank": eDNA collection blank - bottel with mol grade water
                                # "pcntrol-blnd" positive controls (fish blend)
                                # "ncntrl-pcr": amplification blank
                                # "ncntrol-xtr"

# remove all species and ASVs found in all controls 
# isolate all control data in unison 
psob_molten_all_controls <- psob_molten %>% filter(SAMPLE.TYPE %in% c("pcntrol-blnd", "ncntrl-pcr", "blank", "ncntrol-xtr"))
capture.output(get_molten_ps_description(psob_molten_all_controls) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-all_controls_summary.txt")


# pull out all species and asvs that have to be removed from everywhere
cont_asvs <- psob_molten_all_controls %>% pull(ASV) %>% unique()
cont_spec <- psob_molten_all_controls %>% pull(SPECIES) %>% unique()

# isolating and inspecting all contamination
psob_molten_all_contamination <-  psob_molten %>% filter(ASV %in% cont_asvs | SPECIES %in% cont_spec[1:9])

show_plate_loading(psob_molten_all_contamination,  ggply = "ASVPRESENT")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-contamination-plate-loading-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

show_plate_loading(psob_molten_all_contamination,  ggply = "ABUNDANCE")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-contamination-plate-loading-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_coverage_plot (psob_molten_all_contamination, taxlev = "PHYLUM", ptitl = "Contamination across all locations before filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "read counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-contamination-location-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 30, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_ocurrence_plot (psob_molten_all_contamination, taxlev = "PHYLUM", ptitl = "Contamination across all locations before filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "ASV counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-contamination-location-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 30, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# subtract contamination and inspect again 
psob_molten_clean <- anti_join(psob_molten, psob_molten_all_contamination, by = "ASV", copy = FALSE)

get_default_coverage_plot (psob_molten_clean, taxlev = "PHYLUM", ptitl = "Phyla across all locations after filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "read counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-eDNA-location-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 90, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_ocurrence_plot (psob_molten_clean, taxlev = "PHYLUM", ptitl = "Phyla across all locations after filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "ASV counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-eDNA-location-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 90, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

capture.output(get_molten_ps_description(psob_molten_clean) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-all_data_summary.txt")



# remove low abundance reads

#   check the lower end of the read count distribution 
summary(psob_molten_clean$ABUNDANCE)
hist(psob_molten_clean$ABUNDANCE, breaks = 300000, xlim = c(0,25))

#   filter out ASV with less then 15 reads
possible_cc <- psob_molten_clean %>% filter(ABUNDANCE <= treshhold) %>% filter(PHYLUM == "Chordata")

#   low abundance filtering gets rid of (some(!!!), not necessarily all)
#   Trichosurus vulpecula
#   Sus scrofa
get_molten_ps_description(possible_cc) 
capture.output(get_molten_ps_description(possible_cc) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__low_abundance_reads__summary.txt")


# subtract possible cross contamination (by abundance) and inspect again 
psob_molten_clean_new <- anti_join(psob_molten_clean, possible_cc, by = "ASV", copy = FALSE)


# get a full species list with Chordates at the top
psob_molten_clean_chordates_top <- psob_molten_clean_new %>% mutate(ABUNDANCE = ifelse(PHYLUM %in% c("Chordata"), ABUNDANCE, 0))

capture.output(get_molten_ps_description(psob_molten_clean_chordates_top) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-all_data_and_chordates_summary.txt")

# isolate fish and marine mammals
psob_molten_clean_marine <- psob_molten_clean_chordates_top %>% filter(
  CLASS %in% c("Actinopteri", "Chondrichthyes") |
  ORDER %in% c("Cetacea") | 
  FAMILY %in% c("Otariidae")) %>% filter(!(GENUS %in% c("Sardinops")))

capture.output(get_molten_ps_description(psob_molten_clean_marine) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-clean_marine_summary.txt")

# VII. Get abundance values for biologically replicate observations in clean data 
# ===============================================================================

# get a suitable grouping variable in the clean data
# ---------------------------------------------------

# define replicates with use of external table
pth <- c("/Users/paul/Documents/OU_eDNA/191213_field_work/201130_sample_overview_updated.xlsx")
sample_order <- readxl::read_excel(pth, range = "F1:F91", .name_repair = "universal")
sample_order %<>% rowid_to_column("ID") %>% print(n = Inf)
sample_order$RepID <- rep(c(1,2,3), times = 2, length.out = length(sample_order$ID), NA, each = 1)
sample_order %<>% group_by(RepID) %>% mutate(SetID = row_number(RepID))
sample_order %<>% filter(!grepl("blk", sample_name)) 
sample_order %>% print(n = Inf)

# merge set IDs with clean data - get overview of sample covergae and test merging of ids 
overview <- full_join(as_tibble(psob_molten_clean_marine$SAMPLE.NAME), as_tibble(sample_order$sample_name), keep = TRUE) %>% 
  distinct() %>% setnames(c("mdata", "tdata")) %>% arrange(tdata) %>% print(n = Inf)

smpl_coverage <- left_join( overview %>% dplyr::select(tdata), 
                            psob_molten_clean_marine %>% dplyr::select(SAMPLE.NAME, LOC.NAME,  INSIDE.RESERVE),
                            by=c("tdata" = "SAMPLE.NAME"),
                            keep=TRUE) %>% 
                            distinct() %>% 
                            arrange(tdata) %>% 
                            print(n = Inf)

# write sample overview
# ---------------------
# save overview for writing and if needed later
write.xlsx(smpl_coverage, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210211_990_r_get_eDNA_phyloseq__eDNA_sampling_success.xlsx", asTable = FALSE)
save(smpl_coverage, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210211_990_r_get_eDNA_phyloseq__eDNA_sampling_success.Rds")


# merge  replicate water samples
# ------------------------------

# copy grouping variable to main object
clean_marine <- left_join(psob_molten_clean_marine, sample_order,  by=c("SAMPLE.NAME" = "sample_name"))


# test sample merging approach - continue here after 12-Feb-2020
#  1) create test data with 344 rows
foo <- clean_marine[ c("ASV", "SAMPLE.NAME", "ABUNDANCE", "ID", "SetID", "RepID")] %>% arrange(ASV, SetID, RepID) %>% print(n = 20)

#  2) erases possibly unique metadata - to yield 251 rows
foo %>% group_by(ASV, SetID) %>% summarise(SET_ABUNDANCE = sum(ABUNDANCE)) 

#  3) keeps all 344 rows and duplicates SET_ABUNDANCE with each group of ASV and SetID
foo %>% group_by(ASV, SetID) %>% mutate(SET_ABUNDANCE = sum(ABUNDANCE))
foo %>% group_by(ASV, SetID) %>% mutate(SET_ABUNDANCE = sum(ABUNDANCE)) %>% pull(SET_ABUNDANCE)

#  keeping all data and duplicated SET_ABUNDANCE values with each group of ASV and SetID - pull out later when required
clean_marine <- clean_marine %>% group_by(ASV, SetID) %>% mutate(SET_ABUNDANCE = sum(ABUNDANCE))


# show plate loading 
show_plate_loading(psob_molten_clean_marine,  ggply = "ABUNDANCE")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-clean_marine-plate-loading-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# plot composition per location
get_default_coverage_plot (psob_molten_clean_marine, taxlev = "FAMILY", ptitl = "Marine families across all locations (read counts)", pxlab = "families (NCBI taxonomy)", pylab =  "read counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-clean_marine-location-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 120, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_ocurrence_plot (psob_molten_clean_marine, taxlev = "FAMILY", ptitl = "Marine families across all locations (sequence variants)", pxlab = "families (NCBI taxonomy)", pylab =  "ASV counts at each location (y scales fixed)")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-clean_marine-location-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 120, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

show_plate_loading(psob_molten_clean_marine,  ggply = "ASVPRESENT")
ggsave("210211_990_r_get_eDNA_phyloseq__psob-clean_marine-plate-loading-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# get text summary - also show alignments judgments
capture.output(get_molten_ps_description(psob_molten_clean_marine) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/210211_990_r_get_eDNA_phyloseq__psob-clean_marine_summary.txt")


# VII. Save eDNA object for further analysis
# ==========================================

# save or load molten state
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210108_990_r_get_eDNA_phyloseq__export_workspace.Rdata")
save(psob_molten_clean_marine, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210210_990_r_get_eDNA_phyloseq__clean-marine-eDNA.Rds")
save(psob_molten_clean_marine, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210210_990_r_get_eDNA_phyloseq__clean-marine-eDNA.Rds")
