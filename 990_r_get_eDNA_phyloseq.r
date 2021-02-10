# **********************************************
# * Create, filter, and write Physloseq object *
# **********************************************
# 12-Dec-2020, 08-Jan-2020, 11-Jan-2020, 20-Jan-2020, 21-Jan-2020,  10-Feb-2020 


# Packages loading 
# ==================
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse") # work using tibbles
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
  ps %>% select(any_of(show_vars)) %>% summarize_all(n_distinct, na.rm = TRUE) %>% print()
  
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
ggsave("211002_990_r_get_eDNA_phyloseq__psob-unfiltered-sample-type-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_ocurrence_plot (psob_molten, facet_var = "SAMPLE.TYPE", taxlev = "PHYLUM", ptitl = "Phyla across all sample types before filtering", pxlab = "phyla (NCBI taxonomy) ", pylab =  "ASV counts in sample category (y scales fixed)")
ggsave("211002_990_r_get_eDNA_phyloseq__psob-unfiltered-sample-type-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# plot composition per location
get_default_coverage_plot (psob_molten, taxlev = "PHYLUM", ptitl = "Phyla across all locations before filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "read counts at each location (y scales fixed)")
ggsave("211002_990_r_get_eDNA_phyloseq__psob-unfiltered-location-read-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_ocurrence_plot (psob_molten, taxlev = "PHYLUM", ptitl = "Phyla across all locations before filtering", pxlab = "phyla (NCBI taxonomy)", pylab =  "ASV counts at each location (y scales fixed)")
ggsave("211002_990_r_get_eDNA_phyloseq__psob-unfiltered-location-asv-counts.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# written summaries
capture.output(get_molten_ps_description(psob_molten) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/211002_990_r_get_eDNA_phyloseq__psob-unfiltered_summary.txt")

# plate loadings - chordates only otherwise plotter crashes
# - asvs

show_plate_loading(psob_molten,  ggply = "ASVPRESENT")
ggsave("211002_990_r_get_eDNA_phyloseq__psob-unfiltered-plate-loading-asv-counts_chordates_only.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# - read counts
show_plate_loading(psob_molten,  ggply = "ABUNDANCE")
ggsave("211002_990_r_get_eDNA_phyloseq__psob-unfiltered-plate-loading-read-counts_chordates_only.pdf", plot = last_plot(), 
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
capture.output(get_molten_ps_description(dna_cntrl) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/211002_990_r_get_eDNA_phyloseq__psob-dna_cntrl_summary.txt")
capture.output(get_molten_ps_description(pos_cntrl) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/211002_990_r_get_eDNA_phyloseq__psob-pos_cntrl_summary.txt")
capture.output(get_molten_ps_description(neg_cntrl) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/211002_990_r_get_eDNA_phyloseq__psob-neg_cntrl_summary.txt")
capture.output(get_molten_ps_description(dna_smpls) , file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/text_summaries/211002_990_r_get_eDNA_phyloseq__psob-dna_smpls_summary.txt")

# Check cross contamination - visualize positive control data
# ----------------------------------------------------------

# keep read counts from  positive controls
# Species likley from controls
#   in controls: Corydoras sterbai, Chromobotia macracanthus, Poecilia reticulata, Poecilia sphenops, and Danio rerio
#   in data :    Crydoras aeneus,   Chromobotia macracanthus, Poecilia reticulata, Poecilia formosa
psob_molten_fish_controls <- get_tidy_molten_ps(psob) %>% mutate(ABUNDANCE = ifelse(GENUS %in% c("Crydoras", "Chromobotia","Poecilia"), ABUNDANCE, 0))

# show distribution of positive controls across plates
show_plate_loading(psob_molten_fish_controls,  ggply = "ASVPRESENT")
# some cross contamination in 3 samples
ggsave("211002_990_r_get_eDNA_phyloseq__psob-pcontrol-plate-loading-asv-counts_chordates_only.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# extremely low read count  
show_plate_loading(psob_molten_fish_controls,  ggply = "ABUNDANCE")
ggsave("211002_990_r_get_eDNA_phyloseq__psob-pcontrol-plate-loading-read-counts_chordates_only.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures",
         scale = 3, width = 125, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

rm(psob_molten_fish_controls)

# ------- continue here after 10-02-2020 ------- 

# III. Remove  cross contamination
# =======================================

# Summary for decontam package input data is not necessary anymore as done above already 
# copy input data for sanity reasons
psob_input <- psob 

# Identification contamination with package `decontam`
# --------------------------------------------------------

# following:
#   "https://benjjneb.github.io/decontam/vignettes/decontam_intro.html"








# plot initial library sizes
psob_df <- as.data.frame(sample_data(psob_input)) # Put sample_data into a ggplot-friendly data.frame
psob_df$LibrarySize <- sample_sums(psob_input)
psob_df <- psob_df[order(psob_df$LibrarySize),]
psob_df$Index <- seq(nrow(psob_df))

ggplot(data = psob_df, aes(x=Index, y=LibrarySize, color=sample.type)) + geom_point() + theme_bw()
ggsave("210121_990_r_get_eDNA_phyloseq__library-coverages-by-type.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 50, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# set type for subsequent code 
sample_data(psob_input)$ng.ul <- as.numeric(sample_data(psob_input)$ng.ul)

# get combined factor identifying negative controls
sample_data(psob_input)$is.neg <- sample_data(psob_input)$sample.type %in% c("ncntrl-pcr", "ncntrol-xtr", "blank")

# accomodating package decontam - negative controls can't have 0 concentration
# as_tibble(sample_data(psob)) %>% filter(is.neg == TRUE) %>% select(is.neg, libconc) %>% print(n = Inf)
# sum(sample_data(psob)$is.neg)
# sample_data(psob)$is.neg [which(sample_data(psob)$is.neg == TRUE & sample_data(psob)$libconc == 0)] <- FALSE
# sum(sample_data(psob)$is.neg)
# as_tibble(sample_data(psob)) %>% filter(is.neg == TRUE) %>% select(is.neg, libconc) %>% print(n = Inf)

# frequency and prevalence based contamination identification
contamdf.freq <- isContaminant(psob_input, method="combined", conc="ng.ul", neg="is.neg", threshold=0.6)
table(contamdf.freq$contaminant) # FALSE 2969  TRUE 143 
which(contamdf.freq$contaminant) # indices for contaminant ASVs

# inspecting and plotting contaminant ASV
plot_frequency(psob_input, taxa_names(psob)[which(contamdf.freq$contaminant)], conc="ng.ul") +
    xlab("DNA Concentration (PicoGreen fluorescent intensity)")

ggsave("210121_990_r_get_eDNA_phyloseq__decontam-likely-contaminats.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 200, height = 200, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# "In this plot the dashed black line shows the model of a noncontaminant
# sequence feature for which frequency is expected to be independent of the input
# DNA concentration. The red line shows the model of a contaminant sequence
# feature, for which frequency is expected to be inversely proportional to input
# DNA concentration, as contaminating DNA will make up a larger fraction of the
# total DNA in samples with very little total DNA. Clearly Seq3 fits the red
# contaminant model very well, while Seq1 does not."


# further inspect likely contaminants
psob_contaminants <- prune_taxa(contamdf.freq$contaminant, psob_input)

(psob_molten <- get_tidy_molten_ps(psob_contaminants))
get_default_coverage_plot (psob_molten, facet_var = "SAMPLE.TYPE", taxlev = "PHYLUM", ptitl = "Contaminant phyla across all sample types before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")
get_default_ocurrence_plot (psob_molten, facet_var = "SAMPLE.TYPE", taxlev = "PHYLUM", ptitl = "Contaminant phyla across all sample types before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")
get_molten_ps_description(psob_molten)

show_plate_loading(psob_molten)
ggsave("210121_990_r_get_eDNA_phyloseq__psob-contaminants-sample-type-overview.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_coverage_plot (psob_molten, taxlev = "PHYLUM", ptitl = "Contaminant phyla across all locations before filtering", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")

ggsave("210121_990_r_get_eDNA_phyloseq__psob-unfiltered-location-overview.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_molten_ps_description(psob_molten)

rm(psob_molten)

# Make phyloseq object of presence-absence in negative controls and true samples
psob.pa <- transform_sample_counts(psob_input, function(abund) 1*(abund>0))

# sample.type %in% c("ncntrl-pcr", "ncntrol-xtr", "blank")

psob.pa.neg <- prune_samples(sample_data(psob.pa)$sample.type %in% c("ncntrl-pcr", "ncntrol-xtr", "blank"), psob.pa)
psob.pa.pos <- prune_samples(sample_data(psob.pa)$sample.type %in% c("eDNA", "pcntrol-blnd", "pcntrol-zebra"), psob.pa)



# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos = taxa_sums(psob.pa.pos), pa.neg = taxa_sums(psob.pa.neg), contaminant=contamdf.freq$contaminant)

# look at the incidences of taxa observations in negative controls and positive samples.
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("prevalence (PCR and Extraction Controls)") +
  ylab("prevalence (eDNA samples and positive controls)") +
  theme_bw()

ggsave("210121_990_r_get_eDNA_phyloseq__psob-cross-contaminants-prevalence-overview.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 50, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# Removal of contamination with package `decontam`
# ------------------------------------------------

# contamination filtering
psob_decontam <- prune_taxa(!contamdf.freq$contaminant, psob_input)


# Summary after contaminant filtering
(psob_molten <- get_tidy_molten_ps(psob_decontam))
get_default_coverage_plot (psob_molten, facet_var = "SAMPLE.TYPE", taxlev = "PHYLUM", ptitl = "Phyla across all sample types after contaminant removal", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")

ggsave("210121_990_r_get_eDNA_phyloseq__psob-filtered-sample-type-overview.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_default_coverage_plot (psob_molten, taxlev = "PHYLUM", ptitl = "Phyla across all locations after contaminant removal", pxlab = "phyla at all locations", pylab =  "read counts at each location (y scales fixed)")

ggsave("210121_990_r_get_eDNA_phyloseq__psob-filtered-location-overview.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files",
         scale = 3, width = 75, height = 60, units = c("mm"),
         dpi = 500, limitsize = TRUE)

get_molten_ps_description(psob_molten, rank_level = "SUPERKINGDOM", rank_name = "Eukaryota")

# remove temporary object
rm(psob_molten)


# IV. Retain eDNA samples with fish only - and summarise 
# ======================================================

# melt cleaned phyloseq object
(psob_molten <- get_tidy_molten_ps(psob_decontam))


# get info on filtering criteria
unique(psob_molten$PHYLUM)

# isolate and summarize interesting phyla
psob_molten_sponges <- psob_molten %>% filter(PHYLUM == "Porifera")
psob_molten_chordat <- psob_molten %>% filter(PHYLUM == "Chordata")

get_molten_ps_description(psob_molten_sponges, rank_level = "PHYLUM", rank_name = "Porifera")
get_molten_ps_description(psob_molten_chordat, rank_level = "PHYLUM", rank_name = "Chordata", prnt_n = Inf) 

# notable mentions: Mamamlia
# Homo sapiens - human
# Trichosurus vulpecula  - bushtail possum  
# Arctocephalus forsteri - New Zealand fur seal
# Cervus elaphus - Red Deer
# Bos taurus - cattle
# Ovis aries - sheep
# Canis lupus - dog
# Bos frontalis - gayal
# Tursiops truncatus - bottlenose dolphin
# Mus musculus - mouse
# Sus scrofa - pig
# Rattus exulans - Polynesian rat

# notable mentions: Aves
# Gallus gallus - chicken

# keep and check classes Chondrichthyes Actinopteri across entire data
psob_molten_fish <- psob_molten_chordat %>% filter(CLASS %in% c("Chondrichthyes", "Actinopteri"))

get_default_coverage_plot(psob_molten_fish, taxlev = "FAMILY", taxlev_fill = "CLASS", ptitl = "Fish families across all locations after contaminant removal", pxlab = "families at all locations", pylab =  "read counts at each location (y scales fixed)")
get_default_ocurrence_plot(psob_molten_fish, taxlev = "FAMILY", taxlev_fill = "CLASS", ptitl = "Fish families across all locations after contaminant removal", pxlab = "families at all locations", pylab =  "read counts at each location (y scales fixed)")

get_molten_ps_description(psob_molten_fish, rank_level = "PHYLUM", rank_name = "Chordata", prnt_n = Inf) 

# The current data set contains 18588722 sequences across 157 samples and 123 ASV's (123 Chordata, as well as 0 non-Chordata, i.e. ).
#  Sample mean (min., med., max.) coverage is 118400 reads (17, 51474, 831439), and ASV mean (min, median, max) coverage is 151128 reads
#  (5, 3580, 4173587).
#  A tibble: 1 x 11
#    ASV ABUNDANCE SAMPLE LOC.NAME SUPERKINGDOM PHYLUM CLASS ORDER FAMILY GENUS SPECIES
#   <int>     <int>  <int>    <int>        <int>  <int> <int> <int>  <int> <int>   <int>
#    123      1661    157       12            1      1     2    29     48    58      65

# remove ASVs of  controls from eDNA data

# A tibble: 124 x 41
psob_molten_fish_controls <- psob_molten_fish %>% filter(SAMPLE.TYPE %in% c("pcntrol-blnd", "pcntrol-zebra", "ncntrl-pcr", "ncntrol-xtr", "blank")) %>% filter(ABUNDANCE != 0)
asvs_molten_fish_controls <- psob_molten_fish %>% filter(SAMPLE.TYPE %in% c("pcntrol-blnd", "pcntrol-zebra", "ncntrl-pcr", "ncntrol-xtr", "blank")) %>% filter(ABUNDANCE != 0) %>% pull(ASV) %>% unique(.)
length(asvs_molten_fish_controls)

get_default_coverage_plot(psob_molten_fish_controls, facet_var = "SAMPLE.TYPE", taxlev = "FAMILY", taxlev_fill = "CLASS", ptitl = "Fish families across controls after contaminant removal", pxlab = "families at all locations", pylab =  "read counts at each location (y scales fixed)")
get_molten_ps_description(psob_molten_fish_controls, rank_level = "PHYLUM", rank_name = "Chordata", prnt_n = Inf) 

# The current data set contains 2830326 sequences across 32 samples and 84 ASV's (84 Chordata, as well as 0 non-Chordata, i.e. ). 
# Sample mean (min., med., max.) coverage is 88448 reads (17, 5668, 568752), and ASV mean (min, median, max) coverage is 336944 
# reads (40, 7725, 10203800).
# 
# In the following summary variable names shown are hard-coded in function "get_molten_ps_stats": 
# # A tibble: 1 x 11
#     ASV ABUNDANCE SAMPLE LOC.NAME SUPERKINGDOM PHYLUM CLASS ORDER FAMILY GENUS SPECIES
#   <int>     <int>  <int>    <int>        <int>  <int> <int> <int>  <int> <int>   <int>
# 1    84       331     32        1            1      1     2    26     40    48      53


# A tibble: 1,782 x 41
psob_molten_fish_eDNA_samples <- psob_molten_fish %>% filter(SAMPLE.TYPE %in% c("eDNA")) %>% filter(ABUNDANCE != 0)

get_default_coverage_plot(psob_molten_fish_eDNA_samples, taxlev = "FAMILY", taxlev_fill = "CLASS", ptitl = "Fish families in eDNA samples before control removal", pxlab = "families at all locations", pylab =  "read counts at each location (y scales fixed)")
get_molten_ps_description(psob_molten_fish_eDNA_samples, rank_level = "PHYLUM", rank_name = "Chordata", prnt_n = Inf) 

# A tibble: 550 x 41
psob_molten_fish_eDNA_samples <- psob_molten_fish %>% filter(SAMPLE.TYPE %in% c("eDNA")) %>% filter(!ASV %in% asvs_molten_fish_controls)

get_default_coverage_plot(psob_molten_fish_eDNA_samples, taxlev = "FAMILY", taxlev_fill = "CLASS", ptitl = "Fish families in eDNA samples after control removal", pxlab = "families at all locations", pylab =  "read counts at each location (y scales fixed)")
get_molten_ps_description(psob_molten_fish_eDNA_samples, rank_level = "PHYLUM", rank_name = "Chordata", prnt_n = Inf) 


# V. Save eDNA object for further analysis
# ==========================================

# save or load molten state
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210108_990_r_get_eDNA_phyloseq__export_workspace.Rdata")
save(psob_molten_fish_eDNA_samples, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210108_990_r_get_eDNA_phyloseq__exported-eDNA-psob.Rdata")
save(psob_molten_fish_eDNA_samples, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/zenodo/210108_990_r_get_eDNA_phyloseq__exported-eDNA-psob.Rds")

