# **********************************************
# * Create, filter, and write Physloseq object *
# **********************************************
# 12-Dec-2020, 08-Jan-2020

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

# Check sample metadata
# ---------------------

# check REDAME.md 08-Jan-2021 and commit `64f61205e54334e9282616b0fa221c5362aa338c`
# Spelling and order should be should
#   be `sampleid,key,ng-ul,lat-dd,long-dd,pool-content,primer-direction,primer-label,adapter-fwd,flupad,index-fwd,template-primer-fwd,complete-primer-fwd,loc-name,sample-time,sample-type,inside-reserve,sample-name,vol-l,depth-m,notes,xtr-date,row,plate,col,primer-direction_rev,primer-label-rev,adapter-rev,index-rev,template-primer-rev,complete-primer-rev`
#   as in `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/850_prep_q2_predictor-tab__metadata.tsv`
#   as written by `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh`


names(sample_data(psob)) # names are spelled correctly but ordered alphabetically
head(sample_data(psob))  # names are changed from "-" to "." and ordered alphabetically in head command - PhyloSeq sloppy coding

# Save or load initial state after import
# ---------------------------------------
save.image(file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210108_990_r_get_eDNA_phyloseq__import_workspace.Rdata")
save(psob, file = "/Users/paul/Documents/OU_eDNA/201028_Robjects/210108_990_r_get_eDNA_phyloseq__imported-psob.Rdata")
save(psob, file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/zenodo/210108_990_r_get_eDNA_phyloseq__imported-psob.Rds")
