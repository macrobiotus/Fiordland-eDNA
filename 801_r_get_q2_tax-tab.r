# Import Blast results, format, get new taxonomy strings, and write out compatible with Qiime 2 
# ==============================================================================================

# last run 7-Dec-2020, then 2-Feb-2021, then after 18-Nov-2022 
# see https://ropensci.org/tutorials/taxize_tutorial/
#  for handling blast data and getting correct taxonomy strings from the net

# load packages
# --------------
rm(list = ls(all.names = TRUE))
gc()

library("Biostrings") # read fasta file 
library("magrittr")   # more pipes
library("furrr")      # parallel purrrs - for loading 
library("blastxml")   # read blast xml - get via `library(devtools); install_github("BigelowLab/blastxml")`
library("tidyverse")  # work using tibbles
library("janitor")    # clean column names
library("taxonomizr") # query taxon names
library("ggpubr")     # stat_regline_equation
# library("purrr")      # dplyr applies

# Part I: Load Blast results
# --------------------------

# define file path components for listing 
blast_results_folder <- "/Users/paul/Documents/OU_eDNA/201126_preprocessing/blast"
blast_results_pattern <- "*.xml$"

# read all file into lists for `lapply()` -  new file: 68M Aug  4 13:11
blast_results_files <- list.files(path=blast_results_folder, pattern = blast_results_pattern, full.names = TRUE)

# enable multithreading - only useful for multiple files
plan(multicore) 

# takes 7-10 hours on four cores - avoid by reloading full object from disk 
blast_results_list <- furrr::future_map(blast_results_files, blastxml_dump, form = "tibble", .progress = TRUE) 

names(blast_results_list) <- blast_results_files # works

# save(blast_results_list, file="/Users/paul/Documents/OU_eDNA/201028_Robjects/221117_get_q2_tax-tab__blast_results_list.Rdata")
load(file="/Users/paul/Documents/OU_eDNA/201028_Robjects/221117_get_q2_tax-tab__blast_results_list.Rdata", verbose = TRUE)

# code below inserted after 17.11.2022 for EDNA revision

# flatten blast list for further use
flat_blast_results_list <- blast_results_list %>% 
                             bind_rows(, .id = "src" ) %>%  # add source file names as column elements
                             clean_names(.) %>%             # clean columns names 
                             group_by(iteration_query_def)  # group by sequence hash    
                             
# add variable and factor indicating highest bit score
flat_blast_results_list  %<>%  mutate(max_hsp_bit_score = max(hsp_bit_score))
flat_blast_results_list  %<>%  mutate(max_hsp_bit_score_lgl = if_else(hsp_bit_score == max_hsp_bit_score, TRUE, FALSE ))

# for reviewers
# investigate Bit score cut-off as per Liu B, Gibbons T, Ghodsi M, and Pop M. 2010. 
# MetaPhyler: Taxonomic profiling for metagenomic sequences. 2010 IEEE International
# Conference on Bioinformatics and Biomedicine (BIBM).

# add variable marking above average bit scores
# - get the relationship between bitscore and query length
bs_vs_qlen <- lm(hsp_bit_score ~ iteration_query_len, data = flat_blast_results_list)
summary(bs_vs_qlen)

# - average bit score ("avgbits") is 1.7480 x query length
avgbits <- coefficients(bs_vs_qlen)[2]

# - define average Bit for given query length
flat_blast_results_list %<>% mutate(avg_hsp_bit_score = iteration_query_len * avgbits)

# - mark above-average maximum-bitscored hsps as TRUE 
flat_blast_results_list %<>% mutate(abvavg_max_hsp_bit_score_lgl = if_else(max_hsp_bit_score >= avg_hsp_bit_score, TRUE, FALSE ))

# inspect Blast-score related variables
# - to verify previous mutates
flat_blast_results_list[c("iteration_query_def", "hsp_bit_score", "avg_hsp_bit_score",  "max_hsp_bit_score", "max_hsp_bit_score_lgl", "abvavg_max_hsp_bit_score_lgl")] %>% print(n = 50)


# inspect bit scores vs sequence lengths graphically 
#  -  to show sequences how are bit score in comparison to the average
ggplot(flat_blast_results_list, aes(x = iteration_query_len, y= hsp_bit_score)) + 
    geom_hex(aes(fill = factor(max_hsp_bit_score_lgl), colour = after_stat(count))) +
    geom_smooth(method = "lm") + 
    stat_regline_equation(label.x.npc = "center") + 
    theme_bw() +
    ggtitle("sequence query length vs. bit scores of raw blast results") + 
    ylab("bit score") + xlab("query sequence length") +
    scale_fill_discrete(name = "highest bit score\nper query")

ggsave("/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/221118_analysis_outputs/221118_801_bash_bitscore_seqlen_regression_plot.pdf")

# keep only highest-scoring bit scores across iteration_query_def groups (n=1,927)
blast_results <- flat_blast_results_list %>%  slice(which.max(hsp_bit_score))      # save subset

# save object and some time by reloading it - comment in if necessary
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nrow(blast_results) # 1927 now, was 1914 after last blasting, from orginally 2171 after denoising

# Part II: Re-annotate Blast results
# ----------------------------------

# prepareDatabase not needed to be run multiple times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 10-Aug-2020: check external hard drive for readily created database files
# prepareDatabase(sqlFile = "accessionTaxa.sql", tmpDir = "/Users/paul/Sequences/References/taxonomizR/", vocal = TRUE) # takes a very long time - avoid by reloading full object from disk

# function for mutate to convert NCBI accession numbers to taxonomic IDs
# - path updated to external drive 17.04.2020 - defunct, but kept this db at "/Users/paul/Sequences/References/"
#   get_taxid <- function(x) {accessionToTaxa(x, "/Volumes/HGST1TB/Users/paul/Sequences/References/taxonomizR/accessionTaxa.sql", version='base')}
# - path updated to newly downloaded databse (from 15.11.22) on 23.11.2022
get_taxid <- function(x) {accessionToTaxa(x, "/Users/paul/Sequences/References/taxonomizR_221115/taxonomizr.sqlite", version = 'base')}

# function for mutate to use taxonomic IDs and add taxonomy strings
get_strng <- function(x) {getTaxonomy(x, "/Users/paul/Sequences/References/taxonomizR_221115/taxonomizr.sqlite")}

# continue here after 23.11.2022

# add tax ids to table for string lookup - probably takes long time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
blast_results_appended <- blast_results %>% mutate(tax_id = get_taxid(hit_accession)) # takes some time... 
length(blast_results_appended$tax_id) 

# save(blast_results_appended, file="/Users/paul/Documents/OU_eDNA/201028_Robjects/210202_get_q2_tax-tab__blast-noenv_sliced_taxonomy.Rdata")
load(file="/Users/paul/Documents/OU_eDNA/201028_Robjects/210202_get_q2_tax-tab__blast-noenv_sliced_taxonomy.Rdata", verbose=TRUE)

# look up taxonomy table - takes a long time, needs external database.
tax_table <- as_tibble(get_strng(blast_results_appended$tax_id), rownames = "tax_id") %>% mutate(tax_id= as.numeric(tax_id))

nrow(tax_table) 

# getting a tax table without duplicates to enable proper join command later
tax_table <- tax_table %>% arrange(tax_id) %>% distinct(tax_id, superkingdom, phylum, class, order, family, genus, species, .keep_all= TRUE)

# checks
head(tax_table)
nrow(tax_table)
all(!duplicated(tax_table)) #        and no duplicated tax ids anymore
lapply(list(blast_results_appended,tax_table), nrow) # first 10302, second deduplicated and with 2761 - ok 

# https://stackoverflow.com/questions/5706437/whats-the-difference-between-inner-join-left-join-right-join-and-full-join
blast_results_final <- left_join(blast_results_appended, tax_table, copy = TRUE) 
nrow(blast_results_final) # 10302 - table has correct length now 

# Part III: Format Blast results for export
# -----------------------------------------

# correcting factors
blast_results_final %>% ungroup(.) %>% mutate(src = as.factor(src)) -> blast_results_final
levels(blast_results_final$src) 

# diagnostic plot - ok 
ggplot(blast_results_final, aes(x = src, y = phylum, fill = phylum)) + 
    geom_bar(position="stack", stat="identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# adjust this if reading in multiple source files and they need to be re-ordered on the plot
# blast_results_final$src <- plyr::revalue(blast_results_final$src, c("/multiple/path/to/files/"))

# diagnostic plot -ok 
# ggplot(blast_results_final, aes(x = src, y = phylum, fill = phylum)) + 
#     geom_bar(position="stack", stat="identity") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

# omitting factor level correction, didn't work as expected
# blast_results_final$src <- factor(blast_results_final$src, levels = c("18S eukaryotes"))

# save object and some time by reloading it
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save(blast_results_final, file="/Users/paul/Documents/OU_eDNA/201028_Robjects/210202_get_q2_tax-tab__blast-noenv_with-ncbi_taxonomy.Rdata")
load(file="/Users/paul/Documents/OU_eDNA/201028_Robjects/210202_get_q2_tax-tab__blast-noenv_with-ncbi_taxonomy.Rdata")

# Part II: Format taxonomy table for export and export  
# -----------------------------------------------------

# Formatting taxonomy table - selecting relevant columns
#   to match `#OTUID`, `taxonomy`, `confidence`
q2taxtable <- blast_results_final %>% 
  select("iteration_query_def", "superkingdom", "phylum", "class", "order", "family",
    "genus", "species", "hsp_bit_score") %>% mutate_all(replace_na, "undefined")

# - Insert further cleanup code here if desirable - 

# formatting `taxonomy` column
q2taxtable <- q2taxtable %>% unite(taxonomy, c("superkingdom", "phylum", "class", "order",
  "family", "genus", "species"), sep = ";", remove = TRUE, na.rm = FALSE)

names(q2taxtable) <- c("#OTUID", "taxonomy", "confidence")

# Extending taxonomy table from fasta file 
# ----------------------------------------

# read fasta used for Blasting
fnapth <- "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-seq.fasta"
fna = readDNAStringSet(fnapth)
length(names(fna)) # back to 2171

# get length of main table without missing hash values
length(q2taxtable$taxonomy) # 1914 as expected

# add missing hash values from fasta to main table to get complete table  
q2taxtable <- left_join(enframe(names(fna), value = '#OTUID'), q2taxtable, by = c('#OTUID'))

# add indicative taxonomy strings
q2taxtable <- q2taxtable %>% mutate_at(vars(taxonomy), ~replace_na(., "nomatch;nomatch;nomatch;nomatch;nomatch;nomatch;nomatch"))
q2taxtable <- q2taxtable %>% mutate_at(vars(confidence), ~replace_na(., "0"))

q2taxtable$name <- NULL

# Part VI: Export tsv  
# -----------------------------------------------------

# qiime2R compatibility, added 4-May-2020, doesn't help
# names(q2taxtable) <- c("Feature.ID", "Taxon", "Confidence")

write_tsv(q2taxtable, file = "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/800_12S_single_end_ee3-seq_q2taxtable.tsv",
  append = FALSE, col_names = TRUE)
