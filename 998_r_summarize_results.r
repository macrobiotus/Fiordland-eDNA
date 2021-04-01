#  **********************************
#  *                                *        
#  *  Get results from long tables  *
#  *                                *
#  **********************************

# I. Load packages
# ================

library("tidyverse") # because we can't stop using it anymore
library("ggrepel")   # to improve plot labels

library("future.apply") # faster handling of large tables
library("data.table")   # faster handling of large tables

library("irr")          # to calculate ICCs between BRIUV and eDNA - see refs below

library("eulerr")       # to compare BRIUV and eDNA
library("ggplotify")

library("vegan")        # for NMDS 
library("indicspecies") # indicator species  - see citation below

library("FactoMineR") # MCA
library("explor")     # check MCA results in browser
library("factoextra") # get MCA results summaries

library("ggpubr") # combine plots -  http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/


# library("nVennR")
# library("UpSetR")    # Conway, J. R., Lex, A. & Gehlenborg, N. 2017 UpSetR: an R package for the
#                      # visualization of intersecting sets and their properties. Bioinformatics 33,
#                      # 2938Ð2940. (doi:10.1093/bioinformatics/btx364)
#                      # 
#                      # documentation at https://rdrr.io/cran/UpSetR/man/upset.html - hard to follow

# II. Read in data
# ================

rm(list = ls(all.names = TRUE))
gc()

long_table <- readRDS(file = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210301_997_r_format_longtables__analysis_input.Rds")
mdl_specs <- read_csv("/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210309_mdl_tablebyspecies.csv")
mdl_genus <- read_csv("/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210309_mdl_tablebygenus.csv")

# III. Format data 
# ================

# properly set factor variables in MdLs data in MdLs table
cols <- c("a.in", "a.out", "b.in", "b.out", "c.in", "c.out")
mdl_specs[cols] <- lapply(mdl_specs[cols], factor) 
mdl_genus[cols] <- lapply(mdl_genus[cols], factor) 

# define BRUV.PRES and set EDNA.PRES and BOTH.PRES
long_table <- long_table %>% mutate( BRUV.PRES = case_when(SAMPLE.TYPE == "BRUV" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table <- long_table %>% mutate( EDNA.PRES = case_when(SAMPLE.TYPE == "eDNA" & ABUNDANCE >= 1 ~ 1, TRUE ~ 0))
long_table <- long_table %>% mutate( BOTH.PRES = case_when(BRUV.PRES == 1 | EDNA.PRES == 1 ~ 1, TRUE ~ 0))


# 16-Mar-2021 add asterisks ("*") to non-NZ species, and ("**") to non-fish (mammals and crustaceans)
#  after checking with list 
#  Roberts, C., Stewart, A., Struthers, C., Barker, J. & Kortet, S. 2019 Checklist of the Fishes of New Zealand. 

nonnz_fish <- c("Asterropteryx", "Banjos", "Benitochromis", "Bostrychus", "Bovichtus", "Caprodon", "Coptodon", "Engraulis", "Gobiesox", "Gymnoscopelus", "Helcogramma", "Microcanthus", "Opistognathus", "Phoxinus", "Sander", "Scobinichthys")
nonnz_othr <- c("Macroctopus", "Jasus", "Arctocephalus", "Balaenoptera", "Tursiops")

long_table <- long_table %>% mutate(GENUS = case_when(GENUS %in% nonnz_fish ~ paste0(GENUS, "*"),
                                        GENUS %in% nonnz_othr ~ paste0(GENUS, "**"),
                                        TRUE                  ~ GENUS)
                                        )

# describe ASV  yield per primer
smpl_eff <- long_table %>% ungroup() %>% select(SET.ID, REP.ID, SAMPLE.TYPE, PRIMER.LABEL, RESERVE.GROUP.LOCATION) %>% filter(SAMPLE.TYPE == "eDNA") %>% arrange(SET.ID, RESERVE.GROUP.LOCATION) %>% print(n = Inf)
smpl_eff_grp <- smpl_eff %>% mutate(PRIMER.LABEL  = gsub(".*Mi", "", PRIMER.LABEL)) %>% group_by(RESERVE.GROUP.LOCATION, PRIMER.LABEL) %>% summarise(n = n()) %>% ungroup()
smpl_eff_grp %>% arrange(PRIMER.LABEL, RESERVE.GROUP.LOCATION)


# IV. Analyse BRUV vs EDNA 
# =========================

# for now using ICC - possibly use something else
# -----------------------------------------------

# show ASV level observations for eDNA and BRUV 
AsvPresByMethd <- long_table %>% select (RESERVE.GROUP, RESERVE.GROUP.LOCATION, SET.ID, BOTH.PRES, BRUV.PRES, EDNA.PRES, SUPERKINGDOM,  PHYLUM,  CLASS,  ORDER,  FAMILY,  GENUS, SPECIES) 

# show species level observations for eDNA and BRUV - loosing ASV level observation
SpcPresByMethd <- long_table %>% select (RESERVE.GROUP, RESERVE.GROUP.LOCATION, SET.ID, BOTH.PRES, BRUV.PRES, EDNA.PRES, SUPERKINGDOM,  PHYLUM,  CLASS,  ORDER,  FAMILY,  GENUS, SPECIES) %>% distinct() %>% print(n = Inf)

# isolate Tibbles for ICC calculation
SpcPresByMethdPhl <- SpcPresByMethd  %>% group_by(PHYLUM)  %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdCls <- SpcPresByMethd  %>% group_by(CLASS)   %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdOrd <- SpcPresByMethd  %>% group_by(ORDER)   %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdFam <- SpcPresByMethd  %>% group_by(FAMILY)  %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES))
SpcPresByMethdGen <- SpcPresByMethd  %>% group_by(GENUS)   %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdSpc <- SpcPresByMethd  %>% group_by(SPECIES) %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 


# calculate ICCS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4913118/
IccPhl <- icc(SpcPresByMethdPhl[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccCls <- icc(SpcPresByMethdCls[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccOrd <- icc(SpcPresByMethdOrd[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccFam <- icc(SpcPresByMethdFam[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccGen <- icc(SpcPresByMethdGen[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccSpc <- icc(SpcPresByMethdSpc[ , 2:3], model = "twoway", type = "agreement", unit = "single")

# build Tibble from ICC results for plotting
icc_temp <- do.call(rbind, Map(data.frame, PHYLUM=IccPhl, CLASS=IccCls, ORDER=IccOrd, FAMILY=IccFam, GENUS=IccGen, SPECIES=IccSpc))
icc_temp_t <- data.table::transpose(icc_temp)
colnames(icc_temp_t ) <- rownames(icc_temp)
rownames(icc_temp_t ) <- colnames(icc_temp)
colnames(icc_temp_t) <- toupper(colnames(icc_temp_t))
iccs <- icc_temp_t %>% rownames_to_column("LEVEL") %>% as_tibble()
iccs <- iccs %>% mutate(LEVEL=factor(LEVEL, levels=LEVEL)) %>% mutate(LEVEL = factor(LEVEL, levels=c("PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS", "SPECIES")))
iccs <- iccs %>% mutate(across(.cols = c("SUBJECTS", "RATERS", "VALUE", "R0", "FVALUE", "DF1", "DF2", "P.VALUE", 
                        "P.VALUE", "CONF.LEVEL", "LBOUND", "UBOUND") , .fns = as.numeric ))

# plot results
#  rating lines https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4913118/
#   greater than 0.90      excellent 
#   between 0.75 and 0.9,  good
#   between 0.5 and 0.75,  moderate
#   values less than 0.5,  poor
p_eb <- ggplot(iccs, aes(x = LEVEL, y = VALUE, ymin = LBOUND, ymax = UBOUND)) +
  geom_point() + geom_errorbar(width = 0.1) + 
  geom_hline(yintercept = 0.9, color="gray", linetype="dashed") +
  geom_hline(yintercept = 0.75, color="gray", linetype="dashed") +
  geom_hline(yintercept = 0.5, color="gray", linetype="dashed") +
  geom_hline(yintercept = 0, color="gray", linetype="dashed") +
  annotate(geom = "text", x = 0.5, y = 1.1, label = "excellent", hjust = 0, color="gray") +
  annotate(geom = "text", x = 0.5, y = 0.85, label = "good", hjust = 0, color="gray") +
  annotate(geom = "text", x = 0.5, y = 0.6, label = "moderate", hjust = 0, color="gray") +
  annotate(geom = "text", x = 0.5, y = 0.4, label = "poor", hjust = 0, color="gray") +
  annotate(geom = "text", x = 0.5, y = -0.1, label = "contradictory", hjust = 0, color="gray") +
  xlab("Taxonomic level") + 
  ylab("ICC (BRUV, eDNA)") +
  scale_x_discrete(labels=paste0(iccs$LEVEL," (n=", iccs$SUBJECTS, ")" )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("210312_998_r_summarize_results_edna_bruv_comp.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 75, height = 75, units = c("mm"),
         dpi = 500, limitsize = TRUE)  

# ICC alternative (still unfinished) 
# ----------------------------------

SpcPresByMethd <- long_table %>% select (RESERVE.GROUP, RESERVE.GROUP.LOCATION, SET.ID, BOTH.PRES, BRUV.PRES, EDNA.PRES, SUPERKINGDOM,  PHYLUM,  CLASS,  ORDER,  FAMILY,  GENUS, SPECIES) %>% distinct() %>% print(n = Inf)

SpcPresByMethdPhl <- SpcPresByMethd  %>% group_by(PHYLUM)  %>% summarise(eDNA = as.logical(sum(EDNA.PRES)), BRUV = as.logical(sum(BRUV.PRES)))
SpcPresByMethdCls <- SpcPresByMethd  %>% group_by(CLASS)   %>% summarise(eDNA = as.logical(sum(EDNA.PRES)), BRUV = as.logical(sum(BRUV.PRES))) 
SpcPresByMethdOrd <- SpcPresByMethd  %>% group_by(ORDER)   %>% summarise(eDNA = as.logical(sum(EDNA.PRES)), BRUV = as.logical(sum(BRUV.PRES))) 
SpcPresByMethdFam <- SpcPresByMethd  %>% group_by(FAMILY)  %>% summarise(eDNA = as.logical(sum(EDNA.PRES)), BRUV = as.logical(sum(BRUV.PRES)))
SpcPresByMethdGen <- SpcPresByMethd  %>% group_by(GENUS)   %>% summarise(eDNA = as.logical(sum(EDNA.PRES)), BRUV = as.logical(sum(BRUV.PRES))) 
SpcPresByMethdSpc <- SpcPresByMethd  %>% group_by(SPECIES) %>% summarise(eDNA = as.logical(sum(EDNA.PRES)), BRUV = as.logical(sum(BRUV.PRES))) 

fitPhl <- euler(SpcPresByMethdPhl[ , 2:3])
fitCls <- euler(SpcPresByMethdCls[ , 2:3])
fitOrd <- euler(SpcPresByMethdOrd[ , 2:3])
fitFam <- euler(SpcPresByMethdFam[ , 2:3])
fitGen <- euler(SpcPresByMethdGen[ , 2:3])
fitSpc <- euler(SpcPresByMethdSpc[ , 2:3])


pPhl <- as.ggplot(plot(fitPhl, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels = list(font=1, cex=0.8))) + labs(subtitle = "Phylum")
pCls <- as.ggplot(plot(fitCls, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels= FALSE)) + labs(subtitle = "Class")
pOrd <- as.ggplot(plot(fitOrd, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels= FALSE)) + labs(subtitle = "Order")
pFam <- as.ggplot(plot(fitFam, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels= FALSE)) + labs(subtitle = "Family")
pGen <- as.ggplot(plot(fitGen, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels= FALSE)) + labs(subtitle = "Genus")
pSpc <- as.ggplot(plot(fitSpc, quantities = list(type = c("counts", "percent"), font=3, round=2, cex=0.8), labels= FALSE)) + labs(subtitle = "Species")

pEuler <- ggarrange(pPhl, pCls, pOrd, pFam, pGen, pSpc,  ncol = 6, nrow = 1)


# V. Show RESERVE.GROUP.LOCATION similarity based on GENUS overlap 
# ===================================================================
 
# check data
print(long_table, n = Inf)
names(long_table)

# copy data to data table
long_table_dt <- data.table(long_table)
setkey(long_table_dt,ASV) 

# aggregate discrete observation of wither method ("BOTH.PRES") per sampling area (RESERVE.GROUP.LOCATION) on GENUS level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
long_table_dt_agg_gen <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]

# reshape to observation matrix digestible by Vegan, discrete observations will be summed per genus
long_table_dt_agg_gen_mat <- as.matrix(data.table::dcast(setDT(long_table_dt_agg_gen), RESERVE.GROUP.LOCATION~GENUS, value.var="BOTH.PRES", sum, fill=0), rownames=TRUE)

jacc_matrix <- vegdist(long_table_dt_agg_gen_mat, distance="jaccard" )
summary(jacc_matrix)

# get a Jaccard distance matrix (distance define by overlap between sites)
#   see https://rpubs.com/CPEL/NMDS
#   see https://peat-clark.github.io/BIO381/veganTutorial.html
#   see https://fromthebottomoftheheap.net/2013/01/12/decluttering-ordination-plots-in-vegan-part-1-ordilabel/
#   see https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo


# run metaMDS - at pressence does just use the presence absence matrix (to keep genus scrores), but distance matrix is possible as well
long_table_dt_agg_gen_mat_jacc_NMS <-  metaMDS(long_table_dt_agg_gen_mat, distance="jaccard",  k = 2, maxit = 5000,  trymax = 5000, wascores = TRUE, shrink = TRUE)

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
   theme_bw()
ggsave("210312_998_r_summarize_results_jaccard.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 75, height = 75, units = c("mm"),
         dpi = 500, limitsize = TRUE)     


# VI. ANOSIM to test wether or not genus composition based on factors are significantly different
# ==========================================================================================

# On Anosim: 
# 1. CLARKE, K. R. 1993 Non-parametric multivariate analyses of changes in
# community structure. Austral Ecol. 18, 117Ð143.
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


# VII. Show indicator genera for inside/outside each/all reserve(s)
# ============================================================================
# indicator species analysis
# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

# Using package indicspecies
# De C‡ceres, M., Legendre, P. & Moretti, M. 2010 Improving indicator species
# analysis by combining groups of sites. Oikos 119, 1674Ð1684.
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
# Signif. codes:  0 Ô***Õ 0.001 Ô**Õ 0.01 Ô*Õ 0.05 Ô.Õ 0.1 Ô Õ 1


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


# VIII. Multiple Correspondence analysis
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
    ylab(paste0("Dim. 2 (", dim_2_perc, "% Variance)"))
    
ggsave("210312_998_r_summarize_results_mca.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1.5, width = 150, height = 150, units = c("mm"),
         dpi = 500, limitsize = TRUE)     

# see https://rpkgs.datanovia.com/factoextra/reference/fviz_contrib.html
p_cntrb <- fviz_contrib(mca1, choice="var", axes = 1, top = 10, fill = "lightgray", color = "grey") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(title = element_blank())
ggsave("210312_998_r_summarize_results_mca_dim1.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1.3, width = 50, height = 25, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# IX. Create a better map for manuscript  
# =====================================

# ** started 29-03-2021 - not finished **

# following 
#  https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
#  https://semba-blog.netlify.app/10/20/2018/genetic-connectivity-in-western-indian-ocean-region/

# aggregate discrete observation of wither method ("BOTH.PRES") per sampling area (RESERVE.GROUP.LOCATION) on GENUS level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
# https://gis.stackexchange.com/questions/243569/simplify-polygons-of-sf-object

library("sf")          # simple feature objects
library("rmapshaper")  # simplify shape file layers
library("ggsflabel")   # label simple feature in ggplot  https://github.com/yutannihilation/ggsflabel - possibly inluded in ggplot
library("ggsn")        # scale bar in ggplot map

# library("tidyverse")
# library("rgdal")
# library("rnaturalearth")
# library("rnaturalearthdata")


# 1.) read the shape file
# -----------------------
nzshp_hires = read_sf("/Users/paul/GIS/NZ_coast/NZ_Coast_isl.shp")


# 2.) re-project shape file (to simple WGS84)
# -------------------------------------------
# https://r-spatial.github.io/sf/reference/st_transform.html
# https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/reproject-vector-data/
# transfrom to simple WGS84, EPSG:4326, "+proj=longlat +datum=WGS84 +no_defs"
nzshp_hires_WGS84 <- st_transform(nzshp_hires, crs = 4326)


# 3.) simplify shape for low resolution map insets
# -------------------------------------------------
# https://gis.stackexchange.com/questions/243569/simplify-polygons-of-sf-object
nzshp_lores_WGS84 <- rmapshaper::ms_simplify(input = as(nzshp_hires_WGS84, 'Spatial')) %>% st_as_sf()


# 4.) define a bounding box around the field work area
# -------------------------------------------------
# https://geocompr.github.io/post/2019/ggplot2-inset-maps/
# mins and max of point coordinates, and 0.1 degree added
bb_fwork = st_as_sfc(st_bbox(c(xmin = (166.5-0.1), xmax = (167.0+0.1), ymax = (-46.04-0.1), ymin = (-45.52+0.1)), crs = st_crs(4326)))


# 5.) draw overview map including  bounding box
# ----------------------------------------
map_inst <- ggplot(data = nzshp_lores_WGS84) +
    geom_sf(fill = "grey93", color = "red", lwd = 0.5) +
    geom_sf(data = bb_fwork, fill = NA, color = "darkred", size = 1) +
    theme_void()
    

# 6.) draw location map including and add data
# --------------------------------------------

# get data to be shown
# ````````````````````
long_table_dt_map <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("MH.GPS.LAT", "MH.PPS.LONG", "RESERVE.GROUP", "RESERVE.GROUP.INSIDE", "RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]


# group data for plotting and manuscript detail expansion
# ```````````````````````````````````````````````````````
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


# get bounding boxes
# ```````````````````

# https://stackoverflow.com/questions/54696440/create-polygons-representing-bounding-boxes-for-subgroups-using-sf
# function calculates angle with respect to polygon centroid.
# we need this to order the polygon correctly
calc_angle <- function(lon,lat) {
  cent_lon <- mean(lon)
  cent_lat <- mean(lat)
  ang <- atan2(lat - cent_lat, lon - cent_lon)

  return(ang)
}

bbox <- long_table_dt_map %>%
  group_by(RESERVE.GROUP.LOCATION) %>%
  summarise(xmin = min(MH.PPS.LONG) -0.01 ,ymin = min(MH.GPS.LAT) -0.01, xmax=max(MH.PPS.LONG) +0.01,  ymax = max(MH.GPS.LAT) +0.01) %>%
  gather(x,lon,c('xmin','xmax')) %>%
  gather(y,lat,c('ymin','ymax')) %>%
  st_as_sf(coords=c('lon','lat'),crs=4326,remove=F) %>%
  group_by(RESERVE.GROUP.LOCATION) %>%
  mutate(angle = calc_angle(lon,lat)) %>%
  arrange(angle) %>%
  summarise(do_union=FALSE) %>%
  st_cast('POLYGON')


# draw main map
# ``````````````
map_main <- ggplot(data = nzshp_lores_WGS84) +
    geom_sf(fill = "lightgrey") +
    geom_sf(data=bbox, fill = NA, color = "red", size = 1) + 
    coord_sf( xlim = c((166.5-0.1), (167.0+0.1)), ylim = c((-46.04-0.1),(-45.52+0.1)), expand = FALSE) +
    geom_point(data = long_table_dt_map, aes(x = MH.PPS.LONG, y = MH.GPS.LAT, shape = RESERVE.GROUP), color = "darkred", size = 4) +
    geom_point(data = long_table_dt_map, aes(x = MH.PPS.LONG, y = MH.GPS.LAT, shape = RESERVE.GROUP), color = "red", size = 3) +
    geom_sf_label(data=bbox, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 0.06, nudge_y = 0.06) + 
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.position=c(.9,.1), 
          legend.background = element_blank(), 
          legend.key=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    annotation_custom(ggplotGrob(map_inst), xmin = 166.35, xmax = 166.7, ymin = -45.62, ymax = -45.45)

ggsave("210401_998_r_summarize_results_fig1_draft.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 125, height = 175, units = c("mm"),
         dpi = 500, limitsize = TRUE)


# draw genus map (not very useful)
# ````````````````````````````````

map_genera <- ggplot(data = nzshp_lores_WGS84) +
    geom_sf(fill = "lightgrey") +
    geom_sf(data=bbox, fill = NA, color = "red", size = 1) + 
    coord_sf( xlim = c((166.5-0.1), (167.0+0.1)), ylim = c((-46.04-0.1),(-45.52+0.1)), expand = FALSE) +
    geom_point(data = long_table_dt_map, aes(x = MH.PPS.LONG, y = MH.GPS.LAT, shape = RESERVE.GROUP), color = "darkred", size = 4) +
    geom_point(data = long_table_dt_map, aes(x = MH.PPS.LONG, y = MH.GPS.LAT, shape = RESERVE.GROUP), color = "red", size = 3) +
    geom_sf_label(data=bbox, aes(label = RESERVE.GROUP.LOCATION), nudge_x = 0.06, nudge_y = 0.06) + 
    geom_label_repel(data = long_table_dt_map, aes(x = MH.PPS.LONG, y = MH.GPS.LAT, label = GENUS), max.overlaps = Inf, point.size = NA) +
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.position=c(.9,.1), 
          legend.background = element_blank(), 
          legend.key=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

ggsave("210401_998_r_summarize_results_map_genara_draft.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 300, height = 500, units = c("mm"),
         dpi = 500, limitsize = TRUE)





# X. Combine plots for manuscript 
# ================================

# arrange plots - version 1 - with ICCs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggarrange(ggarrange(p_eb, p_cntrb, p_nmds, ncol = 3, labels = c("(a)", "(c)", "(d)")),
          ggarrange(p_mca, ncol = 1, labels = c("(b)")),
          nrow = 2, heights = c(3, 7))
ggsave("210312_998_r_summarize_results_fig2_draft_ICC.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 250, height = 300, units = c("mm"),
         dpi = 500, limitsize = TRUE)

# arrange plots - version with Venn diagrams - version 1
ggarrange(ggarrange(ggarrange(pPhl, pCls, pOrd, pFam, pGen, pSpc, ncol = 3, nrow = 2), p_nmds, p_cntrb, ncol = 3, nrow = 1, labels = c("(a)", "(c)", "(d)")),
          ggarrange(p_mca, ncol = 1, labels = c("(c)")),
          nrow = 2, heights = c(3, 7))

# arrange plots - version 2 -  with Venn diagrams
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# unifinished

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


# arrange plots - version 3 -  with maps and Venndiagrams
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggarrange(
  ggarrange(map_main, ncol = 1, nrow = 1, labels = c("(a)")),
  ggarrange(map_inst, p_nmds, ncol = 1, nrow = 2, labels = c("", "(b)"), heights = c(3, 7)),
  ggarrange(pPhl, pCls, pOrd, pFam, pGen, pSpc, ncol = 1, nrow = 6, labels = c("(c)")),
  ncol = 3, nrow = 1)
