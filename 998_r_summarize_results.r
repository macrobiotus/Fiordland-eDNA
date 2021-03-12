#  **********************************
#  *                                *        
#  *  Get results from long tables  *
#  *                                *
#  **********************************
#  10-Mar-2022

# I. Load packages
# ================

library("tidyverse") # because we can't stop using it anymore
library("data.table")   # faster handling of large tables
library("future.apply") # faster handling of large tables
library("irr")

library("vegan")
library("ggrepel")


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


# IV. Analyse BRUV vs EDNA 
# =========================

# for now using ICC - possibly use something else
# -----------------------------------------------

# show ASV level observations for eDNA and BRUV 
AsvPresByMethd <- long_table %>% select (RESERVE.GROUP, RESERVE.GROUP.LOCATION, SET.ID, BOTH.PRES, BRUV.PRES, EDNA.PRES, SUPERKINGDOM,  PHYLUM,  CLASS,  ORDER,  FAMILY,  GENUS, SPECIES) 

# show species level observations for eDNA and BRUV - loosing ASV level observation
SpcPresByMethd <- long_table %>% select (RESERVE.GROUP, RESERVE.GROUP.LOCATION, SET.ID, BOTH.PRES, BRUV.PRES, EDNA.PRES, SUPERKINGDOM,  PHYLUM,  CLASS,  ORDER,  FAMILY,  GENUS, SPECIES) %>% distinct() %>% print(n = Inf)

# isolate Tibbles for ICC calculation
SpcPresByMethdPhl <- SpcPresByMethd  %>% group_by(PHYLUM) %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdCls <- SpcPresByMethd  %>% group_by(CLASS) %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdOrd <- SpcPresByMethd  %>% group_by(ORDER) %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdFam <- SpcPresByMethd  %>% group_by(FAMILY) %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES))
SpcPresByMethdGen <- SpcPresByMethd  %>% group_by(GENUS) %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 
SpcPresByMethdSpc <- SpcPresByMethd  %>% group_by(SPECIES) %>% summarise(SPEC.OBSCNT.EDNA = sum(EDNA.PRES), SPEC.OBSCNT.BRUV = sum(BRUV.PRES)) 


# calculate ICCS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4913118/
IccPhl <- icc(SpcPresByMethdPhl[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccCls <- icc(SpcPresByMethdCls[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccOrd <- icc(SpcPresByMethdOrd[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccFam <- icc(SpcPresByMethdFam[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccGen <- icc(SpcPresByMethdGen[ , 2:3], model = "twoway", type = "agreement", unit = "single")
IccSpc <- icc(SpcPresByMethdSpc[ , 2:3], model = "twoway", type = "agreement", unit = "single")

# continue this code if Pearson is desirable
# foo1 <- cor.test(SpcPresByMethdPhl$SPEC.OBSCNT.EDNA,  SpcPresByMethdPhl$SPEC.OBSCNT.BRUV,  method=c("pearson"))
# foo2 <- cor.test(SpcPresByMethdCls$SPEC.OBSCNT.EDNA,  SpcPresByMethdCls$SPEC.OBSCNT.BRUV,  method=c("pearson"))
# foo3 <- cor.test(SpcPresByMethdOrd$SPEC.OBSCNT.EDNA,  SpcPresByMethdOrd$SPEC.OBSCNT.BRUV,  method=c("pearson"))
# foo4 <- cor.test(SpcPresByMethdFam$SPEC.OBSCNT.EDNA,  SpcPresByMethdFam$SPEC.OBSCNT.BRUV,  method=c("pearson"))
# foo5 <- cor.test(SpcPresByMethdGen$SPEC.OBSCNT.EDNA,  SpcPresByMethdGen$SPEC.OBSCNT.BRUV,  method=c("pearson"))
# foo6 <- cor.test(SpcPresByMethdSpc$SPEC.OBSCNT.EDNA,  SpcPresByMethdSpc$SPEC.OBSCNT.BRUV,  method=c("pearson"))

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
ggplot(iccs, aes(x = LEVEL, y = VALUE, ymin = LBOUND, ymax = UBOUND)) +
  geom_point() + geom_errorbar(width = 0.1) + 
  geom_hline(yintercept = 0.9, color="gray", linetype="dashed") +
  geom_hline(yintercept = 0.75, color="gray", linetype="dashed") +
  geom_hline(yintercept = 0.5, color="gray", linetype="dashed") +
  geom_hline(yintercept = 0, color="gray", linetype="dashed") +
  annotate(geom = "text", x = 0.5, y = 1.1, label = "excellent", hjust = 0, color="gray") +
  annotate(geom = "text", x = 0.5, y = 0.8, label = "good", hjust = 0, color="gray") +
  annotate(geom = "text", x = 0.5, y = 0.6, label = "moderate", hjust = 0, color="gray") +
  annotate(geom = "text", x = 0.5, y = 0.4, label = "poor", hjust = 0, color="gray") +
  xlab("Taxonomic level") + 
  ylab("Intra Class Correlation (BRUV, eDNA)") +
  scale_x_discrete(labels=paste0(iccs$LEVEL," (n=", iccs$SUBJECTS, ")" )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("210312_998_r_summarize_results_edna_bruv_comp.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 100, height = 100, units = c("mm"),
         dpi = 500, limitsize = TRUE)  
  
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

ggplot(long_table_dt_agg_gen_mat_jacc_NMS.scores, aes(x = NMDS1, y = NMDS2)) +
   geom_point(size = 6, colour = "red", shape = c(16,16,17,17,18,18)) +
   geom_point(size = 6, colour = "red", shape = c(16,16,17,17,18,18)) +
   geom_label_repel(aes(label=RESERVE.GROUP.LOCATION), point.padding = 0.5) +
   theme_bw()
ggsave("210312_998_r_summarize_results_jaccard.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components",
         scale = 1, width = 100, height = 100, units = c("mm"),
         dpi = 500, limitsize = TRUE)     


# VI. ANOSIM Test show that RESERVE.GROUP.LOCATIONs are significantly different
# ==================================================================================

# A. Show that RESERVE.GROUP.LOCATIONs are significantly different 
# -----------------------------------------------------------------
# https://jkzorz.github.io/2019/06/11/ANOSIM-test.html
# - To test if there is a statistical difference between the fish communities of two or more groups of samples.
# - Null Hypothesis: there is no difference between the microbial communities of your groups of samples.

# aggregate discrete observation of wither method ("BOTH.PRES") per sampling area (RESERVE.GROUP.LOCATION) on GENUS level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
long_table_dt_agg_gen_sets <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("SET.ID", "RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("BOTH.PRES") ]

# rename SET.ID to circumvent naming snafu with package data.table
setnames(long_table_dt_agg_gen_sets, "SET.ID", "SET_ID")

# reshape to observation matrix digestible by Vegan, discrete observations will be summed per genus
long_table_dt_agg_gen_mat_sets <- as.matrix(data.table::dcast(setDT(long_table_dt_agg_gen_sets), SET_ID~GENUS, value.var="BOTH.PRES", sum, fill=0), rownames=TRUE)

# get grouping variable of  RESERVE.GROUP.LOCATION
groupings <- as_tibble(long_table_dt_agg_gen_sets %>% select( SET_ID, RESERVE.GROUP.LOCATION)) %>% distinct( )

long_table_dt_agg_gen_mat_sets_ano <-  anosim(long_table_dt_agg_gen_mat_sets, groupings$RESERVE.GROUP.LOCATION, distance = "jaccard", permutations = 9999)

# Dissimilarity: jaccard 
# 
# ANOSIM statistic R: 0.1992 
#       Significance: 0.0244 
# 
# Permutation: free
# Number of permutations: 9999

# B. Show that [inside/outside MR] are significantly different 
# -----------------------------------------------------------------

# [not done yet]


# VII. Show indicator genera for inside/outside each/all reserve(s)
# ============================================================================
# indicator species analysis
# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

# A. Show that RESERVE.GROUP.LOCATIONs are significantly different 
# -----------------------------------------------------------------

# [not done yet]

# B. Show that [inside/outside MR] are significantly different 
# -----------------------------------------------------------------

# [not done yet]


# VIII. Multiple Correspondence analysis
# =======================================

# Construction site: testing venn diagrams
# -----------------------------------------

# four dimension venn plot
library("ggVennDiagram")
genes <- paste("gene",1:1000,sep="")
set.seed(20190708)
x <- list(A=sample(genes,300),B=sample(genes,525),C=sample(genes,440),D=sample(genes,350))
ggVennDiagram(x)
