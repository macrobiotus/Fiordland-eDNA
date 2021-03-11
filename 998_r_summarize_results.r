#   **********************************
#   *                                *        
#   *  Get results from long tables  *
#   *                                *
#   **********************************
#   10-Mar-2022


# I. Load packages
# ================

library("tidyverse") # because we can't stop using it anymore
library("data.table")   # faster handling of large tables
library("future.apply") # faster handling of large tables
library("vegan")
library("nVennR")
library("UpSetR")    # Conway, J. R., Lex, A. & Gehlenborg, N. 2017 UpSetR: an R package for the
                     # visualization of intersecting sets and their properties. Bioinformatics 33,
                     # 2938–2940. (doi:10.1093/bioinformatics/btx364)
                     # 
                     # documentation at https://rdrr.io/cran/UpSetR/man/upset.html - hard to follow


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


# IV. Analyse observations using MDS 
# =================================

# check data
print(long_table, n = Inf)
names(long_table)

# copy data to data table
long_table_dt <- data.table(long_table)
setkey(long_table_dt,ASV) 

# aggregate ABUNDANCE per sampling area (RESERVE.GROUP.LOCATION) down to GENUS level  
#   https://stackoverflow.com/questions/16513827/summarizing-multiple-columns-with-data-table
long_table_dt_agg <- long_table_dt[, lapply(.SD, sum, na.rm=TRUE), by=c("RESERVE.GROUP.LOCATION", "SUPERKINGDOM",  "PHYLUM",  "CLASS",  "ORDER",  "FAMILY",  "GENUS"), .SDcols=c("ABUNDANCE") ]

# define a column with presence absence (ASVPRESENT) if any ABUNDANCE larger then 0
long_table_dt_agg <- long_table_dt_agg[ , ASVPRESENT :=  fifelse(ABUNDANCE == 0 , 0, 1, na=NA)]
print(long_table_dt_agg, n = Inf)

# define a presence absence matrix of genera by area
area_genera_pres <- as.matrix(data.table::dcast(setDT(long_table_dt_agg), RESERVE.GROUP.LOCATION~GENUS, value.var="ASVPRESENT", fill=0), rownames=TRUE)

# get a Jaccard distance matrix (distance define by overlap between sites)
#   see https://rpubs.com/CPEL/NMDS
#   see https://peat-clark.github.io/BIO381/veganTutorial.html
#   see https://fromthebottomoftheheap.net/2013/01/12/decluttering-ordination-plots-in-vegan-part-1-ordilabel/
#   see https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo


area_genera_pres_jacc <- vegan::vegdist(area_genera_pres, method="jaccard", binary=FALSE, diag=TRUE, upper=TRUE, na.rm = FALSE)
class(area_genera_pres_jacc) # for dist objects metaMDS skips transformation and distance calculation

# run metaMDS - at pressence does just use the presence absence matrix (to keep genus scrores), but distance matrix is possible as well
area_genera_pres_jacc_NMS <-  metaMDS(area_genera_pres, distance="jaccard",  k = 4, maxit = 1000,  trymax = 1000, wascores = TRUE, autotransform = FALSE, shrink = TRUE)

# basic plots - variant A
ordiplot(area_genera_pres_jacc_NMS, type = "none") 
orditorp(area_genera_pres_jacc_NMS, display =  "sites", cex = 1.25,air = 0.01)
orditorp(area_genera_pres_jacc_NMS, display = "species", col= "red", air=0.01)

# basic plots - variant B
ordipointlabel(area_genera_pres_jacc_NMS)


# basic plots - variant C
plot(area_genera_pres_jacc_NMS, type = "n", scaling = 3)
ordilabel(area_genera_pres_jacc_NMS, display = "species", font = 2)
ordilabel(area_genera_pres_jacc_NMS, display = "sites", font = 3, fill = "hotpink", col = "blue")

# basic plots - variant D
area_genera_pres_jacc_NMS_plot <- ordipointlabel(area_genera_pres_jacc_NMS, cex = 2)
orditkplot(area_genera_pres_jacc_NMS_plot)
# result saved at /Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/210311_998_r_summarize_results_area_genera_pres_jacc_NMS_plot.pdf


# V. Get other results information 
# =============================

# Construction site: testing venn diagrams
# -----------------------------------------

# four dimension venn plot
library("ggVennDiagram")
genes <- paste("gene",1:1000,sep="")
set.seed(20190708)
x <- list(A=sample(genes,300),B=sample(genes,525),C=sample(genes,440),D=sample(genes,350))
ggVennDiagram(x)
