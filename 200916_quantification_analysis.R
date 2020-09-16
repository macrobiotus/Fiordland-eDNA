# *****************************************************************************
# * associating fluorescence signal with DNA concentrations, lab book page 48 *
# *****************************************************************************
# 16-Sep-2020

# load packages
# =============
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse")   # work using tibbles and ggplot
library("readxl")      # read excel sheets 
library("openxlsx")    # write Excel tables
library("cowplot")     # multiple plots per page

# read-in and format data 
# =======================

# measured qbit values from lab book for plate positions H1  - H8 
qbit_values <- c(0, 9.1, 4.64, 2.21, 1.15, 0.58, 0.282, 0.148)

# pathnames to Excel exports of plate quantifications
paths <- list(
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Eplate1_20200914 234405.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Eplate2_20200914 234524.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Uplate1_20200914 234210.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Uplate2_20200914 234305.xlsx")

tib_msr <- read_excel(paths[[4]], sheet = 2, col_names = TRUE, skip = 22, .name_repair = "unique") %>%
  rename ("well" = "Well", "pos" = `Well Position`, "cycl" = `Cycle Number`, "picg" =  "SYBR")

# correct and set columns 
tib_msr$cycl <- NULL 
tib_msr$conc <- NA

# add qbit values 
tib_msr$conc[which(tib_msr$pos %in% c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8"))] <- qbit_values

# calculate variables  
# ==================

# see https://stats.stackexchange.com/questions/176595/simple-log-regression-model-in-r/176612

# plot quantification data
plot(conc ~ picg, data = tib_msr[complete.cases(tib_msr), ])

# fit model
fit <- lm(conc ~ picg, data = tib_msr[complete.cases(tib_msr), ] )

# check fit - good fit
summary(fit)

# predict new values and plot 
new_concs = predict(fit, newdata <-  tib_msr[!complete.cases(tib_msr), ], interval="confidence")
plot(new_concs[ ,1] ,  tib_msr[!complete.cases(tib_msr), ]$picg, 
  main = "predicted relationship pico-green to DNA concentration",
  ylab = "fluorescence",
  xlab = "[ng/µl]")
  matlines(new_concs, tib_msr[!complete.cases(tib_msr), ]$picg , lwd=1) 
    
# backtransform concentration values - not used  - formerly contained function call 
tib_imp <- cbind( tib_msr[!complete.cases(tib_msr), ],  data.frame( conc = c( new_concs[ ,1])))
tib_imp <- tib_imp[, !duplicated(colnames(tib_imp), fromLast = TRUE)] 

# combine predicted and measured data for verification purposes
# ============================================================

tib_new <- rbind(tib_imp , tib_msr[complete.cases(tib_msr), ])

plot(picg ~conc, data = tib_new, col=ifelse(tib_msr$pos %in% c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8"), 'red', 'black'),
    main = "predicted relationship Picogreen to DNA concentration", 
    sub = "(red dots are qbit quantified standards)",
    ylab = "fluorescence",
    xlab = "[ng/µl]")

# print concentrations in plate format 
# ===================================

tib_new$prow <- as.character(substring(tib_new$pos, 1, 1))
tib_new$pcol <- as.numeric(substring(tib_new$pos, 2, nchar(tib_new$pos)))
tib_new <- tib_new %>% arrange(prow) %>% arrange(pcol)

# matrix(tib_new$pos, nrow = 8, ncol = 12, dimnames = list( c("A","B","C","D","E","F","G","H"), c(seq(1, 12, 1))))
plate <- matrix( round(tib_new$conc, digits = 2), nrow = 8, ncol = 12, dimnames = list( c("A","B","C","D","E","F","G","H"), c(seq(1, 12, 1))))
plate <- replace(plate, which(plate < 0), NA)
plate

write.xlsx(plate, "/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate2_ng-ul.xlsx", 
  rowNames = TRUE, colNames = TRUE, overwrite = FALSE)

