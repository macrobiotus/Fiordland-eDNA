# *****************************************************************************
# * associating fluorescence signal with DNA concentrations, lab book page 48 *
# *****************************************************************************
# 27-Oct-2022

# load packages
# =============
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse")   # work using tibbles and ggplot
library("readxl")      # read excel sheets 
library("openxlsx")    # write Excel tables

# read-in and format data 
# =======================

# HS Standards
# **measured** qbit values from lab book for plate positions H1  - H8 
#   qbit_values <- c(0, 9.1, 4.64, 2.21, 1.15, 0.58, 0.282, 0.148)
# **targeted** qbit values from lab book for plate positions H1  - H8 
qbit_values <- c(0, 10, 5, 2.5, 1.25, 0.625, 0.3125, 0.156)

# BR Standards
# **targeted** qbit values from lab book for plate positions H1  - H8 
qbit_values <- c(0, 100, 50, 25, 12.5, 6.25, 3.125, 1.56)

# pathnames to Excel exports of plate quantifications
paths <- list(
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Eplate1_20200914 234405.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Eplate2_20200914 234524.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Uplate1_20200914 234210.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/Uplate2_20200914 234305.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200917_Eplate1_2nd_20200917 041057.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200917_Eplate2_20200917 041000.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200917_Uplate1_20200917 041151.xlsx",
  "/Users/paul/Documents/OU_eDNA/200128_lab_work/200917_Uplate2_20200917 041123.xlsx"
  )

tib_msr <- read_excel(paths[[8]], sheet = 2, col_names = TRUE, skip = 22, .name_repair = "unique") %>%
  rename ("well" = "Well", "pos" = `Well Position`, "cycl" = `Cycle Number`, "picg" =  "SYBR")

# correct and set columns 
tib_msr$cycl <- NULL 
tib_msr$conc <- NA

# add qbit values 
tib_msr$conc[which(tib_msr$pos %in% c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8"))] <- qbit_values

# calculate variables
# ===================

# see https://stats.stackexchange.com/questions/176595/simple-log-regression-model-in-r/176612

# plot quantification data
plot(conc ~ picg, data = tib_msr[complete.cases(tib_msr), ])

# fit model
fit <- lm(conc ~ picg, data = tib_msr[complete.cases(tib_msr), ] )

# check fit - good fit
summary(fit)

# plot quantification data - with ggplot2
ggplot( tib_msr[complete.cases(tib_msr), ] , aes(x = conc, y = picg)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))

# predict new values and plot 
new_concs = predict(fit, newdata <- tib_msr[!complete.cases(tib_msr), ], interval="confidence")
plot(new_concs[ ,1] ,  tib_msr[!complete.cases(tib_msr), ]$picg, 
  main = "predicted relationship pico-green to DNA concentration",
  ylab = "fluorescence",
  xlab = "[ng/µl]")
  matlines(new_concs, tib_msr[!complete.cases(tib_msr), ]$picg , lwd=1)
    
# backt-transform concentration values - not used  - formerly contained function call 
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

# write.xlsx(plate, "/Users/paul/Documents/OU_eDNA/200128_lab_work/200917_Uplate2_ng-ul.xlsx", 
#   rowNames = TRUE, colNames = TRUE, overwrite = FALSE)

write.xlsx(plate, "/Users/paul/Documents/OU_eDNA/200128_lab_work/221027_Uplate2_ng-ul.xlsx", 
  rowNames = TRUE, colNames = TRUE, overwrite = FALSE)
