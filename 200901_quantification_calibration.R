# *****************************************************************************
# * associating fluorescence signal with DNA concentrations, lab book page 48 *
# *****************************************************************************
# 1-Sep-2020

# load packages
# =============
rm(list = ls(all.names = TRUE))
gc()

library("tidyverse")   # work using tibbles and ggplot
library("readxl")      # read excel sheets 
library("cowplot")     # multiple plots per page

# read-in data 
# ============


msr <- "/Users/paul/Documents/OU_eDNA/200128_lab_work/200901_qpcr_qbit_test_results_formatted_data_with_qbit.xlsx"
tib_msr <- read_excel(msr) %>%
  mutate_at(vars(Well, STCK_NG_UL, ASM_DIL_STCK_NG_UL, MSR_DIL_STCK_NG_UL, DIL_FCT, ASM_DIL_ASSY_NG_UL, MSR_DIL_ASSY_NG_UL, Delta_Rn), as.numeric) %>%
  mutate_at(vars(Well_Position, Content), as.factor)


# plot data 
# =========

p1 <- ggplot(tib_msr,aes(y = Delta_Rn, x = ASM_DIL_ASSY_NG_UL, colour = Content, shape = Content)) +
   geom_point() + geom_smooth(method = "lm") +
   coord_cartesian(xlim = c(0,0.5),ylim= c(0,1500000)) +
   ggtitle("fluorescence as function of \n assumed [c] (ng/µl) of 1:25 dilutions (actual assay of 50 µl) \n ") +
   xlab("ng/µl") + 
   ylab("∆R") +
   theme_bw() 
   
p2 <- ggplot(tib_msr,aes(y = Delta_Rn, x = MSR_DIL_ASSY_NG_UL, colour = Content, shape = Content)) +
   geom_point() + geom_smooth(method = "lm") +
   coord_cartesian(xlim = c(0,0.5),ylim= c(0,1500000)) + 
   ggtitle("fluorescence as function of \n measured [c] (ng/µl) of 1:25 dilutions (actual assay of 50 µl)") +
   xlab("ng/µl") + 
   ylab("∆R") +
   theme_bw() 


p3 <- ggplot(tib_msr,aes(y = Delta_Rn, x = ASM_DIL_STCK_NG_UL, colour = Content, shape = Content)) +
   geom_point() + geom_smooth(method = "lm") +
   coord_cartesian(xlim = c(0,25),ylim= c(0,1500000)) + 
   ggtitle("fluorescence as function of \n assumed [c] (ng/µl) of stock solution (1:2 dilution series)") +
   xlab("ng/µl") + 
   ylab("∆R") +
   theme_bw() 
   
p4 <- ggplot(tib_msr,aes(y = Delta_Rn, x = MSR_DIL_STCK_NG_UL, colour = Content, shape = Content)) +
   geom_point() + geom_smooth(method = "lm")  + 
   coord_cartesian(xlim = c(0,25),ylim= c(0,1500000)) + 
   ggtitle("fluorescence as function of \n measured [c] (ng/µl) of stock solution (1:2 dilution series)") +
   xlab("ng/µl") + 
   ylab("∆R") +
   theme_bw() 

plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
