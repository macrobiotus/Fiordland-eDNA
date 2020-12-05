library("tidyverse")
library("janitor")
library("reshape2")

denoised <- read_tsv("/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-sta.tsv")[-1, ]
denoised <- denoised %>% clean_names()
denoised <- denoised %>% mutate_at(vars(-sample_id), as.numeric)


denoised <-  melt(denoised, id="sample_id")
denoised <- denoised %>% dplyr::select(-one_of("percentage_of_input_passed_filter", "percentage_of_input_non_chimeric"))

ggplot (denoised, aes(x = sample_id, y = value, fill = variable)) + 
  geom_bar (stat="identity", position = position_dodge(width = 0.5)) + 
  theme_bw()