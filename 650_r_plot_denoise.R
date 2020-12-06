library("tidyverse")
library("janitor")
library("reshape2")

denoised <- read_tsv("/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-vis.tsv")[-1, ]
denoised <- denoised %>% clean_names()
denoised <- denoised %>% mutate_at(vars(-sample_id), as.numeric)
denoised <- denoised %>% dplyr::select(-one_of("percentage_of_input_passed_filter", "percentage_of_input_non_chimeric"))
denoised <-  denoised %>% mutate(plate=substr(sample_id, 0, 3))

denoised <-  melt(denoised, id=c("sample_id", "plate"))

ggplot (denoised, aes(x = sample_id, y = value, fill = variable)) + 
  geom_bar (stat="identity", position = position_dodge()) +
  xlab("library name") + 
  ylab("read count") +
  facet_grid(rows = "plate") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5))

ggsave("200412_650_denoised_libraries.pdf", plot = last_plot(), 
         device = "pdf", path = "/Users/paul/Documents/OU_eDNA/201204_display_items",
         scale = 3, width = 200, height = 50, units = c("mm"),
         dpi = 500, limitsize = TRUE)
