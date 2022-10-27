# In-silico PCR to check MiFish primer specificty of mammal primers
# ----------------------------------------------------------

# Cannon, M. V. et al. In silico assessment of primers for eDNA studies
#   using PrimerTree and application to characterize the biodiversity
#   surrounding the Cuyahoga River. Sci. Rep. 6, 22908; doi: 10.1038/srep22908 (2016).

# load packages
library("curl")
library("primerTree") # https://github.com/MVesuviusC/primerTree
                      # compile failed from source as per https://github.com/MVesuviusC/primerTree/blob/master/README.md
                      # using precompiled version http://clustal.org/omega/
library("tidyverse")


# as from /Users/paul/Documents/OU_eDNA/200901_scripts/200_r_metadata_management.R line 107

# Mifish U fwd / rev
# GTCGGTAAAACTCGTGCCAGC CATAGTGGGGTATCTAATCCCAGTTTG

# Mifish E fwd / rev
# GTTGGTHAATCTCGTGCCAGC CATAGTAGGGTATCTAATCCTAGTTTG


# run remote primer BLAST - **chnage API key prior to publication
spp_result_mu <- search_primer_pair(name='Mifish_U', 'GTCGGTAAAACTCGTGCCAGC', 'CATAGTGGGGTATCTAATCCCAGTTTG',
                                 api_key = "a634c6e9c96c3859bca27a2771f6d2872f08",
                                 clustal_options=c(exec='/Users/paul/Biobin/clustalo/clustal-omega-1.2.3-macosx'),
                                 num_aligns=500, total_primer_specificity_mismatch=3
                                )
saveRDS(spp_result_mu, "/Users/paul/Documents/OU_eDNA/201028_Robjects/210421_in_silico_Mifish_U.Rds")
saveRDS(spp_result_mu, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210421_in_silico_Mifish_U.Rds")


spp_result_me <- search_primer_pair(name='Mifish_E_mod', 'GTTGGTHAATCTCGTGCCAGC', 'CATAGTAGGGTATCTAATCCTAGTTTG',
                                 api_key = "a634c6e9c96c3859bca27a2771f6d2872f08",
                                 clustal_options=c(exec='/Users/paul/Biobin/clustalo/clustal-omega-1.2.3-macosx'),
                                 num_aligns=500, total_primer_specificity_mismatch=3
                                )
saveRDS(spp_result_me, "/Users/paul/Documents/OU_eDNA/201028_Robjects/210421_in_silico_Mifish_U.Rds")
saveRDS(spp_result_me, "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210421_in_silico_Mifish_U.Rds")

# plot primer BLASt results
plot(spp_result_mu, ranks = "class", type = "radial", size = 1, legend_cutoff = 100)
ggsave("210421_Mifish_U_in_silico.pdf", plot = last_plot(),
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures/",
       scale = 1, width = 300, height = 300, units = c("mm"),
       dpi = 500, limitsize = TRUE)

plot(spp_result_me, ranks = "class", type = "radial", size = 1, legend_cutoff = 100)
ggsave("210421_Mifish_E_in_silico.pdf", plot = last_plot(),
       device = "pdf", path = "/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/figures/",
       scale = 1, width = 300, height = 300, units = c("mm"),
       dpi = 500, limitsize = TRUE)
