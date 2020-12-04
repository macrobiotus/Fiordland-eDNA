#!/usr/bin/env bash

# 03.12.2020 - Paul Czechowski - paul.czechowski@gmail.com 
# ========================================================
# Denoise 2 sequencing runs using DADA2

# safety settings 
#   "-u" makes it an error to reference a non-existent environment variable 
#   "pipefail" makes things like misspeled-command raise errors.

set -eu
set -o pipefail


# Adjust base paths
# -----------------
if [[ "$HOSTNAME" != "Pauls-MacBook-Pro.local" ]] && [[ "$HOSTNAME" != "macmini-fastpost-1.staff.uod.otago.ac.nz" ]]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on remote...\n"
    trpth="/workdir/pc683/OU_eDNA"
    cores="$(nproc --all)"
elif [[ "$HOSTNAME" == "Pauls-MacBook-Pro.local" ]]  || [[ "$HOSTNAME" == "macmini-fastpost-1.staff.uod.otago.ac.nz" ]]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on local...\n"
    trpth="/Users/paul/Documents/OU_eDNA"
    cores="2"
fi

# define input locations
# ----------------------
inpth[1]='201126_preprocessing/qiime/500_12S_single_end_import.qza'

# define output locations
# -----------------------
otpth_seq[1]='201126_preprocessing/qiime/600_12S_single_end_ee3-seq.qza'

otpth_tab[1]='201126_preprocessing/qiime/600_12S_single_end_ee3-tab.qza'

otpth_stat[1]='201126_preprocessing/qiime/600_12S_single_end_ee3-sta.qza'

otpth_vis[1]='201126_preprocessing/qiime/600_12S_single_end_ee3-vis.qzv'

# trimming parameters 18S - reads already filtered for Phred 25 
# --------------------------------------------------------------
# amplicon should be no more then 325bp 
trnc[1]='0'

# allowed no more then "1" errors in sequence **check log file for last applied value!!!**, default is "2"
eerr[1]='3'

# run script
# ----------


for ((i=1;i<=1;i++)); do

   # denoising
   qiime dada2 denoise-single \
      --i-demultiplexed-seqs "$trpth"/"${inpth[$i]}" \
      --p-trunc-len "${trnc[$i]}" \
      --p-max-ee "${eerr[$i]}" \
      --p-n-threads "$cores" \
      --o-table "$trpth"/"${otpth_tab[$i]}" \
      --o-representative-sequences "$trpth"/"${otpth_seq[$i]}" \
      --o-denoising-stats "$trpth"/"${otpth_stat[$i]}" \
      --verbose

    # export stats file for manual inspection and gnuplot
    qiime metadata tabulate\
      --m-input-file "$trpth"/"${otpth_stat[$i]}" \
      --o-visualization "$trpth"/"${otpth_vis[$i]}"
done
