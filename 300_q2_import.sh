#!/usr/bin/env bash

# 26.11.2020 - Paul Czechowski - paul.czechowski@gmail.com 
# ========================================================
# Import sequence files into Qiime 2.

# safety settings 
#   "-u" makes it an error to reference a non-existent environment variable 
#   "pipefail" makes things like misspeled-command raise errors.

set -eu
set -o pipefail


# Adjust base paths
# -----------------
if [[ "$HOSTNAME" != "Pauls-MacBook-Pro.local" ]] && [[ "$HOSTNAME" != "macmini-fastpost-1.staff.uod.otago.ac.nz" ]]; then
    
    # adjust terminal colours
    bold=$(tput bold)
    normal=$(tput sgr0)

    # location message 
    printf "${bold}$(date):${normal} Execution on remote is not implemnted...\n"

elif [[ "$HOSTNAME" == "Pauls-MacBook-Pro.local" ]]  || [[ "$HOSTNAME" == "macmini-fastpost-1.staff.uod.otago.ac.nz" ]]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on local...\n"
    trpth="/Users/paul/Documents/OU_eDNA"
    cores="2"
fi

# input file array - 12S
inpth[1]=`/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-02-0-1_L001-ds.12af8aae24bb42c3b49d3dd022926985/6102-02-0-1_S2_L001_R1_001.fastq.gz`
inpth[2]=`/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-02-0-1_L001-ds.12af8aae24bb42c3b49d3dd022926985/6102-02-0-1_S2_L001_R2_001.fastq.gz`
inpth[3]=`/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-01-0-1_L001-ds.30ba3924733f460c8b1b94602f109226/6102-01-0-1_S1_L001_R1_001.fastq.gz`
inpth[4]=`/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-01-0-1_L001-ds.30ba3924733f460c8b1b94602f109226/6102-01-0-1_S1_L001_R2_001.fastq.gz`

# output file array - 18S
otpth[1]=`/201126_preprocessing/qiime/6102-02-0-1_S2_L001_R1_001.qza`
otpth[2]=`/201126_preprocessing/qiime/6102-02-0-1_S2_L001_R2_001.qza`
otpth[3]=`/201126_preprocessing/qiime/6102-01-0-1_S1_L001_R1_001.qza`
otpth[4]=`/201126_preprocessing/qiime/6102-01-0-1_S1_L001_R2_001.qza`






# Run import script - adjust `i` starting number! 
# -----------------------------------------------
for ((i=1;i<=2;i++)); do
    
    # call import only if output file isn't already there
    if [ ! -f "$trpth"/"${otpth[$i]}.qza" ]; then
    
      # diagnostic message
      printf "${bold}$(date):${normal} Starting importing from \"$(basename "$trpth"/"${inpth[$i]}")\"...\n"
      qiime tools import \
        --type "SampleData[SequencesWithQuality]" \
        --input-path  "$trpth"/"${inpth[$i]}" \
        --output-path "$trpth"/"${otpth[$i]}" \
        --input-format "SingleEndFastqManifestPhred33" 2>&1 | tee -a "$trpth"/"Zenodo/Processing"/"$(basename ${otpth[$i]} .qza)_import_log.txt" || \
      printf "Import failed at "$(date)" on \"${otpth[$i]}\". \n"
    
    else
    
      # diagnostic message
      printf "${bold}$(date):${normal} Import available for \"$(basename "$trpth"/"${inpth[$i]}")\", skipping.\n"
  
    fi
done
