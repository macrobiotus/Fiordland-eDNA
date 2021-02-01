#!/usr/bin/env bash

# 03.12.2020 - Paul Czechowski - paul.czechowski@gmail.com 
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
    
    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on remote...\n"
    trpth="/workdir/pc683/OU_eDNA"
    cores="2"

elif [[ "$HOSTNAME" == "Pauls-MacBook-Pro.local" ]]  || [[ "$HOSTNAME" == "macmini-fastpost-1.staff.uod.otago.ac.nz" ]]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on local...\n"
    trpth="/Users/paul/Documents/OU_eDNA"
    cores="2"
fi

manifest='201126_preprocessing/metadata/400_create_qiime_manifest__manifest_v2.txt'
otpth='201126_preprocessing/qiime/500_12S_single_end_import.qza'


# Run import script
# -----------------
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path  "$trpth"/"$manifest" \
  --output-path "$trpth"/"$otpth" \
  --input-format SingleEndFastqManifestPhred33
