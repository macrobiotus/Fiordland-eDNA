#!/usr/bin/env bash

# 01.02.2021 - Paul Czechowski - paul.czechowski@gmail.com 
# ========================================================
# Count sequences in demultiplexed fastq files - move if empty 

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
    printf "${bold}$(date):${normal} Execution not implemneted...\n"
    exit
    
    trpth="/workdir/pc683/OU_eDNA"
    cores="$(nproc --all)"
    
    export PYTHONPATH=/programs/cutadapt-2.1/lib/python3.6/site-packages:/programs/cutadapt-2.1/lib64/python3.6/site-packages
    export PATH=/programs/cutadapt-2.1/bin:$PATH

elif [[ "$HOSTNAME" == "Pauls-MacBook-Pro.local" ]]  || [[ "$HOSTNAME" == "macmini-fastpost-1.staff.uod.otago.ac.nz" ]]; then

    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on local...\n"

    trpth="/Users/paul/Documents/OU_eDNA"
    cores="2"
    
fi

# define working directory for portability and easy editing
wd="$trpth/201126_preprocessing/cutadapt"


# as fastest solution taken from:
#   https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-
#   number-of-reads-and-number-of-bases-in-a-fastq-file

# loop over files - copy files without reads and corresponding logs to sub directory
for filename in $wd/*.fastq.gz; do

  [ -e "$filename" ] || continue
  
  
  # diagnostic message
  printf "${bold}$(date):${normal} Checking \"$(basename $filename)\"...\n"

  
  fix_base_count() {
    local counts=($(cat))
    echo "${counts[0]} $((${counts[1]} - ${counts[0]}))"
  }
  
  read_count=$(gzip -dc "$filename" | awk 'NR % 4 == 2' |  wc -cl | fix_base_count | cut -d ' ' -f 1)
  
  # printf "...read count is $read_count.\n"
  
  if [ "$read_count" -eq "0" ]; then
    
    printf "...moving file with $read_count read count.\n"
    
    log_file=$(echo "$filename" | sed 's+.fastq.gz+.log+g')
    
    mkdir -p "$wd/350_empty_fastqs"
    
    mv "$filename" "$wd/350_empty_fastqs"
    mv "$log_file" "$wd/350_empty_fastqs"
  fi

done
