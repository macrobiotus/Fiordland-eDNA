#!/usr/bin/env bash

#SBATCH --job-name      BLAST
#SBATCH --time          05:00:00
#SBATCH --mem           200G  # 50 GB for running - 150 GB for database
#SBATCH --ntasks        1
#SBATCH --cpus-per-task 36    # half a node

# 29.10.2022 - Paul Czechowski - paul.czechowski@gmail.com 
# ========================================================

# * Blast `fasta` against NCBI's nt database. See
#   https://stackoverflow.com/questions/45014279/running-locally-blastn-against-nt-db-thru-python-script.
# * For array handling refer to
#   https://stackoverflow.com/questions/23356779/how-can-i-store-find-command-result-as-arrays-in-bash
#   https://stackoverflow.com/questions/7442417/how-to-sort-an-array-in-bash)
# * For execution on NESI see
#   https://support.nesi.org.nz/hc/en-gb/articles/208619807-BLAST
#
#   For jobs which need more than 24 CPU-hours, eg: those that use large
#   databases (> 10 GB) or large amounts of query sequence (> 1 GB), or slow
#   BLAST searches such as classic blastn (blastn -task blastn).
# 
#   This script copies the BLAST database into the per-job temporary
#   directory $TMPDIR before starting the search. Since compute nodes do not
#   have local disks, this database copy is in memory, and so must be
#   allowed for in the memory requested by the job. As of mid 2021 that is
#   125 GB for the nt database, 124 GB for refseq_protein. 
# * For parssing use MEGAN and `xml` format, '-outfmt 5' 
# 
# * safety settings 
#   "-u" makes it an error to reference a non-existent environment variable 
#   "pipefail" makes things like misspeled-command raise errors.

set -eu
set -o pipefail

# * adjust terminal colour --- 

bold=$(tput bold)
normal=$(tput sgr0)

# * set machine-specific environment ---

if [[ "$HOSTNAME" != "Pauls-MacBook-Pro.local" ]] && [[ "$HOSTNAME" != "mbpro.local" ]]; then
    
  # diagnostic message
  printf "\n${bold}$(date):${normal} Blasting on remote...\n"
  
  # load Blast Modules on NESI 
  module purge
  module load BLAST/2.10.0-GCC-9.2.0
  module load BLASTDB/2021-05
  
  # check available db versions using the command below and alter the date in the previous line accordingly
  ls "$BLASTDB"

  # copy BLASTDB to RAM
  cp "$BLASTDB"/taxdb* "$TMPDIR"/ 
  cp "/nesi/nobackup/uoo03176/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ" "$TMPDIR"/ 
  
  # export path and other variables
  export BLASTDB="$TMPDIR"
  export BASEPATH="/nesi/project/uoo03176/OU_eDNA"
  export NOBACKUP="/nesi/nobackup/uoo03176/OU_eDNA"
  export CPUS="$SLURM_CPUS_PER_TASK"
    
elif [[ "$HOSTNAME" == "mbpro.local" ]]  || [[ "$HOSTNAME" == "macmini-fastpost.local" ]]; then
    
  # diagnostic message
  printf "\n${bold}$(date):${normal} Blasting on local...\n"
  
  # export path and other variables
  export BLASTDB="/Volumes/HGST1TB/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ"
  export BASEPATH="/Users/paul/Documents/OU_eDNA"
  export NOBACKUP="/Volumes/HGST1TB/Users/paul/Documents/OU_eDNA"
  export CPUS="1"
  
fi

# * Copy files to scratch ---
  
if [ -d "$BASEPATH" ]; then
  
  # diagnostic message
  printf "\n${bold}$(date):${normal} Syncing files from \"$BASEPATH\" to \"$NOBACKUP\"...\n"
  
  rsync -azi --delete-after --progress "$BASEPATH" "$(dirname "$NOBACKUP")" 
  
fi 

# * find all .fasta files incl. their paths and put them into an array

# diagnostic message
printf "\n${bold}$(date):${normal} Finding fasta files...\n"

inpth_seq_unsorted=()
while IFS=  read -r -d $'\0'; do
    inpth_seq_unsorted+=("$REPLY")
done < <(find "$NOBACKUP/201126_preprocessing/qiime" -type f \( -name "*.fasta.gz" \) -print0)

# * sort array 
IFS=$'\n' inpth_seq=($(sort <<<"${inpth_seq_unsorted[*]}"))
unset IFS

# * print sorted input filenames
# printf '%s\n' "${inpth_seq[@]}"

# * loop over array of fasta files, create result file, call Blast

for fasta in "${inpth_seq[@]}";do
  
  # create Blast results file name
  filename=$(dirname "$fasta")
  src_dir=$(basename "$fasta")
  tmp_file="752${src_dir:3}"
  tgt_file="${tmp_file%%.*}_blast-221030_MetaFishLibNZ.xml"

  # test Blast results file name creation 
  # printf "$fasta\n"
  # printf "$NOBACKUP"/201126_preprocessing/blast/"$tgt_file.gz\n"
  # exit
  
  if [ ! -f "$NOBACKUP"/201126_preprocessing/blast/"$tgt_file".gz ]; then
    
    # Diagnostic message
    printf "\n${bold}$(date):${normal} Querying \"$fasta\" against \"$BLASTDB\"...\n" && \
    gzip -dc "$fasta" | \
      blastn \
        -db "221030_MetaFishLibNZ" \
        -task blastn \
        -evalue 1e-10 \
        -max_hsps 5 \
        -outfmt 5 \
        -max_target_seqs 5 \
        -out "$NOBACKUP"/201126_preprocessing/blast/"$tgt_file" \
        -num_threads "$CPUS" \
        -qcov_hsp_perc 95 \
        -perc_identity 75 && \
      printf "\n${bold}$(date):${normal} Blast finished writing to \"$NOBACKUP/201126_preprocessing/blast/$tgt_file\".\n" || \
      { printf "\n${bold}$(date):${normal} Blastn failed on \"$fasta\". \n" ; exit 1; }
    printf "\n${bold}$(date):${normal} Compressing \"$NOBACKUP/201126_preprocessing/blast/$tgt_file\".\n"
    
    gzip "$NOBACKUP"/201126_preprocessing/blast/"$tgt_file"
  
  fi
  
done

# * Copy files from scratch ---

if [ -d "$NOBACKUP" ]; then

  # diagnostic message
  printf "\n${bold}$(date):${normal} Syncing files from \"$NOBACKUP\" to \"$BASEPATH\"...\n"
  
  rsync -avzui --progress --delete-after --delete-after --progress "$NOBACKUP" "$(dirname "$BASEPATH")" 
  
fi

exit
