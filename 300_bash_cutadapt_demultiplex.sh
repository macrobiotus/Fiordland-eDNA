#!/usr/bin/env bash

# 05.12.2020 - Paul Czechowski - paul.czechowski@gmail.com 
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

# barcode file generated by `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R` on or after 1-12-2020
barcodes="/201126_preprocessing/metadata/200_cutadapt_barcode_input.txt" 

# input file array - 12S
inpth[1]="$trpth/201120_sequence_data/OG6304-211308098/FASTQ_Generation_2020-11-19_01_33_21Z-343990651/6304-P1-0-1_L001-ds.21f7e23247f540d4b3dfe0004841c2b5/6304-P1-0-1_S1_L001_R1_001.fastq.gz"


# Run import script - adjust `i` starting number! 
# -----------------------------------------------


for ((i=1;i<=1;i++)); do
    
  # creating_partial_output file name
  outdir="$trpth"/201126_preprocessing/cutadapt
  namest="$outdir"/"$(basename "$trpth"/"${inpth[$i]}" .fastq.gz)"
   
  # for debugging only
  # echo "$namest"
  
  # for debugging only
  # echo "$outfile"
      
  # echo "$revbc"
  # echo "$revbcrc"
    
  # diagnostic message
  printf "${bold}$(date):${normal} Demultiplexing \"$(basename "${inpth[$i]}")\"...\n"

  # step through barcode file
  while IFS=" " read -r smpl_nm fwd_bc rev_bc rev_bc_rc
    do
  
      touch "$trpth"/201126_preprocessing/cutadapt/300_bash_cutadapt_demultiplex_input.fastq.gz
      touch "$trpth"/201126_preprocessing/cutadapt/300_bash_cutadapt_demultiplex_ouput.fastq.gz
  
      # generate output file name
      outfile="$namest"__"$smpl_nm".fastq.gz

      # for debugging only
      # echo "$outfile"
        
      # call import only if output file isn't already there
      if [ ! -f "$outfile" ]; then

        FILESIZE=$(wc -c  "${inpth[$i]}" | awk '{print $1}')
        
        # diagnostic message
        printf "${bold}$(date):${normal} ... demultiplexing \"$(basename "${inpth[$i]}")\" with $FILESIZE bytes...\n"
        
        cutadapt -a "$fwd_bc"';required...'"$rev_bc_rc"';required' -o "$namest"__"$smpl_nm".fastq.gz "${inpth[$i]}" --untrimmed-output "$trpth"/201126_preprocessing/cutadapt/300_bash_cutadapt_demultiplex_ouput.fastq.gz  -e 0 --no-indels -m 85 -M 250 --revcomp \
        2>&1 | tee "$trpth"/201126_preprocessing/cutadapt/"$smpl_nm"_cutlog.txt
        
        # move output file to input file so as to not to read and write from the same file at once 
        mv "$trpth"/201126_preprocessing/cutadapt/300_bash_cutadapt_demultiplex_ouput.fastq.gz "$trpth"/201126_preprocessing/cutadapt/300_bash_cutadapt_demultiplex_input.fastq.gz
        
        # WILL BREAK IF [1] CHNAGES: redefine input path so that left over sequences in temp file are used instead of input files
        # to fix copy inpth[1] at loop begin into new variable and use that ne variable throughout script
        inpth[1]="$trpth"/201126_preprocessing/cutadapt/300_bash_cutadapt_demultiplex_input.fastq.gz
        
      else
     
        # diagnostic message
        printf "${bold}$(date):${normal} Import available for \"$(basename "${inpth[$i]}")\", skipping.\n"
   
      fi
    
  done < "$trpth"/"$barcodes"

done

#     # call cutadapt with parallel    
#     parallel --jobs 0 -a "$trpth"/"$barcodes" --colsep ' ' \
#       cutadapt -a {2}'\;required...'{4}'\;required' -o "$namest"__{1}.fastq.gz "$trpth""${inpth[$i]}" --discard-untrimmed -e 0 --no-indels -m 1
  
