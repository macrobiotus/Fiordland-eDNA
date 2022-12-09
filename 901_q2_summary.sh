#!/usr/bin/env bash

# 09.12.2022 - Paul Czechowski - paul.czechowski@gmail.com 
# ========================================================

# set -x
set -eu
set -o pipefail

# Adjust base paths
# -----------------
if [[ "$HOSTNAME" != "Pauls-MacBook-Pro.local" ]] && [[ "$HOSTNAME" != "macmini-fastpost.local" ]]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on remote...\n"
    trpth="/workdir/pc683/OU_eDNA"
    cores="$(nproc --all)"
elif [[ "$HOSTNAME" == "Pauls-MacBook-Pro.local" ]]  || [[ "$HOSTNAME" == "macmini-fastpost.local" ]]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
    printf "${bold}$(date):${normal} Execution on local...\n"
    trpth="/Users/paul/Documents/OU_eDNA"
    cores="2"
fi

# define relative input locations - Qiime files
# --------------------------------------------------------
inpth_map='201126_preprocessing/metadata/850_prep_q2_predictor-tab__metadata.tsv'
# new blast results from en-of-year 2022:
inpth_tax='201126_preprocessing/qiime/801_12S_single_end_ee3-seq_q2taxtable.qza'

# define relative input locations - sequence files
# -----------------------------------------------------------

# (https://stackoverflow.com/questions/23356779/how-can-i-store-the-find-command-results-as-an-array-in-bash)

# Fill sequence array using find 
inpth_seq_unsorted=()
while IFS=  read -r -d $'\0'; do
    inpth_seq_unsorted+=("$REPLY")
done < <(find "$trpth/201126_preprocessing/qiime" \( -name '600_12S_single_end_ee3-seq.qza' \) -print0)

# for debugging - print unsorted sequences
# printf '%s\n'
# printf '%s\n' "${inpth_seq_unsorted[@]}"

# Sort array 
IFS=$'\n' inpth_seq=($(sort <<<"${inpth_seq_unsorted[*]}"))
unset IFS

# for debugging  - print sorted sequences - ok!
# printf '%s\n'
# printf '%s\n' "${inpth_seq[@]}"

# define relative input locations - feature tables
# ------------------------------------------------

# Fill table array using find 
inpth_tab_unsorted=()
while IFS=  read -r -d $'\0'; do
    inpth_tab_unsorted+=("$REPLY")
done < <(find "$trpth/201126_preprocessing/qiime" \( -name '600_12S_single_end_ee3-tab.qza' \) -print0)

# for debugging -  print unsorted tables
# printf '%s\n'
# printf '%s\n' "${inpth_tab_unsorted[@]}"

# Sort array 
IFS=$'\n' inpth_tab=($(sort <<<"${inpth_tab_unsorted[*]}"))
unset IFS

# for debugging -  print sorted tables - ok!
# printf '%s\n'
# printf '%s\n' "${inpth_tab[@]}"

# define relative output locations - feature tables
# otpth_tabv='Zenodo/Qiime/080_18S_denoised_tab_vis.qzv'
# otpth_seqv='Zenodo/Qiime/080_18S_denoised_seq_vis.qzv'
# otpth_bplv='Zenodo/Qiime/080_18S_denoised_tax_vis.qzv'

# loop over filtering parameters, and corresponding file name names additions
for i in "${!inpth_seq[@]}"; do

  # check if files can be mathced otherwise abort script because it would do more harm then good
  seqtest="$(basename "${inpth_seq[$i]//-seq/}")"
  tabtest="$(basename "${inpth_tab[$i]//-tab/}")"
  
  # for debugging
  # echo "$seqtest"
  # echo "$tabtest"
  
  
  if [ "$seqtest" == "$tabtest" ]; then
    echo "Sequence- and table files have been matched, continuing..."
  
    # get input sequence file name - for debugging 
    # echo "${inpth_seq[$i]}"
    
    # get input table file name  - for debugging
    # echo "${inpth_tab[$i]}"
    
    directory="$(dirname "$inpth_seq[$i]")"
    seq_file_tmp="$(basename "${inpth_seq[$i]%.*}")"
    seq_file_name="${seq_file_tmp:4}"
    
    tab_file_tmp="$(basename "${inpth_tab[$i]%.*}")"
    tab_file_name="${tab_file_tmp:4}"
    
    plot_file_temp="$(basename "${inpth_seq[$i]//_seq/}")"
    plot_file_temp="${plot_file_temp:4}"
    plot_file_name="${plot_file_temp%.*}"
    
    extension=".qzv"
    
    # check string construction - for debugging
    # echo "$seq_file_name"
    # echo "$tab_file_name"
    # echo "$plot_file_name"
    
    seq_file_vis_path="$directory/901_$seq_file_name""$extension"
    tab_file_vis_path="$directory/901_$tab_file_name""$extension"
    plot_file_vis_path="$directory/901_$plot_file_name"_barplot"$extension"
    
    # check string construction - for debugging
    # echo "$seq_file_vis_path"
    # echo "$tab_file_vis_path"
    # echo "$plot_file_vis_path"
    
    # Qiime calls
    qiime feature-table tabulate-seqs \
      --i-data "${inpth_seq[$i]}" \
      --o-visualization "$seq_file_vis_path" \
      --verbose

    qiime feature-table summarize \
      --m-sample-metadata-file "$trpth"/"$inpth_map" \
      --i-table "${inpth_tab[$i]}" \
      --o-visualization "$tab_file_vis_path" \
      --verbose
 
    qiime taxa barplot \
      --m-metadata-file "$trpth"/"$inpth_map" \
      --i-taxonomy "$trpth"/"$inpth_tax" \
      --i-table "${inpth_tab[$i]}" \
      --o-visualization "$plot_file_vis_path" \
      --verbose

  else
  
    echo "Sequence- and table files can't be matched, aborting."
    exit
  
  fi
  
done
