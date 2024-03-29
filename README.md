# Fiordland project utility scripts

## **1-Sep-2020**

  * creating `/Users/paul/Documents/OU_eDNA/200901_scripts/200901_quantification_calibration.R`
    * associating fluorescence signal with DNA concentrations, lab book page 48
    * using file `/Users/paul/Documents/OU_eDNA/200128_lab_work/200901_qpcr_qbit_test_results_formatted_data_with_qbit.xlsx`
    * created file `/Users/paul/Documents/OU_eDNA/200128_lab_work/200901_qpcr_qbit_test_results_regressions.pdf`
  * commit `e02097aac4f7cb65d800f267d2e7b7cf397aff7`

## **16-Sep-2020**

  * creating `/Users/paul/Documents/OU_eDNA/200901_scripts/200916_quantification_analysis.R`

## **17-Sep-2020**

  * some unneeded edits on `/Users/paul/Documents/OU_eDNA/200901_scripts/200916_quantification_analysis.R`
  * needs to be updated for new standards
  * using HR reads, as they seem to have performed better:
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate1_HS_ng-ul.xlsx`
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate2_HS_ng-ul.xlsx`
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate1_HS_ng-ul.xlsx`
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate2_HS_ng-ul.xlsx`
  * commit `e277b04a54b5599ad3ce94ad0a7be375335e211d`

## **28-Oct-2020** - creating `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`

  * writing objects to `/Users/paul/Documents/OU_eDNA/201028_Robjects`
  * Read cells from
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200907_plate_layouts.xlsx`
    * `/Users/paul/Documents/OU_eDNA/191031_primers/200302_Bunce_et_al_0000_MiFishEmod_single_step_primers.xlsx`
  * output appropriate format for
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/201028_Otago_Genomics_sample_info.xlsx`
    * and later possibly the mapping file
  * next: spread and subset for Otago Genomics - **pending**
  * commit `c51356eb489093cc333cd4c5024d3d535ce362a2`

## **29-Oct-2020** - continuing `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`

  * finished compiling sample data for Otago Genomics
  * commit `64b9cf9f69114a0b3f6d038368e7a6cc5bcb64aa`

## **26-Nov-2020** - getting metadata for sequence deconvolution

  * use script `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`
  * use data files as found via `find /Users/paul/Documents/OU_eDNA/ -name "*.fastq.gz"`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG5525-205458253/FASTQ_Generation_2020-10-22_16_08_01Z-332354023/5525-01-00-01_L001-ds.0ca0214f7d8c472e803c372dc2541ff6/5525-01-00-01_S1_L001_R2_001.fastq.gz`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG5525-205458253/FASTQ_Generation_2020-10-22_16_08_01Z-332354023/5525-01-00-01_L001-ds.0ca0214f7d8c472e803c372dc2541ff6/5525-01-00-01_S1_L001_R1_001.fastq.gz`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG5525-205458253/FASTQ_Generation_2020-10-22_16_08_01Z-332354023/5525-02-00-01_L001-ds.43cc4eb2a382470d85a8773b0b895cba/5525-02-00-01_S2_L001_R2_001.fastq.gz`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG5525-205458253/FASTQ_Generation_2020-10-22_16_08_01Z-332354023/5525-02-00-01_L001-ds.43cc4eb2a382470d85a8773b0b895cba/5525-02-00-01_S2_L001_R1_001.fastq.gz`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-02-0-1_L001-ds.12af8aae24bb42c3b49d3dd022926985/6102-02-0-1_S2_L001_R1_001.fastq.gz`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-02-0-1_L001-ds.12af8aae24bb42c3b49d3dd022926985/6102-02-0-1_S2_L001_R2_001.fastq.gz`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-01-0-1_L001-ds.30ba3924733f460c8b1b94602f109226/6102-01-0-1_S1_L001_R1_001.fastq.gz`
    * `/Users/paul/Documents/OU_eDNA/201120_sequence_data/OG6102-192920732/FASTQ_Generation_2020-08-28_16_46_04Z-304919615/6102-01-0-1_L001-ds.30ba3924733f460c8b1b94602f109226/6102-01-0-1_S1_L001_R2_001.fastq.gz`
  * installing `qiime2-2020.8` in `conda` environment
  * importing as per `https://docs.qiime2.org/2020.8/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq`
    * section `Multiplexed single-end FASTQ with barcodes in sequence`
  * copied for adjusting `/Users/paul/Documents/OU_eDNA/200901_scripts/100_q2_import.sh`
    * next - check in Qiime documentation and forum - to develop demultiplexing strategy
      * reverse complementing
      * dual barcodes single end reads
  * commit `425584e6393954bd7f209d40121ac5cb8b7345fd`

## **30-Nov-2020** - getting metadata for sequence deconvolution

  * tentatively using Qiime plugin for deconvolution 
    * see ` https://docs.qiime2.org/2020.8/plugins/available/cutadapt/demux-paired/`
  * first needed is metadata file
    * example metadata file at `/Users/paul/Documents/OU_eDNA/201126_script_scratch/sample-metadata_example.tsv`
  * continuing work on `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`
    * included pooling information from lab records
    * requesting new BRUV observation file as alternative to `/Users/paul/Documents/OU_eDNA/191213_field_work/200520_MH_bruv_data.xlsx`
    * as metadata reading in `/Users/paul/Documents/OU_eDNA/191213_field_work/201130_sample_overview_updated.xlsx` 
      * as alternative to `/Users/paul/Documents/OU_eDNA/191213_field_work/200321_sample_overview.xlsx`, in which primer allocations may be outadted
   * draft version of metadata file is generated by `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`
      * not saved yet - **needs formatting of barcode column matching template file** found at `/Users/paul/Documents/OU_eDNA/201126_script_scratch/sample-metadata_example.tsv`
      * also check concentration values in `big_table` object   
   * for work beyond this Qiime script templates can be found in `/Users/paul/Documents/OU_pcm_eukaryotes/Github`
   * all processing script now located at `/Users/paul/Documents/OU_eDNA/200901_scripts`
   * commit `33bc86967048d00c99ad35db5234b3e979a64d5`

## **01-Dec-2020** - renaming files - getting metadata for sequence deconvolution

  * settling on `cutadapt`, version 3.0 with linked adapters for demultiplexing
    * see `https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters-combined-5-and-3-adapter`
    * see `https://unix.stackexchange.com/questions/260840/read-columns-from-file-into-separate-variables`
  * get three column file from  `/Users/paul/Documents/OU_eDNA/200901_scripts/200_sample_matadata_managment.R` - **beta**
    * write to `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata` - **beta**
  * parse three column file with cutadapt for deconvolution - **in progress**
    * wrote `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/200_cutadapt_barcode_input.txt` - **done**
      * checked with primer order sheet - all there as in `/Users/paul/Documents/OU_eDNA/191031_primers/200504_idt_template_plate_filled_fjordland.xls`
    * use cutadapt for deconvolution - **pending**
      * starting `/Users/paul/Documents/OU_eDNA/200901_scripts/300_conda_cutadapt_demultiplex.sh` - **running**
      * using cutadapt parameters `--revcomp -e 0` - **running**
      * if unsuccessful:
        * trim only by index and keep reads pooled be primer U/E - **check results first**
        * requires adjustment of `/Users/paul/Documents/OU_eDNA/200901_scripts/200_sample_matadata_managment.R`  - **check results first**
  * commit `acc6ce02a13594827313fbce8c71f5c6539a7842`
  * corrected for correct reverse complemneting
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/200_sample_matadata_managment.R`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/300_conda_cutadapt_demultiplex.sh`
  * commit `0d17a836ef51aa9c2fa3f0f68f109ffd3696226`
  * added length filtering - **check files sizes, they are suspiciously similar** - **remove PhiX?** - **do not run multiple instances?**

## **03-Dec-2020** - creating manifest file for Qiime imports

  * demultiplexing has finished `/Users/paul/Documents/OU_eDNA/200901_scripts/300_conda_cutadapt_demultiplex.sh` - has finished
    * files at: `/Users/paul/Documents/OU_eDNA/201126_preprocessing/cutadapt`
  * updated `/Users/paul/Documents/OU_eDNA/200901_scripts/200_sample_matadata_managment.R`
    * for manifest file generation exports all data as R objcets, notably also
      * all metadata combined: `/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__big_table.Rdata`
      * object used for demultiplexing information `/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__demux_table.Rdata`
        * which is R object corresponding to `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/200_cutadapt_barcode_input.txt`
  * created script to create manifest file for Qiime: `/Users/paul/Documents/OU_eDNA/200901_scripts/400_create_qiime_manifest.R`
    * manifest at: `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/400_create_qiime_manifest__manifest.txt`
  * commit `36b58b2bb250f2ec0a6a8d1816303ba24eba946e`
  * copied and updated transport scripts for denoising and blasting using Cornell cluster (including commit)
    * `/Users/paul/Documents/OU_eDNA/201203_cu_transport`
  * pushing files to Cornell cluster

## **04-Dec-2020** - denoising finished, continuing processing

  * exported denoising stats using `qiime tools view /Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-vis.qzv`
  * adjusting and running `/Users/paul/Documents/OU_eDNA/200901_scripts/650_gnu_plot_denoise.gnu` - **aborted**
  * coding and running new plotting script `/Users/paul/Documents/OU_eDNA/200901_scripts/650_plot_denoise.R` -  **ok**
  * after denoising, to save space, erasing superflous `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/500_12S_single_end_import.qza`
  * exporting denoised sequences for blasting, before renaming them:
    * `qiime tools export  --input-path /Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-seq.qza --output-path /Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/`
    * `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-seq.fasta.gz`
  * expected peak length post primer and adapter trim: ~310 bp - 75 bp (fwd) - 58 bp (rev) = 177
  * observed length in Geneious import: mode about 255 bp for `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-seq.fasta.gz`
  * re-trimming primers in Geneious: shift distribution mode to about 240 bp - but wrong trimming
    * starting but unfinished `/Users/paul/Documents/OU_eDNA/200901_scripts/5550_q2_cutadapt.sh`
  * rechecking visually - only spurious remnants of barcode and primmer, leaving original file untouched
  * adjusting blast script for cluster and testing locally: `/Users/paul/Documents/OU_eDNA/200901_scripts/750_bash_fasta_blast.sh`
  * commit `49e02fd2147ca9a671a46e1114b11eef5d97238` and upload to cluster

## **05-Dec-2020** - redo cutadapt

  * taxonomy assignment yielded mostly off traget sequences longer then expected 177bp
  * restarting from `cutadapt`, without length filtering, using no filter at all to as to properly parametrize `dada2`
  * `cutadapt` and `parallel` combined in new script revisions
  * commit `d251d90228afe9bd2dc7c793fa95741114a01a75`
  * pulling data from cluster
  * forgot to remove 0 length sequences in cutadapt - call need to be corrected for import
  * redo on cluster - adjusted all path names
  * commit before cluster and push

## **06-Dec-2020**

  * length distribution post denoise finally ok
  * checked on cluster via `cat file.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,10) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sed '/^$/d' | awk '{print $2}' | sort | uniq -c | awk '{print $1,$2}'`
  * starting Blast before copying out files
  * erasing qiime import files (13 GB) and cutadapt output (18 GB)
  * pulling to local
  * checked files and corrected filenames
  * commit `cdd78aa2c32f31b04d258ef1ba4977eca53b81e`

## **07-Dec-2020**

  * adjusted and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/800_r_get_q2_tax-tab.r`
    * wrote file `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/800_12S_single_end_ee3-seq_q2taxtable.tsv`
  * ***next** - get a suitable metadata file
    * read in file `/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__big_table.Rdata` - **ok**
    * by adjusting `/Users/paul/Documents/OU_eDNA/200901_scripts/850_r_prep_q2_predictor-tab.r` - **pending**
    * commit `68f58558ca4b3f5100ea9e2814f10175b0860fc`

## **08-Dec-2020**

  * finished first version of `/Users/paul/Documents/OU_eDNA/200901_scripts/850_r_prep_q2_predictor-tab.r`
  * exported `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/850_prep_q2_predictor-tab__metadata.tsv`
  * commit `9cbef38cc355336abb8e2c4cddfb1f88e3fe8c19`
  * importing taxonomy file:
        `qiime tools import \
           --input-path  "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/800_12S_single_end_ee3-seq_q2taxtable.tsv" \
           --output-path "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/800_12S_single_end_ee3-seq_q2taxtable.qza" \
           --type 'FeatureData[Taxonomy]' \
           --input-format HeaderlessTSVTaxonomyFormat || { echo 'Taxonomy import failed' ; exit 1; }`
  * Qiime taxonomy file with Blast taxonomy available at `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/800_12S_single_end_ee3-seq_q2taxtable.qza`    
  * adjusted and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/900_q2_summary.sh`
    * see for length statistics of raw data `qiime tools view /Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/900_12S_single_end_ee3-seq.qzv`
    * see for metadata-wise statistics of raw data `qiime tools view /Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/900_12S_single_end_ee3-tab.qzv`
  * to export objects, adjusted and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh`
    * exported to `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/980_12S_single_end_ee3-tab_q2_export`
    * correct duplicate column names in `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/980_12S_single_end_ee3-tab_q2_export/taxonomy.tsv`
      * to `Feature ID	Taxon	Confidence`
  * commit `bd8978a72804f8c03dc4c5420f9492d2af0efc6f`

## **12-Dec-2020** - started to work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

    * **next**:
      * sort columns of imported data to match columns
      * possibly chase why columns are not properly sorted

## **08-Jan-2021** - continued work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

  * backtracking pipeline
  * updated `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh` to match
    * `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/850_prep_q2_predictor-tab__metadata.tsv` as written by
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/850_r_prep_q2_predictor-tab.r`
  * checking `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/900_12S_single_end_ee3-tab.qzv`
    * metadata headers are ok as generated by:
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/900_q2_summary.sh`
  * with updated headers (see above) attempting to re-run `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh`
    * compressing `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/980_12S_single_end_ee3-tab_q2_export.zip`
    * erasing contents of `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/980_12S_single_end_ee3-tab_q2_export`  
    * re-running `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh` - **ok**
  * re-running `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
    * understood, commit `64f61205e54334e9282616b0fa221c5362aa338c`
  * import completed - next adjust Phyloseq object processing
    * in `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
    * akin to `/Users/paul/Documents/OU_pcm_eukaryotes/Github/200_r_get_phyloseq.r`
    * commit `91f947c9878c81582b6f4da3e41c318d6426412b`
    * initial plot done, keep on stepping through later ...
    * commit `d674b4a5027346f69233362f6cb9957f63034a51`

## **11-Jan-2021** - continued work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

  * started coding text summary
  * commit `e9c426b970f49a5c45c9ba81e6d464f832a06e9d`  

## **20-Jan-2021** - continued work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

  * can't get package decontam to work
    * can't work around error `conc must be positive numeric`, stays in despite `conc` being positive numeric
    * commit before not using package `c5f23c1a8b14fc4a703b3a40f2834adcee2fddd5`
    * commit with purely subtractive filetring draft `cfe0e80d9592dc03fc872f6691f8f0f89993460b`

## **21-Jan-2021** - continued work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

  * restored `990_r_get_eDNA_phyloseq.r` from commit `c5f23c1a8b14fc4a703b3a40f2834adcee2fddd5`
  * implemented new filtering strategy - filtering is ok, script may need to be structured again for better readability
  * commit `d935013849472b7cc72761648819a64cb86625e`

## **22-Jan-2021** - continued work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

  * below line `531` subtracting all controls instead of only positive controls, not subtracting blanks
  * commit `90e090a85eae12f05eb68acaedfaeff38e57b806`
  * added ASV number plotting function - improved plotting speed by aggregation in both plotting functions
  * cleaned structure
  * started to add function to display vars in plate format
  * commit ` 0b7c76e3fddab69d644c90ec7af3cdda29dc4792`

## **22-Jan-2021** - continued work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

  * finished to add function to display vars in plate format - still slow for unsubset data could be sped up by aggregation
  * commit `faa88147a078e1c2e3738a573d605a83fcbb3275`

## **26-Jan-2021** - continued work on `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`

  * isolated barcodes in which positive controls amples were found
  * commit ` 31a87cc92cef16b6793b89684f8d496e2ead2816`

## **27-Jan-2021** - continued work on checking primer assignments

  * checked barcode labels
    * using: `/Users/paul/Documents/OU_eDNA/191031_primers/200504_idt_template_plate.xls`
    * using: `210126__primers_to_spot_check__200_r_get_phyloseq.xlsx`
  * checked well assignments using lab book
    * using  `200616_idt_plate_filled_fjordland_primers`
  * assuming logically erroneous demultiplexing or primer contamination
    * commit before new demultiplexing `ca5be2f24b46218d54fa4c0c2759a46024b7154a`
  * created backup copy in `/Users/paul/Documents/OU_eDNA_backup`
  * restarting from script `/Users/paul/Documents/OU_eDNA/200901_scripts/200_r_metadata_management.R`
    * finished re-run and increased saving of intermediary objects to `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata`
    * commit `be9756b77616c8024a777ceb4266260b6de41157`
  * re-writing `/Users/paul/Documents/OU_eDNA/200901_scripts/300_bash_cutadapt_demultiplex.sh` 
    * demultiplexing via temp file, subtracting reads from file, rather then using raw data several times
    * running draft version on local
    * commit `3abf0ace62d3d77d6da90265131616f2785dc934`

## **28-Jan-2021** - re-running pipeline

  * got working version of, and running: `/Users/paul/Documents/OU_eDNA/200901_scripts/300_bash_cutadapt_demultiplex.sh`
  * commit `dc9e4a3ebc2dd0662202c512e8e711e26916258d`
  * needed to re-test several length filtering options
  * running latest version again after several tries, script was edited extensively
  * commit: `1bfe5018350cc09308cbe2d192188f8476003a82`

## **01-Feb-2021** - re-run finished of `/Users/paul/Documents/OU_eDNA/200901_scripts/300_bash_cutadapt_demultiplex.sh`

  * used only about half the data at high quality filtering
  * the rest is left here (1.68 of originally 3.8 GB): `/Users/paul/Documents/OU_eDNA/201126_preprocessing/cutadapt/300_bash_cutadapt_demultiplex_input.fastq.gz `
  * commit `cb50223b82534652a0f38d312640d6c94c3228b2`
  * created and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/350_bash_count_reads_and_mv_empty_fastqs.sh`
  * commit `af41faf4482f8f863e2ba300ac060ffc5684601a`
  * adjusted and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/400_r_qiime_manifest.R`
  * commit `d278de85b9180bc19c6b6f0b2b5690a0ba265c8b`
  * uploading files to Cornell cluster

## **02-Feb-2021** - on cluster using `qiime2-2020.6`

  * running `./500_q2_import.sh` - **ok**
  * running `./600_q2_denoise.sh` - **ok**
    * `The filter removed all reads: /tmp/tmpdr8mtydj/U-1-H-1-empty_139_L001_R1_001.fastq.gz not written.`
    * `The filter removed all reads: /tmp/tmpdr8mtydj/U-2-E-5-ncntrl-pcr_160_L001_R1_001.fastq.gz written.`
    * `766815368 total bases in 4196353 reads from 162 samples will be used for learning the error rates`
  * manually exporting ASV set for blasting: contains 2171 sequences as per `grep ">" 600_12S_single_end_ee3-seq.fasta | wc -l`
  * starting Blasting of ASV set with script `750_bash_fasta_blast.sh` - **ok**
  * afterwards compressed environmental GI list again - arrived on local- **ok**
  * in MEGAN check texonomy see `/Users/paul/Documents/OU_eDNA/201126_preprocessing/megan/750_12S_single_end_ee3-seq_blast-noenv.xml.rma6`
    * `/Users/paul/Documents/OU_eDNA/201126_preprocessing/blast/750_12S_single_end_ee3-seq_blast-noenv.xml.gz`
    * `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-seq.fasta.gz`
  * run `/Users/paul/Documents/OU_eDNA/200901_scripts/650_r_plot_denoise.R`- **ok**
    * saved `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_si_auxilliary_files/210202_650_denoised_libraries.pdf`
  * started running `/Users/paul/Documents/OU_eDNA/200901_scripts/800_r_get_q2_tax-tab.r` - **pending**
    * creating `blast_results_list` - **pending**
    * after import erase `/Users/paul/Documents/OU_eDNA/201126_preprocessing/blast/750_2S_single_end_ee3-seq_blast-noenv.xml` - **pending**
    * committing ` 64a63bd8f7564eddad24ea14155323dd35135d19`
  * done running `/Users/paul/Documents/OU_eDNA/200901_scripts/800_r_get_q2_tax-tab.r`
    * commit `32d01ef1f47ca0711c90dd0dcf7f6a18156d5b99`

## **03-Feb-2021** - continuing second processing iteration

  * stepping through `/Users/paul/Documents/OU_eDNA/200901_scripts/850_r_prep_q2_predictor-tab.r`
  * running  `/Users/paul/Documents/OU_eDNA/200901_scripts/900_q2_summary.sh`
    * using `qiime2-2020.8` - **ok**
    * importing taxonomy file:
      `qiime tools import \
         --input-path  "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/800_12S_single_end_ee3-seq_q2taxtable.tsv" \
         --output-path "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/800_12S_single_end_ee3-seq_q2taxtable.qza" \
         --type 'FeatureData[Taxonomy]' \
         --input-format HeaderlessTSVTaxonomyFormat || { echo 'Taxonomy import failed' ; exit 1; }`
  * running `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh`  - **ok**
  * stepping through `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
    * looking much better now - work on publishable version next
    * commit `d6e15bd6e25c18e91f56540bf8a0d43bcee7b7f2`

## **10-Feb-2021** - continuing second processing iteration

  * re-working processing script `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
    * got to line 505 - prior to contamination removal
    * commit `1335c18362eb60828bb52d28def5bfbea281b394`

## **11-Feb-2021** - continued second processing iteration

  * finished processing script `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
  * commit `6d74175af74d98c4069995294201e4e06e00d564`

## **12-Feb-2021** - continued second processing iteration

  * continued processing script `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
    * combined sequences counts across biological replicates - **ok**
    * next:
      * erase old files, update names, and re-write
      * tweak length filtering
      * indicate quality of alignments in clean data object
  * commit `1381ff1ff616e0f3cb7b59240e7e075a7520dd26`

## **19-Feb-2021** - eDNA data check 

  * implemented and checked low abundance filtering in `~/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
    * based on Poisson distribution of cross-contaminated positive control read counts
  * commit `dd61faf0ead4e372ae5c92b3699421ac97505ed2`

## **22-Feb-2021** - eDNA data work and subsequent

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_phyloseq.r`
    * added Blast information to clean eDNA data 
    * updated all output files
  * commit `abb0ec2bf56e40dee0318250c337963da0baeb24`  

## **24-Feb-2021** - eDNA data work and subsequent

  * started working on `~/Documents/OU_eDNA/200901_scripts/995_r_get_BRUV_phyloseq.r`
    * looks up taxonomy strings from local data base by means of manually matched tax ids
    * writes random strings for NCBI tax ids to 
      * `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210224_995_r_get_BRUV_phyloseq_mh bruv_taxa_db_lookup_stable.xlsx`
  * **next**: 
    * check sample locations and coordinate format and manually match up before reading into above script:
      * `/Users/paul/Documents/OU_eDNA/191213_field_work/220224_MH_bruv_data_machine_readable_preformat_long.csv`
      * `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210211_990_r_get_eDNA_phyloseq__eDNA_sampling_success.xlsx`
    * taxonomy id key is:
      * `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210224_995_r_get_BRUV_phyloseq_mh bruv_taxa_db_lookup_stable.xlsx`

## **25-Feb-2021** - BRUV observations received from MH

  * created `/Users/paul/Documents/OU_eDNA/191213_field_work/210225_MH_bruv_data_machine_readable_long.csv`
    * using `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210211_990_r_get_eDNA_phyloseq__eDNA_sampling_success.Rds`
    * using `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210211_990_r_get_eDNA_phyloseq__eDNA_sampling_success.xlsx`
  * looking up tax information from DB in script `/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_BRUV_phyloseq.r`

## **25-Feb-2021** - merging BRUV and eDNA data

  * ready to merge BRUV with eDNA object and send off
    * finished `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_long_table.r`
    * finished `/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_BRUV_long_table.r`
  * commit `a2646b91ba1f1d13050a35c49a64090b803c013e`

## **25-Feb-2021** - continue merging BRUV and eDNA data

  * starte `~/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r` 
  * next:
    * check for proper joining in: `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_long_table.r`
    * continue working on: `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
  * commit `a2646b91ba1f1d13050a35c49a64090b803c013e`

## **01-Mar-2021** - finished formatting file data

  * stored at:
    * `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/210301_997_r_format_longtables__analysis_input.Rds`
    * `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210301_997_r_format_longtables__analysis_input.xlsx`
  * commit `d80b20cc53f5942fb322383523e157a8d26e3911`

## **07-Mar-2021** - working on maps

  * create map using `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210307_sample_map.qgz`
  * exporting data files for analysis and mapping from `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
    * add export for mapping: `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/210301_997_r_format_longtables__analysis_input.csv`

## **09-Mar-2021** - working on results summaries

  * created `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`, commit `df7577353c5baef344b0f502d75b5d07d2625e0b`
    * to summarise long table exported by `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
    * to summarize M.d.L.'s derivative of above file `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210309_mdl_tablebygenus.csv`
    * to summarize M.d.L.'s derivative of above file `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210309_mdl_tablebyspecies.csv`

## **11-Mar-2021** - working on results summaries

  * altered `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
    * renamed the cryptic A, B, C locations in RESERVE.GROUP with something more meaningful (WJ, FF, LS)
    * defined a new variable RESERVE.GROUP.LOCATION 
  * implemented MDS and basic plot in: `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * commit `81d0142db89dfed924a4097efedc761b56354cad`

## **12-Mar-2021** - working on results recode

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * implemented BRIUV vs eDNA comparison - **drafted**
    * implemented NMDS on Jaccards - **started**
    * commit `d3f766a72163f746394092c583c4cad437ea21fd`
    * improvef NMDS plot to GGplot
    * commit `1bc5379ca2d37ed7d042f9d95026fc07b639a7d1`

## **13-Mar-2021** - worked on results recode

    * implements indicator species analysis  - as outlined
    * implement MCS - across SET.ID to complement smallish plot
    * commit `d95a49974e3f70a387c1b2b9f468832864429238` 

## **16-Mar-2021** - worked on results recode

  * completed ANOSIM  analysis - as outlined
  * modified figure labels
  * added asterisk ("*") for non-NZ fish
* **next**
    * check if ICC usage is appropriate
  * commit ` 9c2b10d79c3bd28e76527ef6038ddbe34ee36145`

## **23-Mar-2021** - finding ICC alternative

  * committing unknown changes - `4e4a1207350954a180597a2223eb324cd6daeb16`
  * implemented Venn diagrams and re-arranged figures
  * commit `1e77045d533b7e64c7afa00fa14d7eec747286f1`

## **24-Mar-2021** - started NMDS overhaul - all points, no grouping

  * committed beforehand: `e668968a27c3e9a5179a2a57855ba0906f47feed`
  * saved code as scratch `/Users/paul/Documents/OU_eDNA/201126_script_scratch/998_r_summarize_results_more_nmds.r`
  * and restored above commit

## **25-Mar-2021** - code checks for doc revision

  *  commit `74f97bba1a22bf54b8ee939ee65b713d6378ec4`

## **29-Mar-2021** - started to work on better map

  * using `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`

## **31-Mar-2021** - continued to work on better map


## **01-Apr-2021** - continued to work on better map

  * commit `3f5d2ad4a60600cf5850cfaecd42158da595431f`

## **07-Apr-2021** - working on MCA replacement

  * removed genus map from code in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * removed depreciated ICC code in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * encoded new barplot and new figure arrangement in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * commit `017aedd5bdb61db66bc97057ce259eb33e8896c2`

## **08-Apr-2021** - working on MCA replacement

  * revised NMFDS code
  * started on numerical summaries
  * commit `bce2874a613beeaf8ec134d5f10a03c78df8bf95`

## **09-Apr-2021** - working on heat-map as component of Fig. 2

  * created heat-map
  * added margin sums in the best possible way
  * ***next** 
    * re-export plot with heat map
    * finish numerical summaries using margin calculations in heat-map code
    * expand text

## **10-Apr-2021** - working on figure assembly and numerical summary 

  * working on `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * temp blue whale image credits: `By NOAA Photo Library - anim1754, Public Domain, https://commons.wikimedia.org/w/index.php?curid=17942391`
  * temp dolphin image credits `By NASA - http://mediaarchive.ksc.nasa.gov/detail.cfm?mediaid=21807, Public Domain, https://commons.wikimedia.org/w/index.php?curid=112006`
  * git commit `92e74fb25446c3be931402de3b2a29fa36c6caaa`

## **21-Apr-2021** - working on manuscript revision

  * in silico PCR and plots
  * advancing numerical summaries
  * modified
    * `000_r_in_silico_pcr.R`
    * `990_r_get_eDNA_long_table.r`
    * `998_r_summarize_results.r`
  * commit `75a1b1d669632ee10677954f06e8e285f0ea5c51`

## **29-Apr-2021** - working on manuscript revision: final analysis code again, also separately for eDNA and BRUV

  * adjusted and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
    * output saved to `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components`
    * report at `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/code_reports`
  * adjusted and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_eDNA.r`
    * output saved to `/Users/paul/Documents/OU_eDNA/200403_manuscript/7_si_auxillary_files/*eDNAonly.pdf`
    * report at `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/code_reports`
  * adjusted and ran `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_BRUV.r`
    * output saved to `/Users/paul/Documents/OU_eDNA/200403_manuscript/7_si_auxillary_files/*BRUVonly.pdf`
    * report at `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/code_reports`
  * commit ` fbeb26f5dd58d7f1df59f284a6b9dcb83425a278`
    * revising manuscript files using report
      * `/Users/paul/Documents/OU_eDNA/200403_manuscript/210430_210421_main_text_changes_accepted_WR_ML_mod-while-withMK.docx`
      * `/Users/paul/Documents/OU_eDNA/200403_manuscript/210430_key_numbers_by_dataset_for_SI.xlsx`
      * `/Users/paul/Documents/OU_eDNA/200403_manuscript/210430_si.lyx`

## **04-May-2021** - working on manuscript revision: final analysis code again, also separately for eDNA and BRUV

  * finished re-running:
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_BRUV.r`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_eDNA.r`
  * updated SI
  * commit `f7a40a484b0ec6a120484373e326367cc3d41cd`

## **13-May-2021** - working on manuscript revision: **include OBIS data**

  * started trails with `~/Documents/OU_eDNA/200901_scripts/001_fetch_format_obis.R`
  * steps
    * fetch OBIS data
    * adjust pipeline design and analysis
    * re-create imagery
    * re-outline text
  * pre-alteration commit `0e9513088fa9d2257b41c936ddd1067777b0a610`

## **14-May-2021** - working on manuscript revision: **include OBIS data**

 * use `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_get_OBIS_long_table.R`
   * as drop-in replacement for `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`

## **24-May-2021** - working on manuscript revision - preparing manuscript for Ecography

  * for higher impacte or if rejected include OBIS data at a later stage
  * moving file to scratch scripts:  `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_get_OBIS_long_table.R`

## **05-July-2021** - working on manuscript revision - adding OBIS data to analysis

  * edited, but keeping old functionality in `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
  * started on `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r`
    * fetches OBIS data
    * **next** - re-implement mapping from `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_BRUV.r`
    * **next** - save OBIS data and thin out and add to long table
  * edited `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_BRUV.r`
    * **next** - remove code duplicated in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r`
      * some table formatting code
      * mapping code
  * commit `91b4453d2c00881bd15d8eb7ab691a48f387daba`
  * continued  `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r`
    * re-implemented mapping from `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_BRUV.r`
    * save OBIS data 
      * **next** thin out data and add to long table - with taxon lookup
      * **next** deal with missing data / reselect complete data?
    * edited `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results_BRUV.r`
    * **next** - remove code duplicated in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r`
      * some table formatting code
      * mapping code
  * commit `151ea60de34334b5cc7f30d302732230dac200f6`
  * continued  `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r`
  * commit `edb67644f7e8e7b881831cfad6a5960458e104b0`

## **07-July-2021** - working on manuscript revision - adding OBIS data to analysis

  * checked `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
  * worked on `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r`
    * OBIS data citations saved to `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210707_OBIS_data_citations.xlsx`
    * assessed data completeness
    * for subsequent analyses saves object `long_table` to
      * `/Users/paul/Documents/OU_eDNA/201028_Robjects/998_r_map_and_add_obiss__full_data_raw.Rds`
      * `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/R_objects/998_r_map_and_add_obiss__full_data_raw.Rds`
  * starting to work on `~/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * commit before (extensive) edits (for earlier versions check earlier commits before 5-July-2021)
  * commit `bb6a1342d022e5098e62655fd1ea14ad9a113b3c`
  * wrote functions to get Euler plots
  * commit `24481cee30c0d9dcd8cebc2ad6a7f971e41b0b1b`
  * got Euler plots
  * commit `262df352764d7f6dc70d59aff5b815d380db1b15`
  * started on bar plots

## **08-July-2021** - working on manuscript revision - coding of new display items

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * tried plot overlays
  * commit `222da8e9862c695aa5c07bea5f87bc79c03b1ee1`

## **08-July-2021** - working on manuscript revision - coding of new display items

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * finished Euler plot after data filtering
    * finished heat map after data filtering
  * commit `0443735b32522ce6912ebebd0f137f6750ca9f04`
  * commit `f0575dae1f5d7dbafa7a3347b397412cea2b08ac`

## **15-July-2021** - working on manuscript revision - of tree-sorted fish list

  * getting NCBI tax ids for all fish so as to be able to use NCBI tree via ETE toolkit (`http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html`)
  * **done** - revised `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r` - eDNA ASVs have NCBI tax ids again.
  * **done** - revise `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r` - to keep values in column `NCBI.TAXID`
  * **pending** get Newick Trees from NCBI ids in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r` via ETE toolkit and Phyton
  * trying to get to work `/Users/paul/Documents/OU_eDNA/200901_scripts/get_tree_for_ncbi_taxid_vector.py`
    * **NCBI data installation running but pending and needed to be set up properly**

## **16-July-2021** - working on manuscript revision - of tree-sorted fish list

  * continued working on `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * **pending** - tree needs to be gotten after database works
    * **pending** - Anosim results need to be obtained **next** - subset to observation type via  function?
    * **pending** - indicator species analysis ?
  * commit `77415cc05c999016401f491d16fc33eb7aa77ce6`

## **19-July-2021** - working on manuscript revision - of tree-sorted fish list

  * continued working on `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * **pending** - tree needs to be gotten after database works 
    * **debug filter** - ANOSIM results need to be obtained
    * **pending** - indicator species analysis
  * commit `f0498f34161835d1f21792862d5f6182a5125a16`

## **20-July-2021** - working on manuscript revision

  * continued working on `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * now looking up trivial names
    * now putting out taxonomically ordered flextable of full biodiversity
    * now putting out taxonomically ordered heat-map of full biodiversity
  * commit `0778b9b979d94347cc5fc217404c9505fb6b0758`
    * printing table with full data
  * commit `f19d6ced14469b458158ae40156b47b454159a6f`
    * testeding ANOSIM function - should be bug free
  * commit `86939ab7513314d2ba5b593265fe3ccd2d594a8f `

## **22-July-2021** - adding species observations from literature and OBIS with a wide-circle buffer

  * started working on `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_include_literature_data.r`

## **23-July-2021** - adding species observations from literature and OBIS with a wide-circle buffer

  * checking `/Users/paul/Documents/OU_eDNA/200901_scripts/800_r_get_q2_tax-tab.r`
  * re-running `/Users/paul/Documents/OU_eDNA/200901_scripts/990_r_get_eDNA_long_table.r`
    * defining `NCBI.TAXDB.INC`
    * checking query coverages - **do again later**: `select(PHYLUM, CLASS, ORDER , FAMILY, GENUS, SPECIES, HSP.IDENTITY.PERC) %>% print(n = Inf)`
  * re-running and checking `/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_BRUV_long_table.r`
  * finishing `~/Documents/OU_eDNA/200901_scripts/997_r_include_literature_data.r`
      * defining `NCBI.TAXDB.INC`
  * finished `~/Documents/OU_eDNA/200901_scripts/995_r_get_PUBL_long_table.r`
  * renamed other script files as per commit message
  * next stepping through `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
  * commit `91f11ad4db6c0d6121201dd1594f6ebb5753c5f5`
  * re-ran `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r`
  * `PUBL` data added, `NCBI.TAXID` set to `numeric()`.
  * commit `6b6fb3f766148d2618fd33dde37e36b6d1d340e6`
  * starting to step through `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_get_OBIS_and_map.r`
    * add OBIS data with wide circle
    * correct / fill / append missing variables
    * wrote new map to `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components`
    * saved workspace file to `/Users/paul/Documents/OU_eDNA/210705_r_workspaces/210723_998_r_get_OBIS_and_map.Rdata`
    * commit `a43812d9e2ac15cbd4f0e98222b0c1c384776b37`
    * **next** continue in line `173`

## **24-July-2021** - adding species observations from literature and OBIS with a wide-circle buffer

  * revising `~/Documents/OU_eDNA/200901_scripts/995_r_get_PUBL_long_table.r` - adding grouping variables - **ok**
  * revising `/Users/paul/Documents/OU_eDNA/200901_scripts/997_r_format_longtables.r` - checking correct execution - **ok**
  * revising `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_get_OBIS_and_map.r`
    * progressed to downloading OBIS data, line 205
    * commit `3d3ea8c09e00af1b3c8bdfa9aabbade4102039c1`

## **26-July-2021** - adding species observations from literature and OBIS with a wide-circle buffer

  * revising `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_get_OBIS_and_map.r`
    * added OBIS data with set ID 99
    * formatted all other data (hopefully to completion)
  * commit `9fdde5710e0d1feb6e6695b77df8313c3e6bef78`
  * to correct errors in taxonomy - revised all scripts again
    * from `/Users/paul/Documents/OU_eDNA/200901_scripts/995_r_get_PUBL_long_table.r`
    * to `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * ***next** - continue with analysis in script above
  * commit `48894396fd709c0796ab3a148d082adb540b7986`
  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * finished Euler plots
    * commit `6dda994ad29940d4af7b91128b17a8e49ef577a6`

## **27-July-2021** - getting new display items

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * before extensive revisions of heat map plotting code, commit `3777c1961324b15aa5f901d60c0c6c10bb40c2cb`
    * everything done - display items saved - but ANOVA stuff, commit `4dd9307c14a8e6723c0f21eaf0725e945b6f8336`

## **28-July-2021** - started on working abstract

  * added various data summary code section to `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * commit `6a0276177db5bfc0a87a68af2819989e73e10302`

## **29-July-2021** - implemented ANOSIM

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * implemented ANOSIM
  * commit `15e8f8d8870c607cdbd6b91e4c6b40c8fa4a44e6`
    * implemented indicator species analysis
    * erased older code
    * saved object fro MdL
  * commit `b7c90e2de647ca75f13c3ca5bb9d0ed6e9eb37e4` before proper save
  * commit `9ffba95df4ff2ae183c793172cddc59d33cb1880`

## **06-Aug-2021** - writing paper

  * added barplot code script from MdL
  * added big circle centerpoint code to `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_get_OBIS_and_map.r`
  * commit `d9e1e2b79c3620267d3a7d68414a239510279e6a`

## **09-Aug-2021** - preparing talk

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * added code to generte individual display items
  * no commit yet

## **12-Aug-2021** - adding brief ASV overlap analysis

  * working on `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
  * commit prior to further changes: `3d742087e425120f60c7dc5c03235f2dea364a52`
  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * started adding code code for `ASV.PER.LOC ~ HSP.GAPS + HSP.IDENTITY.PERCENT + NOT.NZ`
    * emailed object to MdL

## **13-Aug-2021** - analysed alignment parameters against species overlap

  * for `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * had other branch open - but not needed
    * continue code at line `773` to complete ASV summaries
  * commit `31496fed427e8c75494cc1ad983cdc6264b1045c`

## **14-Aug-2021** - finished main text and mailed off, started on SI

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * improved model plotting code
    * included code for results section of main text
  * commit `5ea60d60c94cc46c869cc812c2f94337949ac980`

## **22-Aug-2021** - starting revisions for submission

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * lot's of changed based on chat with MdL
    * removing regression code based on averages
    * rescaling coverage fractions to percentage
    * recalculating confidence intervals

## **23-Aug-2021** - starting revisions for submission

  * in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * implemented binomial regression and many other things
  * **commit** before erasing Poisson regression
  * adjusted regressions
  * ran through code completely
    * table formatting is different - results of other computer?
    * p-values are slightly different in ANOSIM - due to permutation testing
    * copying files with DIs from macmini
  * **commit** at code end with checked script

## **26-Aug-2021** - continuing manuscript revisions

  * "finished" main text
  * file `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_summarize_results.r`
    * likely in final state - but needs to be re-run on macmini for correct display items 

## **27-Aug-2021** - implementing additions requested by MK

  * checking `~/Documents/OU_eDNA/200901_scripts/998_r_get_OBIS_and_map.r`

## **10-Sep-2021** - collating data for online upload

  * use `find . -type f -newermt '05/13/2021 0:00:00'`
    * include folders `201028_Robjects`
    * include folders `210705_r_workspaces`
    * include folders `200901_scripts`
    * include folders `201126_preprocessing`
  * find slahes in script files - get a list of files that possibly need archiving
    * `find . -type f -exec grep -H '/' {} +* >> /Users/paul/Documents/OU_eDNA/200403_manuscript/210910_file_names.txt`
  * see script `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/210910_collate_files.sh`
  * commit `b3c4f7b703fb8ff9c43d68c8b4109f94b45ae7b0`
  * Zenodo upload and temporary file cleaning
    * compressing for Zenodo upload: `tar -czvf 210910_SI_for_review.tar.gz /Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/collated_files`
    * erasing duplicated files that can be re-collated: `rm -rfd /Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/collated_files/*`
    * after upload to Zenodo: `rm /Users/paul/Documents/OU_eDNA/200403_manuscript/210910_SI_for_review.tar.gz`

## **26-Apr-2022** - preparing re-submission to **ACS Environmental Science & Technology**

  * as per guidelines `/Users/paul/Documents/OU_eDNA/200403_manuscript/220425_ACS-EST_guidelines.pdf`
    * create a TOC graphic

## **29-Apr-2022** - preparing re-submission to **ACS Environmental Science & Technology**

  * finishing manuscript, upload to ACS, and distribution of submission proof
  * soft-linking `/Users/paul/Documents/OU_eDNA/200403_manuscript/7_si_auxillary_files` to 
    * `/Users/paul/Documents/OU_eDNA/200403_manuscript/220219_CONL-21-0339_rejection`
    * `/Users/paul/Documents/OU_eDNA/200403_manuscript/220407_ACS/7_si_auxillary_files`

## **20-Sep-2022** - preparing second submission to **Environmental DNA**

  * slightly re-sorting files and updating links
  * notes on first rejection are can be found at `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/220916_response_text.docx`
  * see `/Users/paul/Documents/OU_eDNA/200403_manuscript/README.md` for more detailed information

## **27-Oct-2022** - preparing second submission to **Environmental DNA**

  * starting re-analysis on MacBook
  * full back up of old project is on MacMini
  * include new Blast results
  * include new fish table from Roberts et al. 2019
  * include Meta-Fish lib
  * commit prior to re-analysis is `e418d2bac877465426dc475f11617e8fe17fa5cb`
  * adjusting and executing `/OU_eDNA/200901_scripts/000_r_in_silico_pcr.R` - ok
  * commit `defc7ac2021d9542a4e2585cad3ca0ac475dd13b`

## **29-Oct-2022** - preparing second submission to **Environmental DNA**

  * checking and running - with partial re-saves
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/100_r_quantification_analysis.R`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/200_r_metadata_management.R`
  * checking  - without re-saves
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/300_bash_cutadapt_demultiplex.sh`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/350_bash_count_reads_and_mv_empty_fastqs.sh`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/400_r_qiime_manifest.R`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/500_q2_import.sh`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/600_q2_denoise.sh`
  * commit `48466a1dc8d8fe71098f4a32258296fd30874ba9`
  * checking and running - with partial re-saves
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/650_r_plot_denoise.R`
    * notes on `/Users/paul/Documents/OU_eDNA/200901_scripts/750_bash_fasta_blast.sh`
      * code is now depreciated:
      * during revision possibly attempting re-blast
      * new reference data has since become available via NCBI as per most recent reviewers literature highlights
      * Blast should now be done on NESI
      * negative GI list is probably outdated
      * MetaFishLib is available as well and should be used
  * **finished** to adjust transport scripts in `/Users/paul/Documents/OU_eDNA/221027_nesi_transport`
    * commit is `c6ab245ee3ede3c8ef6c819af76197422ba9df96`
  * **starting** to work on `/Users/paul/Documents/OU_eDNA/200901_scripts/751_bash_fasta_blast_nesi.sh`
  * see `https://support.nesi.org.nz/hc/en-gb/articles/208619807-BLAST`
  * downloading new negative GI list
    * 20008447 entries
    * link was `https://www.ncbi.nlm.nih.gov/nuccore/?term=%22environmental%20samples%22[organism]%20OR%20metagenomes[orgn]`
    * new file will be `/Users/paul/Documents/OU_eDNA/201126_preprocessing/blast/221027_gi_list_environmental.txt`

## **30-Oct-2022** - preparing second submission to **Environmental DNA**

  * finished `/Users/paul/Documents/OU_eDNA/200901_scripts/751_bash_fasta_blast_nesi.sh`
  * creating custom BLAST database with sequences created by MetaFishLibNZ 
  * in `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ`
  * see `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/README.md`
  * adjusted script names - one for full `nt` db, one for `221030_MetaFishLibNZ`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/751_bash_fasta_blast_fullnt.sh`
    * `/Users/paul/Documents/OU_eDNA/200901_scripts/752_bash_fasta_blast_metafishlib.sh`
  * running on local: `/Users/paul/Documents/OU_eDNA/200901_scripts/752_bash_fasta_blast_metafishlib.sh`
  * checked results in Megan - no taxonomy visible
  * to try added files below downloaded 31-Oct-2022 - doesn't work  
    * `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ/taxdb.btd`
    * `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ/taxdb.bti`
  * also added `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ/nucl_gb.accession2taxid.gz`
  * if usable MEAGAN will use NCBI taxonomy for MetaFishLib, not Fishbase taxonomy, needs to be appended afterwards
  * `makeblastdb` needs to be repeated - see `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/README.md`  

## **31-Oct-2022** Building Blast reference data base of MetaFishLibNZ sequences **with NCBI taxonomy**

  * **fish base taxonomy needs to be added later to both datasets**
  * for details see `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/README.md`
  * **unsuccessful** - try other approach - learn how to use MEGAN and building reference databases for BLAST

## **2-Nov-2022** Attempting first BLAST on NESI

  * will likley fail

## **4-Nov-2022** Attempting second BLAST on NESI

  * may not fail if script is otherwise error-free and "$TMPDIR" gets set by job managementq
  * failed due to incorrectly set `BLASTDB` variable - re-attempting

## **7-Nov-2022** Attempting third  BLAST on NESI 

 * edited `751_bash_fasta_blast_fullnt.sh` attempt to re-run

## **11-Nov-2022** - BLAST run failed

 * I do not understand the error message, and can't find informetion on it, filing issue
 * upping memory, updating dates, committing, re-attempt

## **12-Nov-2022** - BLAST run failed again

  * pulling files to MacMini
  * commit `d93005ac6e090b70140c73e8b012805aeb58e4b1`

## **14-Nov-2022** - re-attempting BLAST on cluster
 
 * checked syntax in `/Users/paul/Documents/OU_eDNA/200901_scripts/751_bash_fasta_blast_fullnt.sh`
 * bumped up requested BLAST version
 * commit and push to cluster

## **15-Nov-2022** - BLAST done on cluster

 * files pulled to local
 * created `/Users/paul/Documents/OU_eDNA/200901_scripts/801_r_get_q2_tax-tab.r`
   * to parse new BLAST results
   * to calculate Bit-score cutoff as per (Liu et al. 2010)
 * installing package `taxonomizr` in R
 * installing `taxonomizr` database in `/Users/paul/Sequences/References/taxonomizr`
   * probably done despite error, for further setup see `https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html`

## **17-Nov-2022** - experimenting with Bit-Score cutoff

  * started in `/Users/paul/Documents/OU_eDNA/200901_scripts/801_r_get_q2_tax-tab.r`
  * next calculate cutoff values
  * commit `193efc06aa60a897b34b33a7fc420c5ab9cfab4f`

## **18-Nov-2022** - experimenting with Bit-Score cutoff

  * continued in `/Users/paul/Documents/OU_eDNA/200901_scripts/801_r_get_q2_tax-tab.r`
  * plotted and exported bit score assignments in comparison to the average values
  * commit `52078808b2dea9f4d35506cc2dd5f255a5a7b8fb`

## **23-Nov-2022** - defined Bit-Score cutoff variables, starting taxonomic annotation from NCBI
 
 * continued in `/Users/paul/Documents/OU_eDNA/200901_scripts/801_r_get_q2_tax-tab.r`
 * see change history, added variables to indicate above or below-average max. hsp bit score
 * re-installing `taxonomizR` data base
   * see `https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html`
   * R: `httr::timeout(5 * 3600) # increase curl's maximum request time`
   * R: `taxonomizr::prepareDatabase("accessionTaxa.sql", tmpDir =  "/Users/paul/Sequences/References/taxonomizR_221115")`
 * commit `74218f36d1e8085ec687fc5281437454fe6fef63`  

## **24-Nov-2022** - defined Bit-Score cutoff variables, starting taxonomic annotation from NCBI
 
 * `prepareDatabase()` fails repeatedly
 * checking `https://github.com/sherrillmix/taxonomizr/issues/33`
 * downloading accession to taxid manually:
   * `~/Sequences/References/taxonomizR_221115$ wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz` 
 * preparing sqlite file:
   * `read.accession2taxid("/Users/paul/Sequences/References/taxonomizR_221115/nucl_gb.accession2taxid.gz","/Users/paul/Sequences/References/taxonomizR_221115/nucl_gb.accession2taxid.sqlite", vocal=TRUE)`
 * continued in `/Users/paul/Documents/OU_eDNA/200901_scripts/801_r_get_q2_tax-tab.r`
   * stopped at line 189 - `Part II: Format taxonomy table for export and export`
   * to continue load `/Users/paul/Documents/OU_eDNA/200901_scripts/ou_edna_workspace.Rdata`
 * commit `0665dab8e626e3133b5a41a623a7262842b6559f`

## **25-Nov-2022** - stepping through script

 * stepped through `/Users/paul/Documents/OU_eDNA/200901_scripts/801_r_get_q2_tax-tab.r`
 * inspected exported `/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/801_12S_single_end_ee3-seq_q2taxtable.tsv`
 * `open -a "Microsoft Excel" /Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/801_12S_single_end_ee3-seq_q2taxtable.tsv`
 * commit `00180b3563a0055d8c103482360891e60288dfa9`

## **09-Dec-2022** - continuing analysis refresh for manuscript resubmission

 * stepped through `850_r_prep_q2_predictor-tab.r`
 * running  `/Users/paul/Documents/OU_eDNA/200901_scripts/900_q2_summary.sh`
    * using `qiime2-2021.4`, not as previously `qiime2-2020.8` - **ok**
    * converting taxonomy file with new Blast results from `tsv` to Qiime as done **03-Feb-2021* :
      `qiime tools import \
         --input-path  "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/801_12S_single_end_ee3-seq_q2taxtable.tsv" \
         --output-path "/Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/801_12S_single_end_ee3-seq_q2taxtable.qza" \
         --type 'FeatureData[Taxonomy]' \
         --input-format HeaderlessTSVTaxonomyFormat || { echo 'Taxonomy import failed' ; exit 1; }`
 
 * running `/Users/paul/Documents/OU_eDNA/200901_scripts/980_q2_export_objects.sh`  - **ok**
 * inspected `900_q2_summary.sh` - skipping
 * inspected `980_q2_export_objects.sh` - skipping
 * creating copy for editing `901_q2_summary.sh` - **ran successfully**
 * creating copy for editing `981_q2_export_objects` - **ran successfully**
 * stepping through copied `991_r_get_eDNA_long_table.r`
   * completed until line 455
 * commit `e80c09e3857aa624e2379e61572b6d8e86b412fe` - line 455
 * commit `50f3ea87db24fc46b45fffa0064213cb13c57335` - line 533
 
## **10-Dec-2022** - continuing analysis refresh for manuscript resubmission

 * finished stepping through `991_r_get_eDNA_long_table.r`
 * finished stepping through `996_r_get_BRUV_long_table.r`
 * commit `2885bbc4ee470765621d97850aeb4f3cf0d79f34`

## **14-Dec-2022** - continuing analysis refresh for manuscript resubmission

 * start stepping through `996_r_get_PUBL_long_table.r`
 * attempt including synonyms lookup via `rfishbase`
 * one fish has been renamed
   * in corrected list as "Forsterygion capito"
   * in literature list as "Grahamina capito"
 * NCBI tax lookup succede for "Hypoplectrodes huntii", found 14.12.2022
 * formerly: 34 Genera, 42 species fully resolved
 * today:    36 Genera, 44 species fully resolved
 * finished stepping through `996_r_get_PUBL_long_table.r`
 * commit `0f0aeb4fa649cf333a6d05b20e0075f72f0306fe`
 * finished walking through copied file `998_r_format_longtables`
 * created copy `999_r_get_OBIS_and_map.r`
 * created copy `999_r_summarize_results.r`
 * commit `de39e0d376d3e3fb1b91b60448e63e95d03357d6`

## **19-Dec-2022** - continuing analysis refresh for manuscript resubmission

 * start stepping through `999_r_get_OBIS_and_map.r` and check for synonyms
   * corrected export file names in `998_r_format_longtables.r`
  * used date "221214"
  * re-exported map with new crs `2193` (NZGD2000 / New Zealand Transverse Mercator 2000 -- New Zealand Transverse Mercator (NZTM))
  * new map at `/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/221219_999_r_summarize_results_fig1_draft.pdf`
  * stared OBIS lookup 
  * commit `9961020d7fec12cf09031b1140792a56ac247cd9`
  * finished OBIS and NCBI lookup, saved workspace
  * commit `feaf2689071ec98186670951370115d53202e8cf` 

## **20-Dec-2022** - continuing analysis refresh for manuscript resubmission

 * finished stepping through `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_get_OBIS_and_map.r`
   * commit `ed5c516ee15276b23f4c146176e66c3da303f55a`
 * step through `999_r_summarize_results.r` and check for synonyms
  * commit `46bfef8c02565406c20e5ec2aece31d0079cf3ca`

## **23-Jan-2023** - continuing analysis refresh for manuscript resubmission

 * continued stepping through `999_r_summarize_results.r`
   * corrected all synonyms that could be corrected
   * started looking up trivial names manually
   * added some notes to new main text copy (with todays date)
   * commit `abe9dcfbee3e25365eda0d4ce0c7ae0d94f30bd9`

## **06-Mar-2023** - continuing analysis refresh for manuscript resubmission

 * finished up trivial names in `999_r_summarize_results.r`
 
## **16-Mar-2023** - continuing analysis refresh for manuscript resubmission
 
 * finished formatting trivial names in `999_r_summarize_results.r`
   * ran `save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__got_more_trivial-names.Rdata")`
   * added readable outline to script for easier navigation
   * commit `foo`
   * for manual lookup writing Excel table writing `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part.xlsx`
   * for manual lookup taking notes in `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part_annotated.xlsx`
   * using `/Users/paul/Documents/OU_eDNA/200224_references/210908_MA_DOC001887_TePapa_Checklist-of-Fishes-of_full.pdf`
   * updated README and script further
   * commit `8282a315a7dd315e92673c0b1fdd6159c5070abd`

## **12-May-2023** - continuing analysis refresh for manuscript resubmission

 * loading `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
 * loading most recent environment `/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__got_more_trivial-names.Rdata`
 * finished NZ native status and references in `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part_annotated.xlsx` 
 * continue stepping through `999_r_summarize_results.r`
    * in line 655
    * sample type variable isn't set correctly 
    * should be `SAMPLE.TYPE %in% c("eDNA","BRUV","OBIS")` not 1,2,3,4
    * bug is in line 531 of `999_r_get_OBIS_and_map.r`
 * commit `b37436ed9d4a6381ed0a1c891f29dc1db1c0a140`
 * worked off more of todu queue - see below
 * commit `4b9c74b483ec828fd60eb3ff27d96b6c8c696f6a`

## **15-May-2023** - continuing analysis refresh for manuscript resubmission

 * running script `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
 * from line 0 to 591 with saving all intermediate results - look for current data if needed
 * committing prior to modyfying from line 858 regressed mapping code
 * commit `9e59780468df78287660d56b9fd4350def4d195e`
 * finished running `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
 * commit `13db500b665c7e42ba71b8569ca63667f540de18`

## **22-May-2023** - continuing analysis refresh for manuscript resubmission

 * started updating SI text and DI paths
 * moving DI's to dedicated folder
 * revised folder structure for manuscript submission

## **24-May-2023** - restarting manuscript work

 * finished updated SI
 * started revising and appending rebuttal letter
 * minor changes in script during collation of data for SI
 * commit `eff46f29c79b9984e86833e484b28c6f59c033c9`

## **31-May-2023** - restarting manuscript work

 * updated SI to contain BLAST regression
 * updated title page
 * copied new raw DIS from `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components`
 * for adjustment to `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_main_text_di_development`
 * erased as errors where detected, only keeping alias
 * finished revising abstract

## **06-June-2023** - resuming manuscript work

 * stepped through rebuttal - to query 35
 * in MEAGAN created 
   * `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_new_analysis_outputs/751_12S_single_end_ee3-seq_blast-noenv-ex__read_to_taxon.txt`
   * `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_new_analysis_outputs/751_12S_single_end_ee3-seq_blast-noenv-ex__taxon_to_count.txt`
 * copying to `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_new_analysis_outputs`
   * added BLAST results table to `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230606_si.lyx`

## **09-June-2023** - resuming manuscript work
 
 * in `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
   * continuing in line 475
   * loading `/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__got_more_trivial-names.Rdata`
   * continue with: 
     * `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__started_megan_integration.Rdata")`
     * adding MEGAN annotation to table
     * commenting on set differences and intersection in manuscript 
 * commit `a70a08c1edebbbb870aacf259734f2eb979df73d`

## **14-June-2023** - re-running `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
 
 * added MEGAN annotations to export items
 * exported raw DI's for manual adjustment
 * progressed to line 1315 and saved environment - continue with section XI line 1317
 * commit `2e53436d971ccae8100c5666fc0c9370883d68d8`

## **16-June-2023** - finished re-running `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`

 * see `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
 * for paths to new display items
 * for new numerical summaries
 * copying `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__geoheat_edna_bruv_obis.pdf`
 * to `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/figure_1.pdf`
 * updated figure 1 in folder and manuscript and results script
 * updated figure 2 as well in main file and manuscript `cp /Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__biodiv_tiles_only.pdf /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/figure_2.pdf`
 * updated figure 3 as well in main file and manuscript `cp /Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__euler_edna_bruv_obis.pdf /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/figure_3_unedited.pdf` 
 * updated table 1 in main text: `cp /Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__spcies_obs_matching_tiles.html /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/table_1.html` 
 * updated table 1 in main text: `cp /Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__spcies_obs_matching_tiles.docx /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/table_1.docx` 
 * overwrote last file after successful lyx import - did all work transforming table 1 in lyx 

## **19-June-2023** - finished re-running `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
  
  * in manuscript text refreshing results reporting
  * thus adding many reporting lines to script `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
  * commit `e21ce694d8f4e332c1a6bd4012db41595247067e`
  * continued updating numerical summaries in manuscript 
  * thus adding even more reporting lines to script `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
  * commit `30220f570e5d839676d276f10cd19f236fcc18aa` 

## **20-June-2023** - working on main text
 
 * updated results section - extended results reporting script
 * commit `69a55986856d2834ed1ed106fa3ea495fbed8d37`

## **23-June-2023** - restructured and re-ran confusing script `999_r_summarize_results.r`
 
 * script `999_r_summarize_results.r` is messy: 
   * restructuring script to be comprehensible
   * done to section `Get Euler plots`
   * saved environment 
   * commit `706aa8b4fe810329ac1d3372075dc2b923f1ea86`
   * updated todo list and code structure
   * commit `ca1bfe8605562ceaed581876364cfaae81589c09`
 * script `999_r_summarize_results.r` is messy: 
   * re-running from section `2e`, 
   * saving all states after cleaning script 
   * checking concordance of results with previous results -looks ok results are consistent
 * re-ran script `999_r_summarize_results.r` from   section `2e` to finish saving all intermediate environments
 * commit `8c7c288ff7f4c8e8268034ddd611c63507666e32`
 * updated ToDo below 
 * commit `d4cfcdeda2aa006e181d8076051463f2314f8a52`


## **16-August-2023** - revising rebuttal and main text
 
 * updated cover letter and cover page
 * updated rebuttal letter
 * updated main text
 * updated important queue
 * commit `0ef9ed513dd6b73c5c4e765b6c5fde27e7559712`
 
## **17-August-2023** - revising  main text

 * working on main text - extending outline document
 
## **18-August-2023** - revising  main text

 * working on main text  - extending outline document

## **21-August-2023** - work day 127 - revising  main text

 * working on main text  - extending results section
 * getting new branch `git checkout -b "results_revision""`
 * erasing superseded scripts `rm 998_r_summarize_results.r 800_r_get_q2_tax-tab.r 990_r_get_eDNA_long_table.r 995_r_get_BRUV_long_table.r 995_r_get_PUBL_long_table.r 997_r_format_longtables.r 998_r_get_OBIS_and_map.r`
 * opening script `999_r_summarize_results.r`
  * loading data from section 4g line 608 `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__finished_megan_integration.Rdata")`
  * re-running to add information to outline
  * started to match up script strucuture to results, marks set in text file and script 
  * -> continue at main text, Results 1, point 2.
  * -> continue at script line 925 ("Augmentation of  OBIS and literature by eDNA and BRUV")
  * commit `b14eb83a123d91984da3b9d7904f907f4956b827` 

## **23-August-2023** - work day 128 - revising  main text, results

 * finished results section `[0]` and `[1]`
 * matching `999_r_summarize_results.r` to results section `[0]` and `[1]`
 * loading data from section 4g line 608 `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__finished_megan_integration.Rdata")`
 * running to subsequent save point - next load from line 874, just prior to `Results for main text`
 * `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping_done.Rdata")`
 * stepped through to line 1013, and updated results
 * commit `9f79855c0f6d8a14c976055128e3ab2fc33588a3`

## **24-August-2023** - work day 128 - revising  main text, results

 * in script `999_r_summarize_results.r`
   * re-ran from `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping_done.Rdata")`
   * saved state at `save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__moidelling_done.Rdata")`
   * saved state at `save.image(file = "/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__summaries_done.Rdata")`
   * saved state at `save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__dis_done.Rdata")`
   * saved state after ANOSIM `save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__anaosim_done.Rdata")`
   * saved state after indicator species analysis `save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__multipatt_done.Rdata")`
 * obtaing stable results from `999_r_summarize_results.r`
 * in word file with todays date (24.08.2023) results are done matching file
 * commit `5ab354b0177127910ea23bcc97fb196fea8bdf0e`

## **25-August-2023** - work day 129 - revising  main text, discussion
 
 * no notes 

## **28-August-2023** - work day 130 - revising all for submission to coauthors
 
 * created main text version for coauthors
 * checked and rendered SI
 * copied current version for coauthors to online staorage lacally availabel at `/Users/paul/Library/CloudStorage/OneDrive-UniversityofOtago/OU_Fiordland_eDNA`
 * commit `29b8533e827d419d5484f497d79cbc9649bfd770`

## **31-August-2023** - work day 131 - sending modelling results to MdK

 * in script `999_r_summarize_results.r` 
   * loading `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping_done.Rdata")`
   * sending regression input object to MdL `/Users/paul/Documents/OU_eDNA/201028_Robjects/230514_999_r_summarize_results__data_spc_distribution_vs_quality.Rds`

## **01-September-2023** - work day 131

 * [x] revising `230905_cover_letter.lyx`
 * [x] revising `230901_response_text.docx`
   * [x] after download from online repository
 * [x] revising `230901_main_text.docx` with web version and removing scaffolding
   * [x] after download from online repository
   * [x] adress Wills commnets
   * [x] adress MdLs comments
   * [x] update MdLs regression results - to main text and SI 
     * [x] `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230816_communication/230901_logistic_regression_fish.pdf`
     * [x] `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230816_communication/230901_logistic_regression_fish.Rmd`
     * [x] integrated MDL's results screenshot:
       * [x] `mv "/Users/paul/Desktop/220422_screenshots/Screenshot 2023-09-01 at 15.42.10.png" "/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_si_di_development/7_model_coefficients_v3.png"`
       * [x] `mv "/Users/paul/Desktop/220422_screenshots/Screenshot 2023-09-01 at 16.07.12.png" "/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_si_di_development/7_model_coefficients_v3_b.png"`
 * [x] mail off to MK with invite to chat about submission 

## **06-September-2023** - work day 132

 * revised regression results in online version as per MdL

## **08-September-2023** - work day 133

 * online - revised rebuttal letter edits made by MK

## **11-September-2023** - work day 134

 * revised all materials
 * uploaded all materials to shared drive
 * started submission process
 * pdf collation calls were:

```
502  2023-09-11 17:55:25 gs -sDEVICE=pdfwrite  -dNOPAUSE -dBATCH -dSAFER -sOutputFile=/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230911_main_text_supplement_for_review.pdf /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230911_main_text.pdf /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230911_si.pdf
503  2023-09-11 18:03:28 gs -sDEVICE=pdfwrite  -dNOPAUSE -dBATCH -dSAFER -sOutputFile=/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230911_cover_letter__response_text.pdf /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230911_cover_letter.pdf  /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230911_response_text.pdf
```
 * provided proof to coauthors

## **12-September-2023** - work day 134

 * received comments on proofs updated todo
 * commit `f838c4d19625f1fc77a3daa755d685086598b103`

## **19-September-2023** - work day 135

 * saved notes on proof correction at `/220826_eDNA_resubmission/230816_communication/230919_on_proof_correction.pdf`
 * updated todo - worked on todo - see correction of proofs
 * adjusting map in main text (Fig. 1) obfuscating sampling locations
   * in `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
   * loading from line 636 `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/998_r_summarize_results__data_filtered.Rdata")`
   * last modification Aug 23 13:33:18 2023 (`stat`)
   * in `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_get_OBIS_and_map.r`
   * loading from line 581 `/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221220_999_r_map_and_add_obis__end.Rdata`
   * last modification "May 12 17:42:50 2023 (`stat`)
   * modifying map generation in line 157
     * minimal changes only, only created
     * only created `999_r_get_OBIS_and_map__mapggplot_redacted.Rds"`
     * in new branch `new_maps` commit `7c78d08b9e9b6278341dc4143f00ab79e539a250`
   * adjusting in `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r`
     * from line 748
     * resaving at `save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping.Rdata")`
     * resaving at `save.image("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__mapping_done.Rdata")`
   * see new files at (circle size in subfigures needs to be adjusted manually - check commit history of main summary script)
     * `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results_map_main.pdf`
     * `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results_map_edna.pdf`
     * `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results_map_bruv.pdf`
     * `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results_map_obis.pdf`
     * `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components/230515_999_r_summarize_results__geoheat_edna_bruv_obis.pdf`
 * new fig 1 at `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/figure_1.pdf`
 * commit `2bea89e5c3bf9bdc899e575a2edaa33719cd322d`
 * cleaning out git repository as per `https://gist.github.com/boywijnmaalen/f052cd70fa6d924d58621b2e9c806adb`

## **19-September-2023** - work day 135

  * moving files to Zenodo and Zenodo folder (`/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230920_Zenodo_contents`)
    * code release 1.2. from GitHub (with own DOI on Zenodo: `10.5281/zenodo.8359890`) `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230920_Zenodo_contents/Fiordland-eDNA-v1.2.zip`
    * MEGAN file: `cp /Users/paul/Documents/OU_eDNA/201126_preprocessing/megan/751_12S_single_end_ee3-seq_blast-noenv.rma6 /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230920_Zenodo_contents/230920_megan.rma6`
    * BLAST results: `cp /Users/paul/Documents/OU_eDNA/201126_preprocessing/blast/751_12S_single_end_ee3-seq_blast-noenv.xml.gz /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230920_Zenodo_contents/230920_blast.xml.gz`
    * Qiime-exported-reads: `cp /Users/paul/Documents/OU_eDNA/201126_preprocessing/qiime/600_12S_single_end_ee3-seq.fasta.gz /Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230920_Zenodo_contents/230920_asv_sequences.fasta.gz`
  * commit

## **13-Nov-2023** - work day 136 - received comments for second revision

  * received comments for second revision
    * see `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/231113_eDNA_revisions`
    * also updated `/Users/paul/Documents/OU_eDNA/200403_manuscript/README.md`

## **16-Nov-2023** - work day 137 - starting on second revision

 * copying `/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230919_main_text.docx`
   * for keeping to `/OU_eDNA/200403_manuscript/9_submissions/231113_eDNA_revisions/230919_main_text.docx`
   * for editing to `/OU_eDNA/200403_manuscript/9_submissions/231113_eDNA_revisions/231116_main_text.docx` 
 * starting to address revisions in `/OU_eDNA/200403_manuscript/9_submissions/231113_eDNA_revisions/231116_response_text.docx`

## **17-Nov-2023** - work day 138

 * continued on second revision

## **13-Dec-2023** - work day 138

 * continued on second revision
 * checked script `998_r_format_longtables.r`
   * from loaded workspace `/Users/paul/Documents/OU_eDNA/210705_r_workspaces/221214_998_r_format_longtables__analysis_input__image.Rdata`
 * working in `/Users/paul/Documents/OU_eDNA/200901_scripts/991_r_get_eDNA_long_table.r`
   * checked and adjusted summary code
   * no commit (yet)

## **18-Dec-2023** - work day 139

 * finished second revision
 * uploaded to manuscript Central
 * see `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/231113_eDNA_revisions`
 * commit

## **Pending important queue** (last updated 19-September-2023):
 
 * [x] finish revisions in folder `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/231113_eDNA_revisions`

## **Done from important queue** (last updated 16-August-2023):

* [x] revise erroneous species tables 
    * [x] by finishing correcting and re-running `/Users/paul/Documents/OU_eDNA/200901_scripts/999_r_summarize_results.r` 
    * [x] by adding MEGAN results from `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/230522_new_analysis_outputs/751_12S_single_end_ee3-seq_blast-noenv-ex__read_to_taxon.txt`
    * [x] by revising duplicates in species table DI and exports
 * [x] get new main display items in `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components`
 * [x] replace old display items in `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission`
 * [X] update percentages in results text, with main Fig. 2
 * [X] decide whether to work on main text or supplemnet next
 * [x] update results, beyond line 206
 * [x] update results, beyond line 210 
 * [x] finish re-structuring results script `999_r_summarize_results.r` 
 * [x] re-run `999_r_summarize_results.r`  to get stable environmental saves
 * [x] establish stable analysis results available from  `999_r_summarize_results.r` , otherwise update
 * [x] update cover letter
 * [x] main text - finish results section
 * [x] finish and re-run from  stable version of `999_r_summarize_results.r` from `load("/Users/paul/Documents/OU_eDNA/210705_r_workspaces/999_r_summarize_results__finished_megan_integration.Rdata")`
 * [x] get notes document - assemble new aspects for improvements
 * [x] main text - explode for restructuring using reviewer queries
 * [x] main text - adjust Introduction
 * [x] main text - adjust Results, possibly while extending `999_r_summarize_results.r`
 * [x] main text - adjust Discussion
 * [x] main text -re-read entirely and check for flow
 * [x] supplement - adjust after main text 
 * [x] rebuttal - adjust after supplement
 * [x] circulate to co-authors
 * [x] all files - check crosslinks
 * [x] all files - revise with co-authors suggestions - **use online version** 
 * [x] correct proofs 
    * [x] change back title
    * [x] remove map from SI
    * [x] remove specific locations from main text map
    * [x] remove specific locations from main text
    * [x] update data availability statement
    * [x] upload files to Manuscript Central
 * [x] re-upload files to Zenodo
   * [x] all files - collate data - obscure location data
   * [x] all files - upload data - request removal of old versions
 * [x] submit

## **Miscellaneous todo (last updated 16-Mar-2023):**

 * for **new manuscript**
 * figures for **supplement** are here (see dates): `/Users/paul/Documents/OU_eDNA/200403_manuscript/9_submissions/220826_eDNA_resubmission/221118_analysis_outputs`
 * figures for **main text** are here (see dates): `/Users/paul/Documents/OU_eDNA/200403_manuscript/3_main_figures_and_tables_components`
 * **for supplement**: keep in mind fish and non-NZ-fish, and references in `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part_annotated.xlsx`
 * **for rebuttal**: keep in mind changed overal counts in script 999...: line ~643 starting with `# Summary: general species counts`
 * keep in mind that ASV counts counts have changed - see script 999
 * keep in mind that Species counts counts have changed - see script 999
 * keep in mind that new plots have been generated e.g. by script 990
 * run finished MetaFishLib Blast`/Users/paul/Documents/OU_eDNA/200901_scripts/752_bash_fasta_blast_nesi.sh`
 * get taxonomy information for MetaFishLib
   * email Ruppert - how to do this? - MEGAN / Blast
 * new files have been saved for images and SI, starting with 27-Oct-2022
 * include those new files into manuscript and main text
 * in supplement cross-reference `/Users/paul/Documents/OU_eDNA/200403_manuscript/5_online_repository/tables/210707_OBIS_data_citations.xlsx`
 * **add comment regarding OBIS data completeness** calculated in `/Users/paul/Documents/OU_eDNA/200901_scripts/998_r_map_and_add_obis.r`
   *  "after data cleaning and assignment of NCBI taxonomy 826 of 1020 OBIS records were retained fro analysis (81%)"
 * delete outdated 
  * `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ/taxdb.btd`
  * `/Users/paul/Documents/OU_eDNA/220928_meta_fish_lib/221030_MetaFishLibNZ/taxdb.bti`
 * provide notes on **noteworthy fish species** in manuscript and their completeness** 
  * found **Fiordichthys slartibartfasti**
  * found **Eldon's galaxias**
 * provide notes on **non-fish species** in manuscript
  * see **file**: `/Users/paul/Documents/OU_eDNA/200403_manuscript/6_analysis_notes/999_r_summarize_results__long_table__part_annotated.xlsx`
  * found **Thalassarche chrysostoma** - **Grey-headed albatross**
  * found **Thalassarche melanophris** - **Black-browed albatross**
  * found **Megaptera novaeangliae** - **humpback whale**                 
  * found **Diomedea exulans** - **Wandering albatross**
  * found **Arctocephalus forsteri** - **New Zealand fur seal**
  * found **Arctocephalus australis** - **South American fur seal**
  * found **Puffinus griseus** - **Sooty shearwater**           
  * found **Aptenodytes patagonicus** - **king penguin**
  * found **Synoicum kuranui** - possibly medicinal tunicate
  * found **Trididemnum shawi** rare tunicate
   
## Notes on depreciated scripts and why they are depreciated

* `/Users/paul/Documents/OU_eDNA/200901_scripts/750_bash_fasta_blast.sh`
  * used Cornell cluster
  * re-implemented and run on 14-Nov-2022 in `/Users/paul/Documents/OU_eDNA/200901_scripts/751_bash_fasta_blast_fullnt.sh`
* `/Users/paul/Documents/OU_eDNA/200901_scripts/752_bash_fasta_blast_metafishlib.sh`
  * unused so far - can be adjusted if MetafishLib can be used for Blasting, if needed
* `/Users/paul/Documents/OU_eDNA/200901_scripts/800_r_get_q2_tax-tab.r`
  * used until 21-Jun-2022 for "eDNA" submission
  * superseded for "eDNA"  resubmission after 26-Aug-2022 by `/Users/paul/Documents/OU_eDNA/200901_scripts/801_r_get_q2_tax-tab.r`
