# Fiordldand project utility scripts

* **1-Sep-2020**
  * creating `/Users/paul/Documents/OU_eDNA/200901_scripts/200901_quantification_calibration.R`
    * associating fluorescence signal with DNA concentrations, lab book page 48
    * using file `/Users/paul/Documents/OU_eDNA/200128_lab_work/200901_qpcr_qbit_test_results_formatted_data_with_qbit.xlsx`
    * created file `/Users/paul/Documents/OU_eDNA/200128_lab_work/200901_qpcr_qbit_test_results_regressions.pdf`
  * commit `e02097aac4f7cb65d800f267d2e7b7cf397aff7`
* **16-Sep-2020**
  * creating `/Users/paul/Documents/OU_eDNA/200901_scripts/200916_quantification_analysis.R`
* **17-Sep-2020**
  * some unneeded edits on `/Users/paul/Documents/OU_eDNA/200901_scripts/200916_quantification_analysis.R`
  * needs to be updated for new standards
  * using HR reads, as they seem to have performed better:
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate1_HS_ng-ul.xlsx`
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate2_HS_ng-ul.xlsx`
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate1_HS_ng-ul.xlsx`
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate2_HS_ng-ul.xlsx`
  * commit `e277b04a54b5599ad3ce94ad0a7be375335e211d`
* **28-Oct-2020** - creating `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`
  * writing objects to `/Users/paul/Documents/OU_eDNA/201028_Robjects`
  * Read cells from
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200907_plate_layouts.xlsx`
    * `/Users/paul/Documents/OU_eDNA/191031_primers/200302_Bunce_et_al_0000_MiFishEmod_single_step_primers.xlsx`
  * output appropriate format for
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/201028_Otago_Genomics_sample_info.xlsx`
    * and later possibly the mapping file
  * next: spread and subset for Otago Genomics - **pending**
  * commit `c51356eb489093cc333cd4c5024d3d535ce362a2`
* **29-Oct-2020** - continuing `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`
  * finished compiling sample data for Otago Genomics
  * commit `64b9cf9f69114a0b3f6d038368e7a6cc5bcb64aa`
* **26-Nov-2020** - getting metadata for sequence deconvolution
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
* **30-Nov-2020** - getting metadata for sequence deconvolution
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
* **01-Dec-2020** - renaming files - getting metadata for sequence deconvolution
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
* **03-Dec-2020** - creating manifest file for Qiime imports
  * demultiplexing has finished `/Users/paul/Documents/OU_eDNA/200901_scripts/300_conda_cutadapt_demultiplex.sh` - has finished
    * files at: `/Users/paul/Documents/OU_eDNA/201126_preprocessing/cutadapt`
  * updated `/Users/paul/Documents/OU_eDNA/200901_scripts/200_sample_matadata_managment.R`
    * for manifest file generation exports all data as R objcets, notably also
      * all metadata combined: `/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__big_table.Rdata`
      * object used for demultiplexing information `/Users/paul/Documents/OU_eDNA/201028_Robjects/201028_sample_managment__demux_table.Rdata`
        * which is R object corresponding to `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/200_cutadapt_barcode_input.txt`
  * created script to create manifest file for Qiime: `/Users/paul/Documents/OU_eDNA/200901_scripts/400_create_qiime_manifest.R`
    * manifest at: `/Users/paul/Documents/OU_eDNA/201126_preprocessing/metadata/400_create_qiime_manifest__manifest.txt`
  * 
        
    
    

* **unfinished**
  * use `/Users/paul/Documents/OU_eDNA/200901_scripts/200_sample_matadata_managment.R` to import file - **pending**
  * use `/Users/paul/Documents/OU_eDNA/200901_scripts/200_sample_matadata_managment.R` to mapping file - **pending**
  * `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`
    * finish for creation of mapping file - **pending**
  * read in new BRUV observations as alternatives to `/Users/paul/Documents/OU_eDNA/191213_field_work/200520_MH_bruv_data.xlsx`

