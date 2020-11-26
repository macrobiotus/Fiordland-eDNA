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
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate1_ng-ul.xlsx`
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Eplate2_ng-ul.xlsx
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate1_ng-ul.xlsx
    * `/Users/paul/Documents/OU_eDNA/200128_lab_work/200916_Uplate2_ng-ul.xlsx`
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
  * copied for adjusting `/Users/paul/Documents/OU_eDNA/201126_preprocessing/scripts/100_q2_import.sh`
    * next - check in Qiime documentation and forum - to develop demultiplexing strategy
      * reverse complementing
      * dual barcodes single end reads
  * commit `b51522c7bb238cf8883e793ff0b565917a1016a` in `/Users/paul/Documents/OU_eDNA/201126_preprocessing/script`
        

* **unfinished**
  * `/Users/paul/Documents/OU_eDNA/200901_scripts/201028_sample_matadata_managment.R`
    * use for creation of mapping file - **pending**
    * read in library concentrations befor pooling for mapping file  - **pending**
    * read in other sample metadata - **pending**

