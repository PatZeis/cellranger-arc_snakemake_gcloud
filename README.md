# Pipeline for single-cell multiomic RNA and ATAC libraries based on cellranger-arc

## Features:
- Demultiplexing of RNA and ATAC bcl files 
- Alignment to reference genome of multiomic reads
- Transcript and fragment counting 

## Software requirements:
- Google Cloud SDK, connection to user account & gcloud service account with credentials set in shell session 
- Snakemake according to https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html

## cellranger-arc requirements  
- download cellranger-arc-2.0.2.tar.gz and store at gcloud storage in dedicated workflow folder
cellranger-arc needs to be installed separately as it can not be build using conda. Therefore donwload cellranger-arc installation tar archive 'cellranger-arc-2.0.2.tar.gz' from 'https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest' and store it a gcloud storage in the folowing path: <insertbucket>/data/mkfastq/. Example path to 'cellranger-arc-2.0.2.tar.gz' at gcloud storage: zeis_workflows/atlas_V1/data/mkfastq/cellranger-arc-2.0.2.tar.gz 
- download or build reference genome for cellranger-arc and store at gcloud storage in dedicated worklow folder
compressed and archived reference genome folder for Human or mouse can be downloaded from 'https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest' or build according to: 'https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/release-notes/references#GRCh38-2020-A-2.0.0'. After genome build or download, store compressed and archived genome folder at gcloud in the following folder: <insertbucket>/data/ref/. Example path to ref genome: zeis_workflows/atlas_V1/data/ref/ref_genome_101_Macaca.gz

## prerequisits to run the pipeline
- archive bcl folders and upload to gcloud
Archive bcl folder of RNA and ATAC: 
``` bash    
tar -cvf bcl_folder_rna.tar bcl_folder_rna
tar -cvf bcl_folder_atac.tar bcl_folder_atac 
```
Store archived bcl folder at gcloud in the dedicated folder: <insertbucket>/data/mkfastq/ for RNA or <insertbucket>/data/mkfastq/ATAC/ for ATAC. Example paths to RNA and ATAC archived bcl folders: zeis_workflows/atlas_V1/data/mkfastq/bcl_folder.tar and zeis_workflows/atlas_V1/data/mkfastq/ATAC/bcl_atac.tar
- generate sample reference file and upload to gcloud storage
After cloning the github repository, change to directory: cellranger-arc_snakemake_gcloud/config/ and edit the 'generate_mksample_file.sh' script with sample names, corresponding RNA and ATAC sample index sets( samples,indexes should have matching order separated by ".") as well as flow cell names for RNA and ATAC flow cell runs (can be derived from RunInfo.xml file).
After editing the script generate the 'mkfastq_samples.tsv' 
``` bash 
./generate_mksample_file.sh
``` 
and upload to gcloud storage into the folder <insertbucket>/data/ref/. Example path to 'mkfastq_samples.tsv' zeis_workflows/atlas_V1/data/ref/mkfastq_samples.tsv'
- generate config.yaml
edit the bucket in the 'generate_config.sh' script, then execute the script:
``` bash 
./generate_config.sh
```   

## run snakemake

``` bash
snakemake --google-lifesciences --default-remote-prefix <bucket_name> --use-conda --google-lifesciences-region <region> -j <num_jobs>
```

