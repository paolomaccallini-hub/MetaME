# file name: MetaME_func
#
library(httr) 
library(R.utils) 
library(data.table)
library(MungeSumstats)
library(yaml) 
library(ggplot2) 
library(patchwork)
#
#-------------------------------------------------------------------------------
# Read configurations from YAML file
#-------------------------------------------------------------------------------
#
config<-read_yaml("MetaME_config.yml")
#
MAFco_uc<-as.numeric(config$filters$maf_uncommon) # lower cut-off for minor allele frequency of uncommon variants
INFOco<-as.numeric(config$filters$info_cutoff) # imputation quality cut-off
#
samples<-config$samples
#
metal_path<-config$METAL$path_metal_exe
#
populations<-c()
#
h<-0
for (i in 1:length(samples)) {
  if (samples[[i]]==1) {
    h<-h+1
    populations[h]<-names(samples)[i]
  }
}
#
#-------------------------------------------------------------------------------
# Add a data folder, if absent
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
#-------------------------------------------------------------------------------
# Build data base for DecodeME
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/DecodeME")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Download summary statistics, if not present
#
url<-"https://osf.io/download/v4w8g/"
file_path<-file.path(current_dir,"Data/DecodeME/gwas_1.regenie.gz")
if(!file.exists(file_path)) {
  print("Downloading DecodeME summary statistics")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
}   
#
# Download filtered variants, if not present
#
url<-"https://osf.io/download/6uj5x/"
file_path<-file.path(current_dir,"Data/DecodeME/gwas_qced.var.gz")
if(!file.exists(file_path)) {
  print("Downloading DecodeME filtered variants")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
}
#
# Download ReadME, if not present
#
url<-"https://osf.io/download/axp4k/"
file_path<-file.path(current_dir,"Data/DecodeME/readme.txt")
if(!file.exists(file_path)) {
  print("Downloading DecodeME readME")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
}
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Build data base for Million Veteran Project (MVP)
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/MVP")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Download summary statistics, if not present
#
url<-"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479178/GCST90479178.tsv.gz"
file_path<-file.path(current_dir,"Data/MVP/GCST90479178.tsv.gz")
if(!file.exists(file_path)) {
  print("Downloading Million Veteran Project summary statistics")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
}   
#
# Download meta data, if not present
#
url<-"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479178/GCST90479178.tsv.gz-meta.yaml"
file_path<-file.path(current_dir,"Data/MVP/GCST90479178.tsv.gz-meta.yaml")
if(!file.exists(file_path)) {
  print("Downloading MVP readME")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
}  
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Build database for UK Biobank (European Institute of Bioinformatics) 
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/EIB")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Download summary statistics, if not present
#
url<-"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038694/GCST90038694_buildGRCh37.tsv"
file_path<-file.path(current_dir,"Data/EIB/GCST90038694.tsv")
file_path_gz<-file.path(current_dir,"Data/EIB/GCST90038694.tsv.gz")
if(!file.exists(file_path_gz)) {
  print("Downloading UK Biobank EIB summary statistics")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
  gzip(file_path,destname=file_path_gz,remove=TRUE) # zip it and remove unzipped version
}  
#
# Download meta data, if not present
#
url<-"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038694/GCST90038694_buildGRCh37.tsv-meta.yaml"
file_path<-file.path(current_dir,"Data/EIB/GCST90038694.tsv-meta.yaml")
if(!file.exists(file_path)) {
  print("Downloading EIB readME")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
}  
gc() # free unused memory
