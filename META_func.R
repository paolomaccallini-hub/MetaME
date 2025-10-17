# file name: META_func
#
library(httr) 
library(R.utils) # for .bgz files
library(data.table)
library(gtexr) # GTEx 
library(otargen) # Open Target
library(MungeSumstats)
library(yaml) # for settings reading
library(corrplot) # for corrplot
#
#-------------------------------------------------------------------------------
# Read configurations from YAML file
#-------------------------------------------------------------------------------
#
config<-read_yaml("META_config.yml")
#
pco_c<-as.numeric(config$filters$p_value_common) # p value cut-off for common variants
pco_uc<-as.numeric(config$filters$p_value_uncommon) # p value cut-off for uncommon variants
pco_HWE<-as.numeric(config$filters$hwe_p_value) # p value for HW equilibrium
MAFco_c<-as.numeric(config$filters$maf_common) # lower cut-off for minor allele frequency of common variants
MAFco_uc<-as.numeric(config$filters$maf_uncommon) # lower cut-off for minor allele frequency of uncommon variants
INFOco<-as.numeric(config$filters$info_cutoff) # imputation quality cut-off
#
samples<-config$samples
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
# Build database for UK Biobank (Neale Lab)
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/NealeLab")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Download phenotype data, version 2, if not present
#
for (sex in c("both_sexes")) {
  url<-paste0("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.",sex,".v2.tsv.bgz")
  destfile<-paste0("phenotypes.",sex,".v2.tsv.bgz")
  file_path<-file.path(current_dir,"Data/NealeLab/",destfile)
  if(!file.exists(file_path)) {
    print("Downloading the list of phenotypes...")
    RETRY(
      verb = "GET",
      url = url,
      write_disk(file_path, overwrite = TRUE),
      times = 5,           # up to 5 attempts
      pause_min = 5,       # wait 5s between attempts
      terminate_on = c(404,403) # don't retry on these errors
    )
  }   
}
#
# Download complete list of imputed variants, if not present
#
url<-"https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
destfile<-"variants.tsv.bgz"
file_path<-file.path(current_dir,"Data/NealeLab",destfile)
if(!file.exists(file_path)) {
  print("Downloading and editing the full set of imputed variants from Neale Lab (UK Biobank)")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  ) 
  #
  # Edit the file keeping only the columns: myvariants and p_hwe (this takes a while)
  #
  all_variants<-fread(file_path,header="auto",sep="\t")
  all_variants<-all_variants[,c("variant","rsid","p_hwe","info")]
  write.table(all_variants,file=file_path,sep="\t",col.names=T,row.names=F)
}  
#
# Download summary statistics for the selected phenotypes
#
for (phenotype in c("20002_1482")) {
  for (sex in c("both_sexes")) {
    url<-paste0("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/",phenotype,".gwas.imputed_v3.",sex,".tsv.bgz")
    destfile<-paste0(phenotype,".gwas.imputed_v3.",sex,".tsv.bgz")
    file_path<-file.path(current_dir,"Data/NealeLab/",destfile)
    if(!file.exists(file_path)) {
      print("Downloading summary statistics")
      RETRY(
        verb = "GET",
        url = url,
        write_disk(file_path, overwrite = TRUE),
        times = 5,           # up to 5 attempts
        pause_min = 5,       # wait 5s between attempts
        terminate_on = c(404, 403) # don't retry on these errors
      )
    }  
  }
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
file_path<-file.path(current_dir,"Data/EIB/GCST90038694.tsv.gz")
if(!file.exists(file_path)) {
  print("Downloading UK Biobank EIB summary statistics")
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
url<-"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038694/GCST90038694_buildGRCh37.tsv-meta.yaml"
file_path<-file.path(current_dir,"Data/EIB/GCST90038694.tsv.gz-meta.yaml")
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
#
#-------------------------------------------------------------------------------
# Build database for FinnGen 
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/FinnGen")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Download summary statistics from FinnGen 
#
url<-paste0("https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_G6_POSTVIRFAT.gz")
destfile<-"summary_stats_release_finngen_R12_G6_POSTVIRFAT.gz"
file_path<-file.path(current_dir,"Data/FinnGen/",destfile)
if(!file.exists(file_path)) {
  print("Downloading summary statistics from FinnGen")
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
# It calculates and plot the correlation table 
#-------------------------------------------------------------------------------
#
Correlations<-function(mydata,data.select,folder_path,n.cex,alpha) {
  #
  NC<-ncol(mydata)
  file.name<-paste0(data.select,"_corr.jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1000,height=1000)
  p_mat<-matrix(nrow=NC,ncol=NC) # build a matrix with p-values for correlations
  for (i in 1:NC) {
    for (j in 1:NC) {
      temp<-cor.test(mydata[,i],mydata[,j],alternative="two.sided",method="spearman")
      p_mat[i,j]<-temp$p.value
    }
  }
  rownames(p_mat)<-colnames(mydata[,c(1:NC)])
  colnames(p_mat)<-colnames(mydata[,c(1:NC)])
  corrplot(cor(mydata[,1:NC],use="pairwise.complete.obs",method="spearman"),is.corr=F, # plot the correlation table
           method="color",addCoef.col="black",number.digits=2,number.cex=n.cex,
           mar=c(0,0,4,0),p.mat=p_mat,sig.level=alpha,insig="blank",
           title=paste("Z scores"))
  dev.off()    
  #
}


