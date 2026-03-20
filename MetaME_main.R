# file name: MetaME_main
#
#-------------------------------------------------------------------------------
# This script performs GWAS meta analysis for CFS using DecodeME (15,579 EUR), 
# MVP (3,891 EUR), UK Biobank EIB ()
#-------------------------------------------------------------------------------
#
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("MungeSumstats", quietly=TRUE)) BiocManager::install("MungeSumstats")
if (!requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly=TRUE)) BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh38", quietly=TRUE)) BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh37", quietly=TRUE)) BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
if (!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly=TRUE)) BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
#
#-------------------------------------------------------------------------------
# Link to module with functions
#-------------------------------------------------------------------------------
#
source("MetaME_func.R",echo=F)
#
#-------------------------------------------------------------------------------
# Add output folder, if absent
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Munged")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
#-------------------------------------------------------------------------------
# Read summary statistics of DecodeME main cohort, filter, munge, and save
# Assembly: GRCh38
# N cases: 15,579 (EUR)
# N controls: 259,909 (EUR)
# Regression: Logistic
# Reference: https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-/
#-------------------------------------------------------------------------------
#
# Read variants that passed quality filter
#
myQCEDvariants<-fread("Data/DecodeME/gwas_qced.var.gz")
#
# Read summary statistics
#
file_name<-"Data/DecodeME/gwas_1.regenie.gz" 
mydata<-fread(file_name)
head(mydata)
colnames(mydata)
#
# Remove columns we wont use
#
mydata<-mydata[,-"EXTRA"]
mydata<-mydata[,-"A1FREQ_CASES"]
mydata<-mydata[,-"A1FREQ_CONTROLS"]
mydata<-mydata[,-"TEST"]
mydata<-mydata[,-"CHISQ"]
#
# Keep only variants that passed INFO quality filter
#
mydata<-mydata[mydata$ID %in% myQCEDvariants$ID, ]
remove(myQCEDvariants)
#
# Filter by MAF (use cut-off of uncommon variants)
#
mydata<-mydata[mydata$A1FREQ<=1-MAFco_uc,]
mydata<-mydata[mydata$A1FREQ>=MAFco_uc,]
#
# Specify effect allele: in format_sumstats A1 is the non-effect allele, A2 is the effect allele,
# FRQ is the frequency of the effect allele 
# (https://www.bioconductor.org/packages/devel/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html).
# In DecodeME, ALLELE0 is	Non effect alleles (reference), ALLELE1	is Effect alleles (alternate),
# A1FREQ is	Effect allele frequency (from DecodeME readme)
#
for (j in 1:ncol(mydata)) {
  if (colnames(mydata)[j]=="ALLELE0") colnames(mydata)[j]<-"other_allele"
  if (colnames(mydata)[j]=="ALLELE1") colnames(mydata)[j]<-"effect_allele"
  if (colnames(mydata)[j]=="A1FREQ") colnames(mydata)[j]<-"effect_allele_frequency"
}
#
# Munge 
#
munge_path<-"Munged/DME_GRCh38.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,ref_genome="GRCh38",
                  compute_z="BETA",
                  save_path=munge_path, # where to save the sumstat
                  bi_allelic_filter=F, # keep non-biallelic SNPs
                  flip_frq_as_biallelic=T, # keep non-biallelic SNPs even when they need flipping 
                  log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                  log_folder="Munged") # were to save the log
} 
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-0
while(test<=30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$CHROM==mymunged$CHR[n1]&mydata$GENPOS==mymunged$BP[n1])
  if (length(n2)==1) {
    if (((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A2",with=F])&
         abs(mydata[n2,"BETA",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
        ((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A1",with=F])&
         abs(mydata[n2,"BETA",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
      test<-test+1
      next
    } else {
      stop("Munging generated an error!")
    }  
  }
}
remove(mymunged,mydata)
#
# Lift over to GRCh37 
#
mydata<-fread(munge_path)
munge_path<-"Munged/DME_GRCh37.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,
                  ref_genome="GRCh38",
                  convert_ref_genome="GRCh37",
                  save_path=munge_path, # where to save the sumstat
                  bi_allelic_filter=F, # keep non-biallelic SNPs
                  flip_frq_as_biallelic=T, # keep non-biallelic SNPs even when they need flipping 
                  log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                  log_folder="Munged") # were to save the log
} 
#
# Check for errors
# 
mydata<-fread("Munged/DME_GRCh38.tsv.gz")
mymunged<-fread("Munged/DME_GRCh37.tsv.gz")
test<-0
while(test<=30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$SNP==mymunged$SNP[n1])
  if (length(n2)==1) {
    if (((mydata[n2,"A2",with=F]==mymunged[n1,"A2",with=F])&
         abs(mydata[n2,"BETA",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
        ((mydata[n2,"A2",with=F]==mymunged[n1,"A1",with=F])&
         abs(mydata[n2,"BETA",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
      test<-test+1
      next
    } else {
      stop("Munging generated an error!")
    }  
  }
}
remove(mymunged,mydata)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Read summary statistics of MVP, filter, munge, and save
# Assembly: GRCh38 (readME)
# N cases: 3,891 (EUR)
# N controls: 439,202 (EUR)
# Regression: Logistic Mixed Model (SAIGE)
# Reference: https://pubmed.ncbi.nlm.nih.gov/39024449/
#-------------------------------------------------------------------------------
#
# Read summary statistics
#
file_name<-"Data/MVP/GCST90479178.tsv.gz" 
#
mydata<-fread(file_name,header=T,sep="\t")
head(mydata)
colnames(mydata)
#
# Remove a few columns we do not need
#
mydata<-mydata[,-"i2"]
mydata<-mydata[,-"case_af"]
mydata<-mydata[,-"control_af"]
mydata<-mydata[,-"r2"]
mydata<-mydata[,-"q_pval"]
#
# Filter by MAF, remove uncommon SNPs
#
mydata<-mydata[mydata$effect_allele_frequency<=1-MAFco_uc,]
mydata<-mydata[mydata$effect_allele_frequency>=MAFco_uc,]
#
# Calculate beta from odds ratio
#
mydata[,beta:=log(odds_ratio)]
#
# Calculate standard error of beta
#
mydata[,standard_error:=(log(ci_upper)-log(ci_lower))/(2*1.96)]
#
# Remove confidence interval, direction, and alt allele 
#
mydata<-mydata[,-"ci_upper"]
mydata<-mydata[,-"ci_lower"]
mydata<-mydata[,-"alt"]
mydata<-mydata[,-"direction"]
#
# Edit 
#
for (j in 1:ncol(mydata)) {
  if (colnames(mydata)[j]=="num_cases") colnames(mydata)[j]<-"N_CAS"
  if (colnames(mydata)[j]=="num_controls") colnames(mydata)[j]<-"N_CON"
}
#
# Munge
#
munge_path<-"Munged/MVP_GRCh38.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,ref_genome="GRCh38",
                  compute_z="BETA",
                  save_path=munge_path,
                  bi_allelic_filter=F, # keep non-biallelic SNPs
                  flip_frq_as_biallelic=T, # keep non-biallelic SNPs even when they need flipping 
                  log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                  log_folder="Munged") # were to save the log  
}
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-0
while(test<=30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$rsid==mymunged$SNP[n1])
  if (length(n2)==1) {
    if (((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A2",with=F])&
         abs(mydata[n2,"beta",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
        ((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A1",with=F])&
         abs(mydata[n2,"beta",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
      test<-test+1
      next
    } else {
      stop("Munging generated an error!")
    }  
  }
}
remove(mydata,mymunged)
#
# Lift over to GRCh37 
#
mydata<-fread(munge_path)
munge_path<-"Munged/MVP_GRCh37.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,
                  ref_genome="GRCh38",
                  convert_ref_genome="GRCh37",
                  save_path=munge_path, # where to save the sumstat
                  bi_allelic_filter=F, # keep non-biallelic SNPs
                  flip_frq_as_biallelic=T, # keep non-biallelic SNPs even when they need flipping 
                  log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                  log_folder="Munged") # were to save the log
} 
#
# Check for errors
# 
mydata<-fread("Munged/MVP_GRCh38.tsv.gz")
mymunged<-fread("Munged/MVP_GRCh37.tsv.gz")
test<-0
while(test<=30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$SNP==mymunged$SNP[n1])
  if (length(n2)==1) {
    if (((mydata[n2,"A2",with=F]==mymunged[n1,"A2",with=F])&
         abs(mydata[n2,"BETA",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
        ((mydata[n2,"A2",with=F]==mymunged[n1,"A1",with=F])&
         abs(mydata[n2,"BETA",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
      test<-test+1
      next
    } else {
      stop("Munging generated an error!")
    }  
  }
}
remove(mymunged,mydata)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Read summary statistics of UK Biobank (EIB), filter, munge, and save
# Assembly: GRCh37
# N cases: 2,092 (EUR)
# N controls: 482,506 (EUR)
# Regression: Linear 
# Reference: https://europepmc.org/article/MED/33959723
#-------------------------------------------------------------------------------
#
# Read summary statistics
#
file_name<-"Data/EIB/GCST90038694.tsv.gz" 
#
mydata<-fread(file_name,header=T,sep="\t")
head(mydata)
colnames(mydata)
#
# Remove a few columns we do not need
#
mydata<-mydata[,-"CHISQ_LINREG"]
mydata<-mydata[,-"CHISQ_BOLT_LMM_INF"]
mydata<-mydata[,-"P_BOLT_LMM_INF"]
mydata<-mydata[,-"CHISQ_BOLT_LMM"]
mydata<-mydata[,-"GENPOS"]
mydata<-mydata[,-"P_LINREG"]
#
# Filter by MAF, remove uncommon SNPs
#
mydata<-mydata[mydata$effect_allele_frequency<=1-MAFco_uc,]
mydata<-mydata[mydata$effect_allele_frequency>=MAFco_uc,]
#
# Filter by INFO
#
mydata<-mydata[mydata$INFO>=INFOco,]
mydata<-mydata[,-"INFO"]
#
# Add N
#
metadata<-read_yaml("Data/EIB/GCST90038694.tsv-meta.yaml")
mydata$N<-rep(metadata$samples[[1]]$sample_size,nrow(mydata))
mydata$N_cases<-rep(2092,nrow(mydata))
mydata$N_controls<-rep(482506,nrow(mydata))
#
# Munge
#
munge_path<-"Munged/UKBEIB_GRCh38.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,ref_genome="GRCh37",
                  convert_ref_genome="GRCh38",
                  compute_z="BETA",
                  bi_allelic_filter=F,
                  flip_frq_as_biallelic=T,
                  save_path=munge_path,
                  log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                  log_folder="Munged") # were to save the log
} 
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-0
while(test<=30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$variant_id==mymunged$SNP[n1])
  if (length(n2)==1) {
    if (((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A2",with=F])&
         abs(mydata[n2,"beta",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
        ((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A1",with=F])&
         abs(mydata[n2,"beta",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
      test<-test+1
      next
    } else {
      stop("Munging generated an error!")
    }  
  }
}
remove(mydata,mymunged)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Meta GWAS with METAL for GRCh38 sumstats
#-------------------------------------------------------------------------------
#
# Generate the script for METAL
#
lines<-c() # it contains the instruction for METAL
h<-0
#
# Input description and analysis description
#
h<-h+1
lines[h]<-"SEPARATOR TAB"
h<-h+1
lines[h]<-"MARKERLABEL SNP"
h<-h+1
lines[h]<-"ALLELELABELS A1 A2"
h<-h+1
lines[h]<-"PVALUELABEL P"
h<-h+1
lines[h]<-"EFFECTLABEL BETA"
h<-h+1
lines[h]<-"STDERRLABEL SE"
h<-h+1
lines[h]<-"FREQLABEL FRQ"
h<-h+1
lines[h]<-"WEIGHTLABEL Neff" # we use the effective size 
h<-h+1
lines[h]<-"SCHEME SAMPLESIZE" # use weighted zeta scores
h<-h+1
lines[h]<-"OVERLAP ON" # correct for samples overlap 
h<-h+1
lines[h]<-"VERBOSE OFF"
h<-h+1
lines[h]<-"AVERAGEFREQ ON"
h<-h+1
lines[h]<-"REMOVEFILTERS" # I want all the available SNPs!
#
# Input files
#
for (i in 1:length(populations)) {
  file_name<-paste0(current_dir,"/Munged/",populations[i],"_GRCh38.tsv")
  file_name<-sub("C:","/mnt/c",file_name)
  h<-h+1
  lines[h]<-paste("PROCESS",file_name)
}
#
# Output file
#
file_name<-paste0(current_dir,"/Output/GWAS_METAL_",paste0(populations,collapse="_"),"_GRCh38_ .tbl")
file_name<-sub("C:","/mnt/c",file_name)
h<-h+1
lines[h]<-paste("OUTFILE",file_name)
#
# Run the analysis
#
h<-h+1
lines[h]<-"ANALYZE"
h<-h+1
lines[h]<-"QUIT"
#
# Save METAL script as txt file
#
writeLines(lines,"metal_script.txt")
#
#-------------------------------------------------------------------------------
# Run METAL on GRCh38 sumstats
#-------------------------------------------------------------------------------
#
# Add output folder
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Output")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Extract files, if they are zipped, and add Neff for METAL
#
for (i in 1:length(populations)) {
  file_name_gz<-paste0(current_dir,"/Munged/",populations[i],"_GRCh38.tsv.gz")
  file_name<-sub(".gz","",file_name_gz)
  if(!file.exists(file_name)) {
    mydata<-fread(file_name_gz)
    mydata[, Neff := 4 / ((1 / N_CAS) + (1 / N_CON))]
    file_name<-sub(".gz","",file_name)
    fwrite(mydata,file_name,sep="\t")
    remove(mydata)
  }
}
#
# Run Metal
#
current_dir_win<-getwd()
current_dir_wsl<-gsub("C:","/mnt/c",current_dir_win)
script_path<-paste0(current_dir_wsl,"/metal_script.txt")
command<-paste("wsl",metal_path,script_path)
#
time.start<-as.numeric(Sys.time())
output<-system(command,wait=T,intern=T) # Run Exomiser
time.end<-as.numeric(Sys.time())
print(paste("METAL ended the analysis in",round((time.end-time.start)/60),"minutes"))
print(output)
#
# Remove unzipped files
#
for (i in 1:length(populations)) {
  file_name<-paste0(current_dir,"/Munged/",populations[i],"_GRCh38.tsv")
  if(file.exists(file_name)) {
    file.remove(file_name)
  }
}
#
#-------------------------------------------------------------------------------
# Read Metal meta-analysis, edit, and munge it for downstream analysis 
#-------------------------------------------------------------------------------
#
# Read summary statistics
#
file_name<-paste0(current_dir,"/Output/GWAS_METAL_",paste0(populations,collapse="_"),"_GRCh38_1.tbl")
mydata<-fread(file_name)
head(mydata)
colnames(mydata)
#
# Filter by MAF
#
mydata<-mydata[mydata$Freq1>=MAFco_uc,]
mydata<-mydata[mydata$Freq1<=(1-MAFco_uc),]
#
# Change names
#
for (j in 1:ncol(mydata)) {
  if (colnames(mydata)[j]=="MarkerName") colnames(mydata)[j]<-"SNP" 
  if (colnames(mydata)[j]=="Allele1") colnames(mydata)[j]<-"A1"
  if (colnames(mydata)[j]=="Allele2") colnames(mydata)[j]<-"A2" 
  if (colnames(mydata)[j]=="Freq1") colnames(mydata)[j]<-"FRQ"
  if (colnames(mydata)[j]=="Zscore") colnames(mydata)[j]<-"Z"
  if (colnames(mydata)[j]=="P-value") colnames(mydata)[j]<-"P"
  if (colnames(mydata)[j]=="Weight") colnames(mydata)[j]<-"Neff"
}
#
# Calculate SE (it is necessary to calculate BETA from Z, which is required for FUMA)
#
mydata[,SE:=(0.5*N*FRQ*(1-FRQ))^(-0.5)]
#
# Compute BETA
#
mydata[,BETA:=Z*SE]
#
# Munge in GRCh38
#
munge_path<-paste0(current_dir,"/Output/GWAS_METAL_",paste0(populations,collapse="_"),"_GRCh38.tsv.gz")
#
if (file.exists(munge_path)) file.remove(munge_path)
format_sumstats(mydata,ref_genome="GRCh38",
                bi_allelic_filter=F,
                flip_frq_as_biallelic=T,
                save_path=munge_path,
                log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                log_folder="Output") # were to save the log
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-0
while(test<=30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$SNP==mymunged$SNP[n1])
  if (length(n2)==1) {
    if (((toupper(mydata[n2,"A2",with=F])==mymunged[n1,"A2",with=F])&
         abs(mydata[n2,"Z",with=F]-mymunged[n1,"Z",with=F])<1e-6)|
        ((toupper(mydata[n2,"A2",with=F])==mymunged[n1,"A1",with=F])&
         abs(mydata[n2,"Z",with=F]+mymunged[n1,"Z",with=F])<1e-6)) {
      test<-test+1
      next
    } else {
      stop("Munging generated an error!")
    }  
  }
}
#
# Lift over to GRCh37
#
munge_path<-paste0(current_dir,"/Output/GWAS_METAL_",paste0(populations,collapse="_"),"_GRCh37.tsv.gz")
#
if (file.exists(munge_path)) file.remove(munge_path)
format_sumstats(mydata,ref_genome="GRCh38",
                convert_ref_genome="GRCh37",
                bi_allelic_filter=F,
                flip_frq_as_biallelic=T,
                save_path=munge_path,
                log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                log_folder="Output") # were to save the log  
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-0
while(test<=30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$SNP==mymunged$SNP[n1])
  if (length(n2)==1) {
    if (((toupper(mydata[n2,"A2",with=F])==mymunged[n1,"A2",with=F])&
         abs(mydata[n2,"Z",with=F]-mymunged[n1,"Z",with=F])<1e-6)|
        ((toupper(mydata[n2,"A2",with=F])==mymunged[n1,"A1",with=F])&
         abs(mydata[n2,"Z",with=F]+mymunged[n1,"Z",with=F])<1e-6)) {
      test<-test+1
      next
    } else {
      stop("Munging generated an error!")
    }  
  }
}
remove(mydata)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Plot Z scores for best associations 
#-------------------------------------------------------------------------------
#
# Read meta-analysis
#
file_name<-paste0(current_dir,"/Output/GWAS_METAL_",paste0(populations,collapse="_"),
                  "_GRCh38.tsv.gz")
mydata<-fread(file_name,sep="\t")
#
# Find best SNPs in meta analysis
#
mybest<-list()
n.best<-9
setorder(mydata,P)
mybest[[1]]<-mydata[1:n.best,]
names(mybest)[1]<-"METAL"
remove(mydata)
#
# Recover the best SNPs from the input GWAS
#
for (i in 1:length(populations)) {
  file_name_gz<-paste0(current_dir,"/Munged/",populations[i],"_GRCh38.tsv.gz")
  mydata<-fread(file_name_gz)
  df<-mybest[[1]][,c("SNP","CHR","BP","A1","A2")]
  mybest[[1+i]]<-merge(mydata,df,by=c("SNP","CHR","BP","A1","A2"),all.x=F)
  mybest[[1+i]]<-mybest[[1+i]][mybest[[1]]$SNP,on="SNP"] # ask for the same order in SNPs
  mybest[[1+i]][, Neff := 4 / ((1 / N_CAS) + (1 / N_CON))]
  names(mybest)[1+i]<-populations[i]
  remove(mydata)
}
#
# Plot Zs for a comparison between meta analysis and input GWASs
#
plots<-list()
for (i in 1:n.best) {
  Zeta<-data.frame(name=names(mybest),Z=rep(NA,length(mybest)),
                   FRQ=rep(NA,length(mybest)),
                   Neff=rep(NA,length(mybest)),P=rep(NA,length(mybest)))
  for (j in 1:length(mybest)) {
    Zeta$Z[j]<-mybest[[j]]$Z[i]
    Zeta$FRQ[j]<-mybest[[j]]$FRQ[i]
    if (j==1) {
      Zeta$Neff[j]<-round(mybest[[j]]$N[i]) # for the meta analysis N is the Neff corrected for overlap
    } else {
      Zeta$Neff[j]<-round(mybest[[j]]$Neff[i])
    }
    Zeta$P[j]<-mybest[[j]]$P[i]
  }
  Zeta$name<-factor(Zeta$name,levels=unique(Zeta$name))
  #
  # prepare the ID of the variant with hg38 coordinates and RSID 
  #
  variant<-paste0(c(mybest[[1]]$CHR[i],mybest[[1]]$BP[i],mybest[[1]]$A1[i],
                  mybest[[1]]$A2[i]),collapse=":")
  variant<-paste0(variant," (",mybest[[1]]$SNP[i],") P = ",Zeta$P[1])
  #
  # build plot 
  #
  plots[[i]]<-ggplot(Zeta,aes(x=name,y=Z)) +
    geom_hline(yintercept=0,linetype="dashed",colour="black",linewidth=0.8) +
    geom_segment(aes(x=name,xend=name,y=0,yend=Z)) +
    geom_point(size=4) +
    labs(x=NULL,y="Z-score") +
    labs(title=variant) +
    #
    # Add FRQ
    #
    geom_text(
      aes(label = sprintf("%.3f", FRQ)),
      vjust = -1.2,
      color = "blue",
      size = 4
    ) +
    #
    # Add Neff
    #
    geom_text(
      aes(label = format(round(Neff), big.mark = ",")),
      vjust = 1.8,
      color = "darkred",
      size = 4
    ) +
    
    theme_grey(base_size=14) +
    coord_flip() +
    theme(
      plot.title=element_text(size=14,hjust=0.5),
      axis.title.x=element_text(size=16),
      axis.text=element_text(size=16)
    )
}
#
# Plot results
#
file_name<-paste0(current_dir,"/Output/GWAS_METAL_",paste0(populations,collapse="_"),
                  "_GRCh38.pdf")
pdf(file_name,width=20,height=20)
wrap_plots(plots,ncol=3)
dev.off()