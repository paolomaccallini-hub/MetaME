# file name: META_main
#
#-------------------------------------------------------------------------------
# This script performs GWAS meta analysis for CFS using DecodeME (15579 EUR), 
# MVP (3891 EUR), UK Biobank NealeLab(1659 EUR), UK Bioban EIB (), Finngen (283 FIN)
#-------------------------------------------------------------------------------
#
install.packages('BiocManager')
BiocManager::install("MungeSumstats",force=TRUE)
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
#
#-------------------------------------------------------------------------------
# Link to module with functions
#-------------------------------------------------------------------------------
#
source("META_func.R",echo=F)
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
# N cases: 15579 (EUR)
# N controls: 259909 (EUR)
# Regression: Logistic
# Reference: https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-
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
# Filter by MAF
#
mydata<-mydata[mydata$A1FREQ<=1-MAFco_uc,]
mydata<-mydata[mydata$A1FREQ>=MAFco_uc,]
#
# Specify effect allele: in format_sumstats A1 is the non-effect allele, A2 is the effect allele,
# FRQ is the frequency of the effect allele 
# (https://www.bioconductor.org/packages/devel/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html).
# In DecodeME ALLELE0	Non effect alleles (reference), ALLELE1	Effect alleles (alternate),
# A1FREQ	Effect allele frequencies (from DecodeME readme)
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
  format_sumstats(mydata,ref_genome="GRCh38",compute_z="BETA",save_path=munge_path)  
} 
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-1
while(test<30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$CHROM==mymunged$CHR[n1]&mydata$GENPOS==mymunged$BP[n1])
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
remove(mydata,mymunged)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Read summary statistics of MVP, filter, munge, and save
# Assembly: GRCh38 (readME)
# N cases: 3891 (EUR)
# N controls: 443093 (EUR)
# Regression: Linear Mixed Model
# Reference: https://pubmed.ncbi.nlm.nih.gov/39024449/
#-------------------------------------------------------------------------------
#
# Read summary statistics
#
file_name<-"DATA/MVP/GCST90479178.tsv.gz" 
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
# Filter by MAF
#
mydata<-mydata[mydata$effect_allele_frequency<=1-MAFco_uc,]
mydata<-mydata[mydata$effect_allele_frequency>=MAFco_uc,]
#
# Calculate beta from odds ratio
#
mydata[,beta:=log(odds_ratio)]
#
# Calculate standard error
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
  format_sumstats(mydata,ref_genome="GRCh38",compute_z="BETA",save_path=munge_path)  
}
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-1
while(test<30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$rsid==mymunged$SNP[n1])
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
remove(mydata,mymunged)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Read summary statistics of UK Biobank main CFS cohort (Neale Lab)
# Assembly: GRCh37
# N cases: 1659 (EUR)
# N controls: 359482 (EUR)
# Regression: Linear 
# Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC9777867
#-------------------------------------------------------------------------------
#
# Read summary statistics
#
file_name<-"Data/NealeLab/20002_1482.gwas.imputed_v3.both_sexes.tsv.bgz" 
mydata<-fread(file_name)
head(mydata)
colnames(mydata)
#
# Read variant annotations
#
all_variants<-fread("Data/NealeLab/variants.tsv.bgz",header=TRUE,sep="\t")
#
# Remove columns we wont use
#
mydata<-mydata[,-"AC"]
mydata<-mydata[,-"ytx"]
mydata<-mydata[,-"expected_case_minor_AC"]
#
# Add annotation
#
mydata<-merge(mydata,all_variants,by="variant",all=T)
remove(all_variants)
#
# Filter by MAF
#
mydata<-mydata[mydata$minor_AF>=MAFco_uc,]
#
# Filter by info
#
mydata<-mydata[mydata$info>=INFOco,]
mydata<-mydata[,-"info"]
#
# Remove low confidence variants
#
mydata<-mydata[low_confidence_variant==F,]
mydata<-mydata[,-"low_confidence_variant"]
#
# Remove variants outside HWE
#
mydata<-mydata[p_hwe>pco_HWE]
mydata<-mydata[,-"p_hwe"]
#
# Split columns with coordinates and alleles
# Note that in this sumstat the alternative allele is the effect allele:
# https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=227859291#gid=227859291
#
mydata[,c("CHR","BP","other_allele","effect_allele"):=tstrsplit(variant,":",fixed=T)]
mydata[,BP:=as.numeric(BP)]
mydata<-mydata[,-"variant"]
#
# Add N
#
phenotypes<-fread("Data/NealeLab/phenotypes.both_sexes.v2.tsv.bgz")
phenotypes<-phenotypes[phenotype=="20002_1482"]
N_cases<-phenotypes$n_cases
N_controls<-phenotypes$n_controls
N<-phenotypes$n_non_missing
mydata$N<-rep(N,nrow(mydata))
mydata$N_cases<-rep(N_cases,nrow(mydata))
mydata$N_controls<-rep(N_controls,nrow(mydata))
#
# Munge 
#
munge_path<-"Munged/UKBNL_GRCh38.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,ref_genome="GRCh37",convert_ref_genome="GRCh38",compute_z="BETA",save_path=munge_path)  
} 
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-1
while(test<30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$rsid==mymunged$SNP[n1])
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
remove(mydata,mymunged)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Read summary statistics of UK Biobank (EIB), filter, munge, and save
# Assembly: GRCh37
# N cases: 2092 (EUR)
# N controls: 482506 (EUR)
# Regression: Linear 
# Reference: https://europepmc.org/article/MED/33959723
#-------------------------------------------------------------------------------
#
# Read summary statistics
#
file_name<-"DATA/EIB/GCST90038694.tsv.gz" 
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
# Filter by MAF
#
mydata<-mydata[mydata$effect_allele_frequency<=1-MAFco_uc,]
mydata<-mydata[mydata$effect_allele_frequency>=MAFco_uc,]
#
# Filter by INFO
#
mydata<-mydata[mydata$INFO>INFOco,]
mydata<-mydata[,-"INFO"]
#
# Add N
#
metadata<-read_yaml("DATA/EIB/GCST90038694.tsv.gz-meta.yaml")
mydata$N<-rep(metadata$samples[[1]]$sample_size,nrow(mydata))
mydata$N_cases<-rep(2092,nrow(mydata))
mydata$N_controls<-rep(482506,nrow(mydata))
#
# Munge
#
munge_path<-"Munged/UKBEIB_GRCh38.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,ref_genome="GRCh37",convert_ref_genome="GRCh38",compute_z="BETA",save_path=munge_path)  
} 
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-1
while(test<30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$variant_id==mymunged$SNP[n1])
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
remove(mydata,mymunged)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Read summary statistics of Finngen, filter, munge, and save
# Assembly: GRCh38
# N cases: 283 (FIN)
# N controls: 663029 (FIN)
# Regression: Logistic
# Reference: https://pubmed.ncbi.nlm.nih.gov/36653562/
#-------------------------------------------------------------------------------
#
# Read summary statistics
#
file_name<-"Data/FinnGen/summary_stats_release_finngen_R12_G6_POSTVIRFAT.gz" # genome_assembly: GRCh38
mydata<-fread(file_name)
head(mydata)
colnames(mydata)
#
# Remove columns we wont use
#
mydata<-mydata[,-"af_alt_cases"]
mydata<-mydata[,-"af_alt_controls"]
mydata<-mydata[,-"nearest_genes"]
mydata<-mydata[,-"mlogp"]
#
# Filter by MAF
#
mydata<-mydata[mydata$af_alt>=MAFco_uc,]
mydata<-mydata[mydata$af_alt<=(1-MAFco_uc),]
#
# Calculate OR
#
mydata$OR<-exp(mydata$beta)
#
# Edit column names
#
for (j in 1:ncol(mydata)) {
  if (colnames(mydata)[j]=="#chrom") colnames(mydata)[j]<-"CHR" 
  if (colnames(mydata)[j]=="sebeta") colnames(mydata)[j]<-"SE"
  if (colnames(mydata)[j]=="af_alt") colnames(mydata)[j]<-"effect_allele_frequency" 
  if (colnames(mydata)[j]=="rsids") colnames(mydata)[j]<-"SNP"
  if (colnames(mydata)[j]=="ref") colnames(mydata)[j]<-"other_allele"
  if (colnames(mydata)[j]=="alt") colnames(mydata)[j]<-"effect_allele"
}
#
# Add sample size
#
mydata$N<-rep(463312,nrow(mydata))
mydata$N_cases<-rep(283,nrow(mydata))
mydata$N_controls<-rep(463029,nrow(mydata))
#
# Munge 
#
munge_path<-"Munged/FG_GRCh38.tsv.gz"
#
if (!file.exists(munge_path)) {
  format_sumstats(mydata,ref_genome="GRCh38",compute_z="BETA",save_path=munge_path)  
} 
#
# Check for errors
# 
mymunged<-fread(munge_path)
test<-1
while(test<30) {
  n1<-sample(c(1:nrow(mymunged)),1)
  n2<-which(mydata$SNP==mymunged$SNP[n1])
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
remove(mydata,mymunged)
gc() # free unused memory
#
#-------------------------------------------------------------------------------
# Read munged summary statistics and collect all data in a single file
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Output")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
GWAS_names<-names(samples)
GWAS_list<-list()
for (i in 1:length(GWAS_names)) {
  file_name<-paste0("Munged/",GWAS_names[i],"_GRCh38.tsv.gz")
  GWAS_list[[i]]<-fread(file_name)
  #
  # Calculate Neff
  #
  GWAS_list[[i]][, Neff := 4 / ((1 / N_CAS) + (1 / N_CON))]
  #
  # Select columns
  #
  GWAS_list[[i]]<-GWAS_list[[i]][,c("SNP","CHR","BP","A1","A2","BETA","SE","P",
                                    "FRQ","N","Neff","Z")]
  #
  # Edit column names
  #
  NI<-GWAS_names[i]
  colnames(GWAS_list[[i]])<-c("SNP","CHR","BP","A1","A2",paste0("BETA_",NI),
                            paste0("SE_",NI),paste0("P_",NI),paste0("FRQ_",NI),
                            paste0("N_",NI),paste0("Neff_",NI),paste0("Z_",NI))
  names(GWAS_list)[i]<-GWAS_names[i]
}
#
# Build a table with number of shared SNPs
#
matSNPs<-matrix(data=NA,nrow=length(GWAS_names),ncol=length(GWAS_names))
rownames(matSNPs)<-GWAS_names
colnames(matSNPs)<-GWAS_names
for (i in 1:length(GWAS_names)) {
  for (j in 1:length(GWAS_names)) {
    index<-which(GWAS_list[[i]]$SNP%in%GWAS_list[[j]]$SNP)
    matSNPs[i,j]<-length(index)
  }
}
write.csv(matSNPs,"Output/Shared_SNPs.csv",row.names=T,quote=F)
#
# Write a file with all the SNPs, even if present in only some of the GWAS
#
GWAS_FULL<-Reduce(function(x,y) merge(x,y,by=c("SNP","CHR","BP","A1","A2"),all=T),GWAS_list)
remove(GWAS_list)
gc() # free unused memory
#
# Correlation table for Z scores
#
Z_cols<-grep("Z_",colnames(GWAS_FULL),value=T)
df<-as.data.frame(GWAS_FULL[sample(seq(1,nrow(GWAS_FULL)),2000000),Z_cols,with=F])
Correlations(df,data.select="Z_scores",folder_path="Output",n.cex=2,alpha=0.05)
remove(df)
#
# Save the result
#
fwrite(GWAS_FULL,"Output/GWAS_FULL.tsv.gz",sep="\t")
remove(GWAS_FULL)
#
#-------------------------------------------------------------------------------
# Meta GWAS 
#-------------------------------------------------------------------------------
#
# Read data
#
GWAS_FULL<-fread("Output/GWAS_FULL.tsv.gz")
#
# Select columns on selected populations
#
index<-c()
for (pop in populations) {
  index<-c(index,which(grepl(pop,colnames(GWAS_FULL))))
}
CLM<-colnames(GWAS_FULL)[index]
CLM<-c("SNP","CHR","BP","A1","A2",CLM)
GWAS_META<-GWAS_FULL[,CLM,with=F]
remove(GWAS_FULL)
#
# Compute N
#
n_cols<-grep("N_",colnames(GWAS_META),value=T)
GWAS_META[,N:=rowSums(GWAS_META[,n_cols,with=F],na.rm=T)]
#
# Remove variants with no data
#
GWAS_META<-GWAS_META[N>0,]
#
# Compute weights (w = sqrt(Neff))
# 
for (i in 1:length(populations)) {
  index<-grep(paste0("Neff_",populations[i]),colnames(GWAS_META),value=TRUE)
  GWAS_META[,paste0("w_",populations[i]):=GWAS_META[,index,with=F]^0.5]
}
#
# Compute weighted Z ( Z = (w1*Z1 + w2*Z2 + w3*Z3 + ...)/(w1^2 + w2^2 + w3^2 + ...)^0.5 )
# 
w_cols<-grep("w_",colnames(GWAS_META),value=T)
Z_cols<-grep("Z_",names(GWAS_META),value=T)
GWAS_META[,Z:=rowSums(GWAS_META[,w_cols,with=F]*GWAS_META[,Z_cols,with=F],na.rm=T)/
            (rowSums(GWAS_META[,w_cols,with=F]^2,na.rm=T))^0.5]
#
# Compute p-value
#
GWAS_META[,P:=2*pnorm(-abs(Z))]
#
# Compute Neff
#
Neff_cols<-grep("Neff_",colnames(GWAS_META),value=T)
GWAS_META[,Neff:=rowSums(GWAS_META[,Neff_cols,with=F],na.rm=T)]
#
# Compute SE
#
GWAS_META[,SE:=Neff^(-0.5)]
#
# Compute BETA
#
GWAS_META[,BETA:=Z*SE]
#
# Save 
#
file_name<-paste0("Output/GWAS_META_",paste0(populations,collapse="_"),"_GRCh38.tsv.gz")
fwrite(GWAS_META[,c("SNP","CHR","BP","A1","A2","BETA","SE","P","N","Z")],file_name,sep="\t")
#
# Lift over from GRCh38 to GRCh37 and save
#
munge_path<-paste0("Output/GWAS_META_",paste0(populations,collapse="_"),"_GRCh37.tsv.gz")
#
if (!file.exists(munge_path)) {
  format_sumstats(GWAS_META[,c("SNP","CHR","BP","A1","A2","BETA","SE","P","N","Z")],
                  ref_genome="GRCh38",convert_ref_genome="GRCh37",save_path=munge_path)  
} 

