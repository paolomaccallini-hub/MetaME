# META_GWAS

## Abstract

META_GWAS is an R-based pipeline for running a genome-wide association study (GWAS) meta-analysis of chronic fatigue syndrome (CFS) summary statistics. The workflow automates the downloading of five public summary statistics: DecodeME, Million Veteran Project, UK Biobank from the Neale Lab, UK Biobank from the European Bioinformatics Institute, and FinnGen. The files are then brought to a common ground by column renaming, lift-over to GRCh38 and allele flipping (when necessary). The pipeline writes the five munged summary statistics along with a merged file (GWAS_FULL.tsv.gz), including all the variants that appear at least in one of them, with the annotations (p-value, Beta, OR, N, etc.). Next, a meta-GWAS can be built from a subset of these five GWASs, specified by the user in a YAML file. In particular, I present here a meta-GWAS of 21,561 cases (European ancestry) built using DecodeME, MVP, and UK Biobank. Given the presence of different regression models across the summary statistics, the sample size-based method was employed to carry out the calculation of weighted Z scores for the variants of the meta-GWAS. The summary statistics are generated with respect to both GRCh38 and GRCh37. The latter was input to FUMA using the UK Biobank as the reference population to select candidate genes and perform tissue-level and cell-type analysis by FUMA's proprietary regression. Post-GWAS analysis revealed six risk loci, mapped to 41 candidate genes. Tissue-level regression associates ME/CFS with several brain regions, spanning the cerebellum, basal ganglia, and frontal cortex. Cell-type analysis reveals a significant regression between the genetic profile of the meta-GWAS (as defined by FUMA pipeline) and gene-expression in human neurons of the hippocampus, the Globus Pallidus, and the Amygdala. The regression with hippocampal neurons is confirmed in mouse brain. The main results of FUMA output are reported here, and the reader can explore and download the complete output at this link: [meta-GWAS analysis](https://fuma.ctglab.nl/snp2gene/666819). Fine-mapping of the meta-GWAS by a custom pipeline is undergoing and will be added to this repository. To date, this represents the biggest GWAS ever performed on ME/CFS patients.

## Methods

### Data sources

`META_main.R` retrieves and saves in the `\Data` folder the following five summary statistics.

| Database                 | Symbol     | cases       | controls    | Trait            | Regression  | Ancestry  | Assembly  | Reference     | Summary Statistics | 
| :----------------------- | :--------- | :---------- | :---------- | :--------------- | :---------- | :-------- | :-------- | :------------ | :----------------- |
| **DecodeME** | DME | 15579 | 259909 |CFS (CCC/IOM)  | Logistic | EUR | GRCh38 | ([Preprint_2025](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-)) | ([GWAS-1](https://osf.io/rgqs3/files/osfstorage))
| **MillionVeteranProject** | MVP | 3891 | 443093 | PheCode_798.1_CFS | Logistic | EUR | GRCh38 | ([Verma_2024](https://pubmed.ncbi.nlm.nih.gov/39024449/)) | ([GCST90479178](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479178/)) |
| **UKBiobank (NealeLab)** | UKBNL | 1659 | 359482 | self_reported CFS | Linear | EUR | GRCh37 | ([NealeLab](https://www.nealelab.is/uk-biobank)) | ([20002_1482](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=178908679#gid=178908679))
| **UKBiobank (EIB)**      | UKBEIB| 2092 | 482506 | self_reported CFS | Linear | EUR | GRCh37 | ([Dönertaş_2021](https://europepmc.org/article/MED/33959723)) | ([GCST90038694](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038694)) |
| **FinnGen** | FG | 283 | 663029 | Post-viral fatigue | Logistic | FIN | GRCh38 | ([Kurki_2023](https://pubmed.ncbi.nlm.nih.gov/36653562/)) | ([R12_G6_POSTVIRFAT](https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/)) |

Consult the respective project documentation for licensing or data use agreements that may apply to the downloaded data.

### Munging and Liftover

The pipeline standardises the five summary statistics using the R package MungeSumstats ([Murphy 2021](https://academic.oup.com/bioinformatics/article/37/23/4593/6380562?login=false)), which munges allele columns and lifts coordinates from GRCh37 to GRCh38, if necessary. Munged sumstats are saved in `\Munged`. The column labels and their meaning are as follows:

| Column | Description |
| :----- | :---------- |
| **SNP** | rs ID |
| **CHR** | chromosome number (GRCh38) |
| **BP** | base pair position (GRCh38) |
| **A1** | reference allele (GRCh38) |
| **A2** | effect allele |
| **Z** | Z-score |
| **BETA** | regression coefficient between trait and effect allele |
| **SE** | standard error |
| **P** | p-value |
| **N** | size (cases + controls) |
| **N_CAS** | number of cases |
| **N_CON** | number of controls |

Each sumstat can have additional columns (like **OR**, for logistic regression, or **LOG10P**). As an example, here you can see the rows of the munged sumstat of DecodeMe (`\Munged\DME_GRCh38.tsv.gz`) that contain the 5 most significant variants:

| SNP | CHR | BP | A1 | A2 | VARIANT_ID | FRQ | N | N_CAS | N_CON | BETA | SE | LOG10P | P | Z |
|----|-----|----|----|----|-------------|-----|------|--------|--------|--------|--------|--------|------------|-----------|
| rs6066909 | 20 | 48913376 | C | T | 20:48913376:C:T | 0.634113 | 275488 | 15579 | 259909 | 0.0904424 | 0.0133438 | 10.91380 | 1.219551e-11 | 6.777859 |
| rs6012555 | 20 | 48911205 | A | C | 20:48911205:A:C | 0.633910 | 275488 | 15579 | 259909 | 0.0901852 | 0.0133428 | 10.85740 | 1.388673e-11 | 6.759091 |
| rs6125539 | 20 | 49087273 | C | A | 20:49087273:C:A | 0.405455 | 275488 | 15579 | 259909 | -0.0807783 | 0.0130701 | 9.19416 | 6.394992e-10 | -6.180389 |
| rs6125576 | 20 | 49160382 | A | T | 20:49160382:A:T | 0.405720 | 275488 | 15579 | 259909 | -0.0806653 | 0.0130677 | 9.17348 | 6.706872e-10 | -6.172877 |
| rs4810909 | 20 | 49004835 | G | A | 20:49004835:G:A | 0.594451 | 275488 | 15579 | 259909 | 0.0804680 | 0.0130685 | 9.13114 | 7.393669e-10 | 6.157401 |

### Merging

The five munged sumstats are merged into a file (`\Output\GWAS_FULL.tsv.gz`) with a row for each variant that is present in at least one of them. The annotations from each original sumstat are also included. As an example, these are the rows of `GWAS_FULL.tsv.gz` corresponding to the five most significant variants of the DecodeME:


| SNP        | CHR | BP       | A1 | A2 | BETA_DME   | SE_DME    | P_DME        | FRQ_DME  | N_DME  | Neff_DME | Z_DME     | BETA_MVP     | SE_MVP     | P_MVP   | FRQ_MVP | N_MVP  | Neff_MVP | Z_MVP      | BETA_UKBNL   | SE_UKBNL    | P_UKBNL   | FRQ_UKBNL | N_UKBNL | Neff_UKBNL | Z_UKBNL    | BETA_UKBEIB  | SE_UKBEIB   | P_UKBEIB | FRQ_UKBEIB | N_UKBEIB | Neff_UKBEIB | Z_UKBEIB    | BETA_FG    | SE_FG     | P_FG     | FRQ_FG   | N_FG   | Neff_FG  | Z_FG       |
| ---------- | --- | -------- | -- | -- | ---------- | --------- | ------------ | -------- | ------ | -------- | --------- | ------------ | ---------- | ------- | ------- | ------ | -------- | ---------- | ------------ | ----------- | --------- | --------- | ------- | ---------- | ---------- | ------------ | ----------- | -------- | ---------- | -------- | ----------- | ----------- | ---------- | --------- | -------- | -------- | ------ | -------- | ---------- |
| rs6066909  | 20  | 48913376 | C  | T  | 0.0904424  | 0.0133438 | 1.219551e-11 | 0.634113 | 275488 | 58792    | 6.777859  | 0.032419891  | 0.02621097 | 0.21670 | 0.6284  | 443093 | 15427.33 | 1.2368827  | 2.13005e-04  | 0.000165199 | 0.1972670 | 0.365005  | 361141  | 6605.516   | 1.2893843  | 1.29759e-04  | 0.000139908 | 0.35     | 0.639675   | 484598   | 8331.876    | 0.92745947  | -0.0303995 | 0.0850934 | 0.720905 | 0.572202 | 463312 | 1131.309 | -0.3572486 |
| rs6012555  | 20  | 48911205 | A  | C  | 0.0901852  | 0.0133428 | 1.388673e-11 | 0.633910 | 275488 | 58792    | 6.759091  | 0.032419891  | 0.02621097 | 0.21690 | 0.6283  | 443093 | 15427.33 | 1.2368827  | 2.08068e-04  | 0.000165116 | 0.2076220 | 0.365159  | 361141  | 6605.516   | 1.2601323  | 1.24720e-04  | 0.000139831 | 0.37     | 0.639536   | 484598   | 8331.876    | 0.89193383  | -0.0329416 | 0.0850883 | 0.698649 | 0.572903 | 463312 | 1131.309 | -0.3871461 |
| rs6125539  | 20  | 49087273 | C  | A  | -0.0807783 | 0.0130701 | 6.394992e-10 | 0.405455 | 275488 | 58792    | -6.180389 | -0.039220713 | 0.02325622 | 0.09509 | 0.4024  | 443093 | 15427.33 | -1.6864614 | -1.53411e-04 | 0.000162279 | 0.3444790 | 0.404623  | 361141  | 6605.516   | -0.9453534 | -7.64157e-05 | 0.000137466 | 0.58     | 0.398730   | 484598   | 8331.876    | -0.55588800 | 0.0183845  | 0.0844899 | 0.827745 | 0.463853 | 463312 | 1131.309 | 0.2175941  |
| rs6125576  | 20  | 49160382 | A  | T  | -0.0806653 | 0.0130677 | 6.706872e-10 | 0.405720 | 275488 | 58792    | -6.172877 | -0.044016885 | 0.02391808 | 0.06681 | 0.3981  | 443093 | 15427.33 | -1.8403185 | -1.55553e-04 | 0.000162327 | 0.3379280 | 0.405388  | 361141  | 6605.516   | -0.9582694 | -7.53679e-05 | 0.000137508 | 0.58     | 0.399673   | 484598   | 8331.876    | -0.54809829 | 0.0189497  | 0.0844585 | 0.822472 | 0.463746 | 463312 | 1131.309 | 0.2243670  |
| rs4810909  | 20  | 49004835 | G  | A  | 0.0804680  | 0.0130685 | 7.393669e-10 | 0.594451 | 275488 | 58792    | 6.157401  | 0.038325114  | 0.02455894 | 0.11740 | 0.5898  | 443093 | 15427.33 | 1.5605363  | 1.26297e-04  | 0.000162224 | 0.4362540 | 0.404779  | 361141  | 6605.516   | 0.7785346  | 5.19962e-05  | 0.000137419 | 0.70     | 0.601067   | 484598   | 8331.876    | 0.37837708  | -0.0197586 | 0.0854544 | 0.815020 | 0.536524 | 463312 | 1131.309 | -0.2339558 |

Note that the suffix `_DME` stands for DecodeME, `_MVP` indicates Million Veteran Project and so forth. So, for instance, `SE_MPV` indicates the standard error of the regression coefficient from the Million Veteran Project sumstat. The pipeline write a table with the number of vairants that are shared by the five summary statistics in a pair-wise comparison (table below). It also builds a regression table for the zeta scores of a random sample of 1 million variants (Figure 1).

|       | DME | MVP | UKBNL | UKBEIB | FG |
|-------|------|------|--------|---------|--------|
| **DME** | 4,412,709 | 4,167,666 | 4,065,502 | 4,080,741 | 3,837,228 |
| **MVP** | 4,167,666 | 4,561,580 | 4,148,745 | 4,214,608 | 3,917,662 |
| **UKBNL** | 4,065,502 | 4,148,745 | 4,265,191 | 4,244,550 | 3,768,635 |
| **UKBEIB** | 4,080,741 | 4,214,608 | 4,244,550 | 4,356,852 | 3,815,038 |
| **FG** | 3,837,228 | 3,917,662 | 3,768,635 | 3,815,038 | 4,611,615 |

<p align="center">
  <img width="400" height="400" alt="Z-score correlations" src="https://github.com/user-attachments/assets/11159551-4de2-4fc5-911c-39b3700d4844" />
</p>

<p align="center"><em>Figure 1. Correlation between Z-scores from one million randomly selected variants across the five summary statistics.</em></p>

### Meta-analysis

While DecodeME, MVP, and FinnGen rely on logistic regression for trait-variant marginal associations, the two UK Biobank GWASs employed a linear regression model. In a case like this, Inverse Variance Weighted (IVW) meta-analysis is not feasible, since BETAs and SEs from different regression models are not directly comparable. We must rely on zeta scores instead, using a sample-size weighted meta-analysis, as described in ([Willer 2010](https://academic.oup.com/bioinformatics/article/26/17/2190/198154)). Briefly, for each cohort *i*, a weight was assigned proportional to the square root of the effective sample size:

$$
w_i = \sqrt{N_{eff}^{(i)}}
$$

with $\ N_{eff}^{(i)}=\frac{4}{\frac{1}{N_{CAS}}+\frac{1}{N_{CON}}} $. The combined Z-score was computed as:

$$
Z = \frac{\sum_i w_i Z_i}{\sqrt{\sum_i w_i^2}}
$$

where $\ Z_1 $, $\ Z_2 $, … are the per-cohort Z statistics and $\ w_1 $, $\ w_2 $, … are their corresponding weights. These calculations are repeated for each variant that is present in a least one of the selected summary statistics. Next, the p-value for each variant of the meta-GWAS was calculated:

$$
P=2\Phi(|Z|)
$$

To calculate the regression coefficients, I estimated the standard error as 

$$
SE=\left(\sum_iN_{eff}^{(i)}\right)^{-\frac{1}{2}}=\left(N_{eff}\right)^{-\frac{1}{2}}
$$

and I then used the well-known relation between BETA, Z, and SE: $\ BETA=SE \cdot Z$. 

### FUMA: gene-mapping and tissue analysis

Once I generated the Meta-GWAS in GRCh38, I transferred it to GRCh37 using MungeSumstats again. Both meta-GWAS are available in this repository, inside the folder `Output` with the following names: `GWAS_META_DME_MVP_UKBEIB_GRCh37.tsv.gz` and `GWAS_META_DME_MVP_UKBEIB_GRCh38.tsv.gz`. I submitted the GRCh37 summary statistics to FUMA (Functional Mapping and Annotation of GWAS), a web service that includes several modules for various stages of GWAS analysis. The SNP2GENE module is the first step: it identifies risk loci and performs gene mapping according to several criteria (positional mapping, eQTL-based and chromatin-based mapping) ([Watanabe 2017](https://academic.oup.com/bioinformatics/article/37/23/4593/6380562)). In SNP2GENE, I requested the UK Biobank reference population (white British), which is in GRCh37. I used GTEx v8 for gene mapping by eQTL, and included positional mapping. I requested MAGMA analysis, which is a necessary step for subsequently running the Cell Type module. I left the default values for the other parameters. Next, tissue-level analysis was assessed using MAGMA gene-property analysis as implemented in FUMA, which performs a regression of gene-level association Z-scores on tissue-specific gene expression levels (from GTEx and other datasets) to identify tissues in which genetic associations are overrepresented. The analysis is available at this link: [meta-GWAS analysis](https://fuma.ctglab.nl/snp2gene/666819).

### FUMA: cell-type analysis

I used the results from the SNP2GENE analysis as input to the Cell Type module ([Watanabe 2019](https://www.nature.com/articles/s41467-019-11181-1)), using as datasets the ensemble of 28 scRNA-seq datasets covering human and mouse brain. I focused on the brain because MAGMA tissue expression analysis over the 53 tissues of GTEx v8 suggests a significant regression of several anatomical regions of the central nervous system. I included steps 1, 2, and 3 of the standard cell type analysis (they are explained below). I requested a Bonferroni multiple comparison test. The analysis consists of three steps. In particular, the first step tests the significance of the estimate of $\ B_{E}$ in the following linear model, which is computed for each one of the cell types (across all the datasets selected, 28 in our case):

$$
Z_{gene_{i}}\,=\,B_{0}+E_{gene_i}^{c}B_{E}+\overline{E}_{gene_{i}}B_{A}+G_i^{(1)}B_1+G_i^{(2)}B_2+\cdot\cdot\cdot+G_i^{(n)}B_n
$$

for i=1, 2,..., N, where N is the total number of genes included in the analysis (usually around 18,000 genes), $\ Z_{gene_{i}}$ is the zeta-score computed for gene i by the SNP2GENE module, $\ E_{gene_i}^{c}$ is the gene expression of gene i in cell type c, $\ \overline{E}_{gene_{i}}$ is the average expression of of gene i across multiple cell types, and $\ G_i^{(n)}$ with j in 1, 2, ..., n are confounders such as gene length and correations between genes calculate from the LD matrices of the reference population ([Watanabe 2019](https://www.nature.com/articles/s41467-019-11181-1)). In step one, p-values are calculated for each cell type from the statistical test on the estimate $\ B_{E}$. These p-values are then corrected for multiple comparisons (Bonferroni) within each dataset. Sep three detects independent associations across all the datasets employed. 

## Results

### Risk loci and candidate genes

The SNP2GENE module of FUMA identifies six risk loci, reported below. You can explore and download the results of the analysis of FUMA on my meta-GWAS at this link: [meta-GWAS analysis](https://fuma.ctglab.nl/snp2gene/666819).

| GenomicLocus | uniqID | rsID | chr | pos | p | start | end | nSNPs | nGWASSNPs | nIndSigSNPs | IndSigSNPs | nLeadSNPs |
|---------------|--------|------|-----|-------------|------------------|----------|----------|--------|-------------|--------------|-------------|-----------|
| 1 | 2:64087270:A:T | rs4671073 | 2 | 64087270 | 2.0066161187e-08 | 63911761 | 64287632 | 71 | 27 | 1 | rs4671073 | 1 |
| 2 | 6:98403022:C:T | rs4363043 | 6 | 98403022 | 1.55029684749e-09 | 98310291 | 98547305 | 246 | 89 | 2 | rs4363043;rs9490512 | 1 |
| 3 | 12:118563202:C:T | rs12579722 | 12 | 118563202 | 3.85593332074e-08 | 118543707 | 118584594 | 18 | 10 | 1 | rs12579722 | 1 |
| 4 | 15:55158922:A:G | rs7165327 | 15 | 55158922 | 3.70443367717e-08 | 55137606 | 55184544 | 226 | 63 | 1 | rs7165327 | 1 |
| 5 | 20:47529913:C:T | rs6066909 | 20 | 47529913 | 5.78184785952e-11 | 47511792 | 47914180 | 209 | 59 | 2 | rs6066909;rs8119308 | 1 |
| 6 | 22:41418229:G:T | rs71327107 | 22 | 41418229 | 3.8457222389e-08 | 41215672 | 41713111 | 144 | 65 | 1 | rs71327107 | 1 |

These chromosomal regions map to 41 genes by positional criteria and eQTLs (GTEx v8):

| ensg | symbol | chr | start | end | strand | type | entrezID | HUGO |
|------|---------|-----|--------|--------|---------|----------------|-------|------|
| ENSG00000115507 | OTX1 | 2 | 63277192 | 63284971 | 1 | protein_coding | 5013 | OTX1 |
| ENSG00000143951 | WDPCP | 2 | 63348518 | 64054977 | -1 | protein_coding | 51057 | WDPCP |
| ENSG00000169764 | UGP2 | 2 | 64068074 | 64118696 | 1 | protein_coding | 7360 | UGP2 |
| ENSG00000143952 | VPS54 | 2 | 64119280 | 64246206 | -1 | protein_coding | 51542 | VPS54 |
| ENSG00000176834 | VSIG10 | 12 | 118501398 | 118573831 | -1 | protein_coding | 101928274 | VSIG10 |
| ENSG00000089220 | PEBP1 | 12 | 118573663 | 118583389 | 1 | protein_coding | 5037 | PEBP1 |
| ENSG00000135090 | TAOK3 | 12 | 118587606 | 118810750 | -1 | protein_coding | 51347 | TAOK3 |
| ENSG00000111707 | SUDS3 | 12 | 118814185 | 118855840 | 1 | protein_coding | 64426 | SUDS3 |
| ENSG00000124126 | PREX1 | 20 | 47240790 | 47444420 | -1 | protein_coding | 57580 | PREX1 |
| ENSG00000124198 | ARFGEF2 | 20 | 47538427 | 47653230 | 1 | protein_coding | 10564 | ARFGEF2 |
| ENSG00000124207 | CSE1L | 20 | 47662849 | 47713489 | 1 | protein_coding | 1434 | CSE1L |
| ENSG00000124214 | STAU1 | 20 | 47729878 | 47804904 | -1 | protein_coding | 6780 | STAU1 |
| ENSG00000124228 | DDX27 | 20 | 47835884 | 47860614 | 1 | protein_coding | 55661 | DDX27 |
| ENSG00000124201 | ZNFX1 | 20 | 47854483 | 47894963 | -1 | protein_coding | 57169 | ZNFX1 |
| ENSG00000128285 | MCHR1 | 22 | 41074754 | 41078818 | 1 | protein_coding | 2847 | MCHR1 |
| ENSG00000100372 | SLC25A17 | 22 | 41165634 | 41215403 | -1 | protein_coding | 10478 | SLC25A17 |
| ENSG00000100380 | ST13 | 22 | 41220539 | 41253026 | -1 | protein_coding | 6767 | ST13 |
| ENSG00000196236 | XPNPEP3 | 22 | 41253081 | 41363838 | 1 | protein_coding | 63929 | XPNPEP3 |
| ENSG00000100393 | EP300 | 22 | 41487790 | 41576081 | 1 | protein_coding | 2033 | EP300 |
| ENSG00000100395 | L3MBTL2 | 22 | 41601209 | 41627275 | 1 | protein_coding | 83746 | L3MBTL2 |
| ENSG00000100399 | CHADL | 22 | 41625517 | 41636938 | -1 | protein_coding | 150356 | CHADL |
| ENSG00000100401 | RANGAP1 | 22 | 41641615 | 41682255 | -1 | protein_coding | 5905 | RANGAP1 |
| ENSG00000269104 | AL035681.1 | 22 | 41685388 | 41685686 | -1 | protein_coding | NA | NA |
| ENSG00000100403 | ZC3H7B | 22 | 41697526 | 41756151 | 1 | protein_coding | 23264 | ZC3H7B |
| ENSG00000167074 | TEF | 22 | 41763337 | 41795330 | 1 | protein_coding | 7008 | TEF |
| ENSG00000100410 | PHF5A | 22 | 41855721 | 41864729 | -1 | protein_coding | 84844 | PHF5A |
| ENSG00000100412 | ACO2 | 22 | 41865129 | 41924993 | 1 | protein_coding | 50 | ACO2 |
| ENSG00000100413 | POLR3H | 22 | 41921808 | 41940610 | -1 | protein_coding | 171568 | POLR3H |
| ENSG00000172346 | CSDC2 | 22 | 41956767 | 41973745 | 1 | protein_coding | 27254 | CSDC2 |
| ENSG00000100417 | PMM1 | 22 | 41972898 | 41985894 | -1 | protein_coding | 5372 | PMM1 |
| ENSG00000100418 | DESI1 | 22 | 41994032 | 42017100 | -1 | protein_coding | 27351 | DESI1 |
| ENSG00000196419 | XRCC6 | 22 | 42017123 | 42060044 | 1 | protein_coding | 2547 | XRCC6 |
| ENSG00000100138 | NHP2L1 | 22 | 42069934 | 42086508 | -1 | protein_coding | 4809 | NHP2L1 |
| ENSG00000167077 | MEI1 | 22 | 42095503 | 42195460 | 1 | protein_coding | 150365 | MEI1 |
| ENSG00000100147 | CCDC134 | 22 | 42196683 | 42222303 | 1 | protein_coding | 79879 | CCDC134 |
| ENSG00000198911 | SREBF2 | 22 | 42229109 | 42303312 | 1 | protein_coding | 6721 | SREBF2 |
| ENSG00000159958 | TNFRSF13C | 22 | 42321045 | 42322822 | -1 | protein_coding | 115650 | TNFRSF13C |
| ENSG00000100162 | CENPM | 22 | 42334725 | 42343168 | -1 | protein_coding | 79019 | CENPM |
| ENSG00000183066 | WBP2NL | 22 | 42394729 | 42454460 | 1 | protein_coding | 164684 | WBP2NL |
| ENSG00000198951 | NAGA | 22 | 42454358 | 42466846 | -1 | protein_coding | 4668 | NAGA |
| ENSG00000100197 | CYP2D6 | 22 | 42522501 | 42526908 | -1 | protein_coding | 101929829 | CYP2D6 |

### Gene-set analysis

MAGMA-proprietary gene-set analysis identifies the Gene Ontology (cellular-component level) term `Glutamatergic synapse` as significant after Bonferroni correction:

| Gene Set | N genes | Beta | Beta STD | SE | P | Pbon |
|-----------|----------|--------|-----------|--------|-----------|-----------|
| GOCC_GLUTAMATERGIC_SYNAPSE | 386 | 0.19153 | 0.027396 | 0.042283 | 2.97×10⁻⁶ | 0.05056 |

### Tissue analysis

MAGMA-proprietary tissue analysis, based on a linear regression between zeta scores assigned to all the human genes from the GWAS analysis and tissue-specific gene-expression profiles, highlights several significant brain associations, spanning the basal ganglia, the cerebellum, and the cortex (Figure 2).

<figure align="center">
  <img width="1074" height="565" alt="MAGMA tissue analysis" src="https://github.com/user-attachments/assets/f6ebff7c-f3ec-4ddd-8e9d-212ac7dd6139">
  <figcaption><em>Figure 2. MAGMA tissue analysis showing significant associations across basal ganglia, cerebellum, and cortex.</em></figcaption>
</figure>

### Cell-type analysis

The results of FUMA cell-type analysis on the human and mouse brain are reported in the following table and in Figure 3. The genetic profile from the meta-GWAS regresses significantly with the gene-expression profile of human neurons from the external segment of the Globus pallidus, with the central nucleus of the Amygdala, with rostral CA3 of proper Hippocampus, in human brain. The regression with neurons of the hippocampus is confirmed in mouse brain (HC neurons). These results and the references to the original RNA-seq studies and to interactive brain atlases are collected in the following table.

| Dataset | Region | Level | Full Name | Species | Paper | Reference Atlas |
|:--------|:-------|:-----:|:---------|:--------|:-------|:-----------|
| Siletti  | GPe    | 1     | Globus_Pallidus_External_Segment| Human | ([Siletti 2023](https://pubmed.ncbi.nlm.nih.gov/37824663/)) | ([Allen Human Brain Atlas](https://atlas.brain-map.org/atlas?atlas=265297125#atlas=265297125&plate=112360888&structure=9001&x=40320&y=46978.1328125&zoom=-7&resolution=124.49&z=3) |
| Siletti  | CeN    | 1     | Amygdala_Central_Nucleus | Human | ([Siletti 2023](https://pubmed.ncbi.nlm.nih.gov/37824663/)) | ([Allen Brain Atlas](https://atlas.brain-map.org/atlas?atlas=265297125#atlas=265297125&plate=112360888&structure=9001&x=40320&y=46978.1328125&zoom=-7&resolution=124.49&z=3) |
| Siletti | CA3 | 1 | Hippocampus_Proper_Rostral_CA3 | Human | ([Siletti 2023](https://pubmed.ncbi.nlm.nih.gov/37824663/)) | ([Allen Brain Atlas](https://atlas.brain-map.org/atlas?atlas=265297125#atlas=265297125&plate=112360888&structure=9001&x=40320&y=46978.1328125&zoom=-7&resolution=124.49&z=3) |
| DropViz | HC | 1 | Hippocampus | Mouse | ([Saunders 2019](https://pubmed.ncbi.nlm.nih.gov/30096299/)) | ([DrpViz Atlas](http://dropviz.org/)) |

<figure align="center">
  <img width="1074" height="565" alt="MAGMA cell type analysis analysis" src=![Immagine 2025-10-13 173100](https://github.com/user-attachments/assets/24125377-d937-4394-864f-20e1b0d3a185)">
  <figcaption><em>Figure 3. MAGMA cell-type analysis showing significant associations with human neurons from the external segment of globus pallidus (GPe), the Central Nucleus of Amygdala (CeN), rostral CA3 of Proper Hippocampus (Hibb.CA3), and hippocampus from mouse (HC Neuron). The matrix shows no complete collinearity betwwen the significant regressions. For an explication of how to read the matrix generated by the third step, please visit [this resource](https://fuma.ctglab.nl/tutorial#workflow). As an introduction, consider that stars indicate the collinear covariates of the regression model discussed above, while the element of row i and column j indicates the proportional significance (PS) of cell type j conditioning on cell type i.</em></figcaption>
</figure>

## About the pipeline

### Repository components

- `META_main.R` – orchestrates package installation, summary statistic munging, and harmonisation steps for each cohort.
- `META_func.R` – helper functions that download input datasets, build the project directory structure, and load configuration values.
- `META_config.yml` – user-editable configuration file for p-value, minor allele frequency, and imputation quality cut-offs as well as cohort selection flags.
- `Output\GWAS_FULL.tsv.gz` - it has a row for each variant that appears in at least one of the five summary statistics used as input, with all the annotations.
- `Output\GWAS_META_DME_MVP_UKBEIB_GRCh37.tsv.gz` - GRCh37 summary statistics of the meta-GWAS obtained from DecodeME, Million Veteran Project, and UKBiobanl (European Institute of Bioinformatics)
- `Output\GWAS_META_DME_MVP_UKBEIB_GRCh37.tsv.gz` - GRCh38 summary statistics of the meta-GWAS obtained from DecodeME, Million Veteran Project, and UKBiobanl (European Institute of Bioinformatics)
- `LICENSE` – licensing information for the project.

### Requirements

- R (≥ 4.0 recommended)
- Access to Bioconductor and CRAN to install dependencies
- Internet connection to retrieve cohort summary statistics from OSF, GWAS Catalog, and UK Biobank distribution endpoints

The main script installs the required packages automatically, including:

- `BiocManager`, `MungeSumstats`, `BSgenome.Hsapiens.NCBI.GRCh38`, `SNPlocs.Hsapiens.dbSNP155.GRCh38`, `SNPlocs.Hsapiens.dbSNP155.GRCh37`, and `BSgenome.Hsapiens.1000genomes.hs37d5`
- `httr`, `R.utils`, `data.table`, `gtexr`, `otargen`, `yaml`, `corrplot`

You may prefer to install these ahead of time to avoid repeated installations when running the pipeline on a cluster.

### Configuration

Edit `META_config.yml` to control filtering thresholds and which cohorts to include. The file exposes:

- `filters`: genome-wide significance thresholds for common and uncommon variants (`p_value_common`, `p_value_uncommon`), Hardy–Weinberg equilibrium cut-off (`hwe_p_value`), minimum allele frequency thresholds (`maf_common`, `maf_uncommon`), and an imputation INFO score cut-off (`info_cutoff`).
- `samples`: binary flags (1 = include, 0 = skip) for DecodeME, MVP, UK Biobank (Neale Lab), UK Biobank (EIB), and FinnGen cohorts.

Update these values before launching the analysis to ensure only the desired datasets are downloaded and processed.

### Running the pipeline

1. Open an R session in the repository root (or ensure the working directory is set to the project).
2. Run the main script:
   
   ```r
   Rscript META_main.R
   ```
   
4. The script creates required directories (`Data/`, `Munged/`) and downloads each cohort selected in `META_config.yml`.
5. Summary statistics are munged via `MungeSumstats` into harmonised `TSV.GZ` files stored under `Munged/`.

Depending on connectivity and download sizes, the first execution can take a while. Temporary outputs and downloaded data are preserved for reuse.

### Outputs

- Harmonised summary statistics per cohort saved in the `Munged/` directory (e.g., `Munged/DME_GRCh38.tsv.gz`).
- Downloaded raw summary statistics and metadata stored beneath `Data/` in cohort-specific subdirectories.
- Meta-GWAS summary statistics in both GrCh37 and GRCh38 stored beneath `Output/`.

These munged files are intended for subsequent meta-analysis steps, which are not included in this repository.

### License

This project is released under the terms described in the included `LICENSE` file.
