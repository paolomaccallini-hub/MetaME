# META_GWAS

## Abstract

META_GWAS is an R-based pipeline for running a genome-wide association study (GWAS) meta-analysis of chronic fatigue syndrome (CFS) summary statistics. The workflow automates the downloading of five public summary statistics: DecodeME, Million Veteran Project, UK Biobank from the Neale Lab, UK Biobank from the European Bioinformatics Institute, and FinnGen. The files are then brought to a common ground by column renaming, lift-over to GRCh38 and allele flipping (when necessary). The pipeline writes the five munged summary statistics along with a merged file (GWAS_FULL.tsv.gz), including all the variants that appear at least in one of them, with the annotations (p-value, Beta, OR, N, etc.). The user can set the samples to use for meta-GWAS analysis from a YAML file. In particular, I present here a meta-GWAS of 21,561 cases (European ancestry) built using DecodeME, MVP, and UK Biobank. Given the presence of different regression models across the summary statistics, the sample size-based method was employed to carry out the calculation of weighted Z scores for the variants of the meta-GWAS. 

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

The pipeline standardises the five summary statistics with [`format_sumstats`](https://github.com/neurogenomics/MungeSumstats), which munges allele columns and lifts coordinates from GRCh37 to GRCh38, if necessary. Munged sumstats are saved in `\Munged`. The column labels and their meaning are as follows:

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

Each sumstat can have additional columns (like **OR**, for logistic regression, or **LOG10P**). As an example, here you can see the rows of the munged sumstat of DecodeMe (`\Munged\DME_GRCh38.tsv.gz`) that contain the 10 most significant variants:

| SNP | CHR | BP | A1 | A2 | VARIANT_ID | FRQ | N | N_CAS | N_CON | BETA | SE | LOG10P | P | Z |
|----|-----|----|----|----|-------------|-----|------|--------|--------|--------|--------|--------|------------|-----------|
| rs6066909 | 20 | 48913376 | C | T | 20:48913376:C:T | 0.634113 | 275488 | 15579 | 259909 | 0.0904424 | 0.0133438 | 10.91380 | 1.219551e-11 | 6.777859 |
| rs6012555 | 20 | 48911205 | A | C | 20:48911205:A:C | 0.633910 | 275488 | 15579 | 259909 | 0.0901852 | 0.0133428 | 10.85740 | 1.388673e-11 | 6.759091 |
| rs6125539 | 20 | 49087273 | C | A | 20:49087273:C:A | 0.405455 | 275488 | 15579 | 259909 | -0.0807783 | 0.0130701 | 9.19416 | 6.394992e-10 | -6.180389 |
| rs6125576 | 20 | 49160382 | A | T | 20:49160382:A:T | 0.405720 | 275488 | 15579 | 259909 | -0.0806653 | 0.0130677 | 9.17348 | 6.706872e-10 | -6.172877 |
| rs4810909 | 20 | 49004835 | G | A | 20:49004835:G:A | 0.594451 | 275488 | 15579 | 259909 | 0.0804680 | 0.0130685 | 9.13114 | 7.393669e-10 | 6.157401 |
| rs755589 | 20 | 49104848 | G | A | 20:49104848:G:A | 0.405675 | 275488 | 15579 | 259909 | -0.0804335 | 0.0130691 | 9.12302 | 7.533209e-10 | -6.154479 |
| rs13044900 | 20 | 49146163 | G | A | 20:49146163:G:A | 0.405640 | 275488 | 15579 | 259909 | -0.0803478 | 0.0130667 | 9.10819 | 7.794890e-10 | -6.149051 |
| rs6125590 | 20 | 49204089 | G | T | 20:49204089:G:T | 0.405764 | 275488 | 15579 | 259909 | -0.0803416 | 0.0130664 | 9.10728 | 7.811240e-10 | -6.148717 |
| rs6063371 | 20 | 49202768 | G | A | 20:49202768:G:A | 0.405742 | 275488 | 15579 | 259909 | -0.0802805 | 0.0130656 | 9.09553 | 8.025461e-10 | -6.144417 |
| rs11697622 | 20 | 49110542 | G | A | 20:49110542:G:A | 0.406386 | 275488 | 15579 | 259909 | -0.0801362 | 0.0130626 | 9.06919 | 8.527270e-10 | -6.134782 |

### Merging

The five munged sumstats are merged into a file (`\Output\GWAS_FULL.tsv.gz`) with a row for each variant that is present in at least one of them. The annotations from each original sumstat are also included. As an example, these are the rows of `GWAS_FULL.tsv.gz` corresponding to the ten most significant variants of the DecodeME:


| SNP        | CHR | BP       | A1 | A2 | BETA_DME   | SE_DME    | P_DME        | FRQ_DME  | N_DME  | Neff_DME | Z_DME     | BETA_MVP     | SE_MVP     | P_MVP   | FRQ_MVP | N_MVP  | Neff_MVP | Z_MVP      | BETA_UKBNL   | SE_UKBNL    | P_UKBNL   | FRQ_UKBNL | N_UKBNL | Neff_UKBNL | Z_UKBNL    | BETA_UKBEIB  | SE_UKBEIB   | P_UKBEIB | FRQ_UKBEIB | N_UKBEIB | Neff_UKBEIB | Z_UKBEIB    | BETA_FG    | SE_FG     | P_FG     | FRQ_FG   | N_FG   | Neff_FG  | Z_FG       |
| ---------- | --- | -------- | -- | -- | ---------- | --------- | ------------ | -------- | ------ | -------- | --------- | ------------ | ---------- | ------- | ------- | ------ | -------- | ---------- | ------------ | ----------- | --------- | --------- | ------- | ---------- | ---------- | ------------ | ----------- | -------- | ---------- | -------- | ----------- | ----------- | ---------- | --------- | -------- | -------- | ------ | -------- | ---------- |
| rs6066909  | 20  | 48913376 | C  | T  | 0.0904424  | 0.0133438 | 1.219551e-11 | 0.634113 | 275488 | 58792    | 6.777859  | 0.032419891  | 0.02621097 | 0.21670 | 0.6284  | 443093 | 15427.33 | 1.2368827  | 2.13005e-04  | 0.000165199 | 0.1972670 | 0.365005  | 361141  | 6605.516   | 1.2893843  | 1.29759e-04  | 0.000139908 | 0.35     | 0.639675   | 484598   | 8331.876    | 0.92745947  | -0.0303995 | 0.0850934 | 0.720905 | 0.572202 | 463312 | 1131.309 | -0.3572486 |
| rs6012555  | 20  | 48911205 | A  | C  | 0.0901852  | 0.0133428 | 1.388673e-11 | 0.633910 | 275488 | 58792    | 6.759091  | 0.032419891  | 0.02621097 | 0.21690 | 0.6283  | 443093 | 15427.33 | 1.2368827  | 2.08068e-04  | 0.000165116 | 0.2076220 | 0.365159  | 361141  | 6605.516   | 1.2601323  | 1.24720e-04  | 0.000139831 | 0.37     | 0.639536   | 484598   | 8331.876    | 0.89193383  | -0.0329416 | 0.0850883 | 0.698649 | 0.572903 | 463312 | 1131.309 | -0.3871461 |
| rs6125539  | 20  | 49087273 | C  | A  | -0.0807783 | 0.0130701 | 6.394992e-10 | 0.405455 | 275488 | 58792    | -6.180389 | -0.039220713 | 0.02325622 | 0.09509 | 0.4024  | 443093 | 15427.33 | -1.6864614 | -1.53411e-04 | 0.000162279 | 0.3444790 | 0.404623  | 361141  | 6605.516   | -0.9453534 | -7.64157e-05 | 0.000137466 | 0.58     | 0.398730   | 484598   | 8331.876    | -0.55588800 | 0.0183845  | 0.0844899 | 0.827745 | 0.463853 | 463312 | 1131.309 | 0.2175941  |
| rs6125576  | 20  | 49160382 | A  | T  | -0.0806653 | 0.0130677 | 6.706872e-10 | 0.405720 | 275488 | 58792    | -6.172877 | -0.044016885 | 0.02391808 | 0.06681 | 0.3981  | 443093 | 15427.33 | -1.8403185 | -1.55553e-04 | 0.000162327 | 0.3379280 | 0.405388  | 361141  | 6605.516   | -0.9582694 | -7.53679e-05 | 0.000137508 | 0.58     | 0.399673   | 484598   | 8331.876    | -0.54809829 | 0.0189497  | 0.0844585 | 0.822472 | 0.463746 | 463312 | 1131.309 | 0.2243670  |
| rs4810909  | 20  | 49004835 | G  | A  | 0.0804680  | 0.0130685 | 7.393669e-10 | 0.594451 | 275488 | 58792    | 6.157401  | 0.038325114  | 0.02455894 | 0.11740 | 0.5898  | 443093 | 15427.33 | 1.5605363  | 1.26297e-04  | 0.000162224 | 0.4362540 | 0.404779  | 361141  | 6605.516   | 0.7785346  | 5.19962e-05  | 0.000137419 | 0.70     | 0.601067   | 484598   | 8331.876    | 0.37837708  | -0.0197586 | 0.0854544 | 0.815020 | 0.536524 | 463312 | 1131.309 | -0.2339558 |
| rs755589   | 20  | 49104848 | G  | A  | -0.0804335 | 0.0130691 | 7.533209e-10 | 0.405675 | 275488 | 58792    | -6.154479 | -0.039220713 | 0.02367336 | 0.09320 | 0.3977  | 443093 | 15427.33 | -1.6567447 | -1.50369e-04 | 0.000162206 | 0.3539150 | 0.405008  | 361141  | 6605.516   | -0.9270249 | -6.96557e-05 | 0.000137215 | 0.61     | 0.399988   | 484598   | 8331.876    | -0.50763911 | 0.0193008  | 0.0844842 | 0.819293 | 0.463646 | 463312 | 1131.309 | 0.2284546  |
| rs13044900 | 20  | 49146163 | G  | A  | -0.0803478 | 0.0130667 | 7.794890e-10 | 0.405640 | 275488 | 58792    | -6.149051 | -0.041141943 | 0.02375610 | 0.08127 | 0.3972  | 443093 | 15427.33 | -1.7318473 | -1.56058e-04 | 0.000162217 | 0.3360360 | 0.405221  | 361141  | 6605.516   | -0.9620323 | -7.70543e-05 | 0.000137407 | 0.57     | 0.399568   | 484598   | 8331.876    | -0.56077420 | 0.0189623  | 0.0844573 | 0.822353 | 0.463746 | 463312 | 1131.309 | 0.2245194  |
| rs6125590  | 20  | 49204089 | G  | T  | -0.0803416 | 0.0130664 | 7.811240e-10 | 0.405764 | 275488 | 58792    | -6.148717 | -0.045928932 | 0.02369368 | 0.05130 | 0.3984  | 443093 | 15427.33 | -1.9384465 | -1.25748e-04 | 0.000162506 | 0.4390470 | 0.405146  | 361141  | 6605.516   | -0.7738053 | -4.86869e-05 | 0.000137682 | 0.72     | 0.399426   | 484598   | 8331.876    | -0.35361848 | 0.0398295  | 0.0846524 | 0.637993 | 0.458601 | 463312 | 1131.309 | 0.4705064  |
| rs6063371  | 20  | 49202768 | G  | A  | -0.0802805 | 0.0130656 | 8.025461e-10 | 0.405742 | 275488 | 58792    | -6.144417 | -0.045928932 | 0.02369368 | 0.05130 | 0.3984  | 443093 | 15427.33 | -1.9384465 | -1.24475e-04 | 0.000162457 | 0.4435570 | 0.405135  | 361141  | 6605.516   | -0.7662027 | -5.02831e-05 | 0.000137632 | 0.71     | 0.399420   | 484598   | 8331.876    | -0.36534454 | 0.0394094  | 0.0846447 | 0.641512 | 0.458675 | 463312 | 1131.309 | 0.4655862  |
| rs11697622 | 20  | 49110542 | G  | A  | -0.0801362 | 0.0130626 | 8.527270e-10 | 0.406386 | 275488 | 58792    | -6.134782 | -0.045928932 | 0.02433473 | 0.05989 | 0.3906  | 443093 | 15427.33 | -1.8873818 | -1.51080e-04 | 0.000162202 | 0.3516310 | 0.406264  | 361141  | 6605.516   | -0.9314312 | -5.99007e-05 | 0.000136560 | 0.66     | 0.405742   | 484598   | 8331.876    | -0.43864016 | 0.0257011  | 0.0844750 | 0.760941 | 0.463756 | 463312 | 1131.309 | 0.3042450  |

Note that the suffix `_DME` stands for DecodeME, `_MVP` indicates Million Veteran Project and so forth. So, for instance, `SE_MPV` indicates the standard error of the regression coefficient from the Million Veteran Project sumstat.

### Meta-analysis

While DecodeME, MVP, and FinnGen rely on logistic regression for trait-variant marginal associations, the two UK Biobank GWASs employed a linear regression model. In a case like this, Inverse Variance Weighted (IVW) meta-analysis is not feasible, since BETAs and SEs from different regression models are not directly comparable. We must rely on zeta scores instead, using a sample-size weighted meta-analysis, as described in ([Willer 2010](https://academic.oup.com/bioinformatics/article/26/17/2190/198154)). Briefly, for each cohort *i*, a weight was assigned proportional to the square root of the effective sample size:

$$
w_i = \sqrt{N_{eff}^{(i)}}
$$

with $\ N_{eff}^{(i)}=\frac{4}{\frac{1}{N_{CAS}}+\frac{1}{N_{CON}}} $. The combined Z-score was computed as:

$$
Z = \frac{\sum_i w_i Z_i}{\sqrt{\sum_i w_i^2}}
$$

where $\ Z_1 $, $\ Z_2 $, … are the per-cohort Z statistics and $\ w_1 $, $\ w_2 $, … are their corresponding weights. These calculations are repeated for each variant that is present in a least one of the selected summary statistics. Next, 


- `META_main.R` – orchestrates package installation, summary statistic munging, and harmonisation steps for each cohort.
- `META_func.R` – helper functions that download input datasets, build the project directory structure, and load configuration values.
- `META_config.yml` – user-editable configuration file for p-value, minor allele frequency, and imputation quality cut-offs as well as cohort selection flags.
- `LICENSE` – licensing information for the project.

## Requirements

- R (≥ 4.0 recommended)
- Access to Bioconductor and CRAN to install dependencies
- Internet connection to retrieve cohort summary statistics from OSF, GWAS Catalog, and UK Biobank distribution endpoints

The main script installs the required packages automatically, including:

- `BiocManager`, `MungeSumstats`, `BSgenome.Hsapiens.NCBI.GRCh38`, `SNPlocs.Hsapiens.dbSNP155.GRCh38`, `SNPlocs.Hsapiens.dbSNP155.GRCh37`, and `BSgenome.Hsapiens.1000genomes.hs37d5`
- `httr`, `R.utils`, `data.table`, `gtexr`, `otargen`, `yaml`, `corrplot`

You may prefer to install these ahead of time to avoid repeated installations when running the pipeline on a cluster.

## Configuration

Edit `META_config.yml` to control filtering thresholds and which cohorts to include. The file exposes:

- `filters`: genome-wide significance thresholds for common and uncommon variants (`p_value_common`, `p_value_uncommon`), Hardy–Weinberg equilibrium cut-off (`hwe_p_value`), minimum allele frequency thresholds (`maf_common`, `maf_uncommon`), and an imputation INFO score cut-off (`info_cutoff`).
- `samples`: binary flags (1 = include, 0 = skip) for DecodeME, MVP, UK Biobank (Neale Lab), UK Biobank (EIB), and FinnGen cohorts.

Update these values before launching the analysis to ensure only the desired datasets are downloaded and processed.

## Running the pipeline

1. Open an R session in the repository root (or ensure the working directory is set to the project).
2. Run the main script:
   ```r
   Rscript META_main.R
   ```
3. The script creates required directories (`Data/`, `Munged/`) and downloads each cohort selected in `META_config.yml`.
4. Summary statistics are munged via `MungeSumstats` into harmonised `TSV.GZ` files stored under `Munged/`.

Depending on connectivity and download sizes, the first execution can take a while. Temporary outputs and downloaded data are preserved for reuse.

## Outputs

- Harmonised summary statistics per cohort saved in the `Munged/` directory (e.g., `Munged/DME_GRCh38.tsv.gz`).
- Downloaded raw summary statistics and metadata stored beneath `Data/` in cohort-specific subdirectories.

These munged files are intended for subsequent meta-analysis steps, which are not included in this repository.



## License

This project is released under the terms described in the included `LICENSE` file.
