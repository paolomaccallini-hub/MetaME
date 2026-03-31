# MetaME

**MetaME** is an R-based pipeline for running GWAS meta-analyses of Myalgic Encephalomyelitis / Chronic Fatigue Syndrome (ME/CFS). It downloads publicly available summary statistics from multiple cohorts, harmonises and quality-controls them with MungeSumstats, runs sample-size-weighted fixed-effect meta-analyses with METAL, and produces per-variant Z-score comparison plots for the top associations. All downloads are conditional (skipped if the file already exists) and performed with automatic retry logic.

---

## Repository structure

```
MetaME/
├── MetaME_main.R        # Main orchestration script
├── MetaME_func.R        # Helper functions and cohort-specific munging routines
├── MetaME_config.yml    # QC filters and path to the METAL executable
├── Cohorts.csv          # Cohort registry with sample sizes
├── Meta_analyses.csv    # Binary matrix defining which cohorts enter each meta-analysis
└── LICENSE              # MIT licence
```

The pipeline creates three directories at runtime:

| Directory | Contents |
|-----------|----------|
| `Data/`   | Raw summary statistics downloaded from public repositories |
| `Munged/` | Harmonised, QC-filtered summary statistics (GRCh37 and GRCh38) |
| `Output/` | METAL results, munged meta-analysis files, Z-score PDF plots, METAL log |

---

## Cohorts

| Cohort | N_cases | N_controls | N_total | Neff | Ancestry | Regression | Source |
|--------|--------:|----------:|--------:|------:|----------|------------|--------|
| DME_1 | 15,579 | 259,909 | 275,488 | 58,792 | EUR | Logistic (REGENIE) | [DecodeME preprint](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-/) |
| DME_2 | 15,579 | 155,790 | 171,369 | 56,651 | EUR | Logistic (REGENIE) | DecodeME |
| DME_1_female | 12,833 | 218,949 | 231,782 | 48,490 | EUR | Logistic (REGENIE) | DecodeME |
| DME_1_male | 2,746 | 40,960 | 43,706 | 10,294 | EUR | Logistic (REGENIE) | DecodeME |
| DME_1_infectious_onset | 9,738 | 259,909 | 269,647 | 37,545 | EUR | Logistic (REGENIE) | DecodeME |
| MVP | 3,891 | 439,202 | 443,093 | 15,427 | EUR | Logistic Mixed (SAIGE) | [GWAS Catalog GCST90479178](https://www.ebi.ac.uk/gwas/studies/GCST90479178) |
| UKBEIB | 2,092 | 482,506 | 484,598 | 8,332 | EUR | Linear (BOLT-LMM) | [GWAS Catalog GCST90038694](https://www.ebi.ac.uk/gwas/studies/GCST90038694) |
| UKBNL_both_sexes | 1,659 | 359,482 | 361,141 | 6,606 | EUR | Linear | [Neale Lab UKB round 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9777867) |
| UKBNL_female | 1,208 | 192,945 | 194,153 | 4,802 | EUR | Linear | Neale Lab UKB round 2 |
| UKBNL_male | 451 | 166,537 | 166,988 | 1,799 | EUR | Linear | Neale Lab UKB round 2 |
| FG | 283 | 463,029 | 463,312 | 1,131 | FIN | Logistic | [FinnGen R12 G6_POSTVIRFAT](https://pubmed.ncbi.nlm.nih.gov/36653562/) |
| LC | 6,450 | 1,093,995 | 1,100,445 | 25,649 | EUR | Logistic | [Long COVID GWAS GCST90454541](https://air.unimi.it/handle/2434/1172544) |

**DME = DecodeME; MVP = Million Veteran Program; UKBEIB = UK Biobank (EBI release); UKBNL = UK Biobank (Neale Lab release); FG = FinnGen; LC = Long COVID.**

> **Note on sample overlap.** DecodeME and UK Biobank share controls drawn from the same population. Long COVID controls partially overlap with DecodeME, UK Biobank, and FinnGen. When cohorts with shared controls are combined, the pipeline automatically sets `OVERLAP ON` in the METAL script; sex-mismatched designs (e.g. DME female + UKB male) are treated as independent and receive no overlap correction.

---

## Meta-analyses

`Meta_analyses.csv` is a binary indicator matrix whose rows define named meta-analyses and whose columns are cohort flags. The pipeline iterates over every row and runs a full meta-analysis for each. The 16 pre-defined meta-analyses are:

| Meta-analysis | Cohorts included | Overlapping Controls |
|---------------|-----------------|--------------------|
| DME_1_MVP | DME_1, MVP | NO |
| DME_2_MVP | DME_2, MVP | NO |
| DME_1_MVP_FG | DME_1, MVP, FG | NO |
| DME_1_infectious_onset_LC | DME_1_infectious_onset, LC | YES |
| DME_1_MVP_LC | DME_1, MVP, LC | YES |
| DME_1_male_MVP | DME_1_male, MVP | NO |
| DME_2_UKBNL_both_sexes_MVP | DME_2, UKBNL_both_sexes, MVP | YES |
| DME_1_female_UKBNL_male_MVP | DME_1_female, UKBNL_male, MVP | NO |
| DME_1_female_UKBNL_male_FG | DME_1_female, UKBNL_male, FG | NO |
| DME_1_female_UKBNL_male | DME_1_female, UKBNL_male | NO |
| DME_1_male_UKBNL_female_MVP | DME_1_male, UKBNL_female, MVP | NO |
| DME_1_male_UKBNL_female | DME_1_male, UKBNL_female | NO |
| UKBEIB_MVP | UKBEIB, MVP | NO |
| UKBNL_both_sexes_MVP | UKBNL_both_sexes, MVP | NO |
| UKBNL_both_sexes_MVP_FG | UKBNL_both_sexes, MVP, FG | NO |
| MVP_LC | MVP, LC | NO |

The primary meta-analysis for the ME/CFS study is **DME_1_MVP**.

---

## Requirements

### Software

- **R** ≥ 4.4.1 (tested)
- **METAL** — the generic-metal binary must be compiled and its path set in `MetaME_config.yml`. On Windows the pipeline calls METAL via WSL; on Linux/macOS point directly to the binary.
- **Python** ≥ 3.12 (required internally by MungeSumstats for liftover)

### R packages

Installed automatically at startup if absent:

```r
# Bioconductor
MungeSumstats
BSgenome.Hsapiens.NCBI.GRCh38
SNPlocs.Hsapiens.dbSNP155.GRCh38
SNPlocs.Hsapiens.dbSNP155.GRCh37
BSgenome.Hsapiens.1000genomes.hs37d5

# CRAN
httr        # conditional downloads with retry logic
R.utils     # gzip compression
data.table  # fast I/O
yaml        # configuration parsing
ggplot2     # Z-score plots
patchwork   # multi-panel plot layout
```

---

## Installation

```bash
git clone https://github.com/paolomaccallini-hub/MetaME.git
cd MetaME
```

Edit `MetaME_config.yml` to set the path to the METAL executable:

```yaml
METAL:
  path_metal_exe: "/path/to/metal"
```

On Windows with WSL the path should be in WSL form, e.g. `/mnt/c/Users/<user>/generic-metal/build/bin/metal`.

---

## Usage

Open `MetaME_main.R` in RStudio (or run it from the command line with `Rscript`) from the repository root:

```bash
Rscript MetaME_main.R
```

The script reads `Meta_analyses.csv` and loops over every defined meta-analysis. For each:

1. Identifies which cohorts are needed.
2. Downloads missing raw summary statistics (conditional, with retry).
3. Munges and harmonises each cohort to both GRCh38 and GRCh37 (skipped if the munged file already exists).
4. Verifies munging integrity via a 30-sample random-check loop.
5. Runs METAL in sample-size-weighted Z-score scheme.
6. Post-processes the METAL output: renames columns, computes SE and BETA from Z, and munges the meta-analysis result in both assemblies.
7. Produces a PDF of Z-score lollipop plots for the top 9 variants, annotated with allele frequencies (blue) and effective sample sizes (dark red).

At the end of the full run, a summary pass over all munging logs populates `Cohorts.csv` (variant counts and genome-wide significant hits per cohort) and `Meta_analyses_2.csv` (the same for each meta-analysis).

---

## Pipeline details

### Quality control filters

Configured in `MetaME_config.yml`:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `maf_uncommon` | 0.01 | Minor allele frequency lower bound |
| `info_cutoff` | 0.90 | Imputation quality lower bound (INFO score) |
| `hwe_p_value` | 1×10⁻⁶ | Hardy–Weinberg equilibrium p-value lower bound |

### Cohort-specific preprocessing

Each cohort requires bespoke column renaming and, in some cases, derived quantities before entering MungeSumstats:

- **DecodeME** — Z-scores computed from BETA; QC-passed variants list (`gwas_qced.var.gz`) applied as an additional filter.
- **UKBEIB** — Z-scores computed from BETA; linear regression summary statistics lifted over from GRCh37 to GRCh38.
- **UKBNL** — Beta/SE from linear regression; phenotype N extracted from Neale Lab phenotype manifest; liftover from GRCh37 to GRCh38.
- **MVP** — Beta and SE derived from odds ratio and 95% confidence intervals (`log(OR)`, `(log(CI_upper)−log(CI_lower)) / (2×1.96)`).
- **FinnGen** — Already in GRCh38; Z-scores computed from BETA.
- **Long COVID** — Already in GRCh38; N cases and controls fixed from publication.

### METAL scheme

METAL is run with `SCHEME SAMPLESIZE` using the effective sample size (`Neff = 4 / (1/N_cases + 1/N_controls)`) as weight. `OVERLAP ON` is set automatically for cohort combinations that share controls (see sample overlap note above). `AVERAGEFREQ ON` is used to compute weighted allele frequencies across studies.

### Munging integrity checks

After every call to `MungeSumstats::format_sumstats`, 30 randomly sampled variants are cross-checked between the pre- and post-munging data frames to verify that effect allele orientation and BETA values are consistent. The pipeline stops with an error if any mismatch is detected.

---

## Output files

For each meta-analysis (example: `DME_1_MVP`):

| File | Description |
|------|-------------|
| `Output/GWAS_METAL_DME_1_MVP_GRCh38_1.tbl` | Raw METAL output (GRCh38) |
| `Output/GWAS_METAL_DME_1_MVP_GRCh38.tsv.gz` | Munged meta-analysis, GRCh38 |
| `Output/GWAS_METAL_DME_1_MVP_GRCh37.tsv.gz` | Munged meta-analysis, GRCh37 |
| `Output/GWAS_METAL_DME_1_MVP_GRCh38.pdf` | Z-score plots for top 9 variants |
| `Output/GWAS_METAL_DME_1_MVP_GRCh38_metal.log` | METAL console output |

Munged summary statistics contain: `SNP`, `CHR`, `BP`, `A1`, `A2`, `FRQ`, `Z`, `BETA`, `SE`, `P`, `N`, `N_CAS`, `N_CON`.

---

## Configuration reference

`MetaME_config.yml`:

```yaml
filters:
  maf_uncommon: 0.01      # lower MAF bound (variants below this are removed)
  info_cutoff:  0.9       # minimum imputation INFO score
  hwe_p_value:  1e-06     # minimum HWE p-value

METAL:
  path_metal_exe: "/path/to/metal"   # full path to the METAL executable
```

`Cohorts.csv` — read-only at pipeline start; updated at the end of the summary pass. Fields: `Cohort`, `N_cases`, `N_controls`, `N`, `Neff`, `rows` (post-QC variant count), `gws` (genome-wide significant variants at P < 5×10⁻⁸).

`Meta_analyses.csv` — binary indicator matrix; also updated at the end of the run with `rows` and `gws` columns per meta-analysis (saved as `Meta_analyses_2.csv`).

---

## Citation

If you use MetaME in your research, please cite the repository and the primary data sources listed in the Cohorts table above.

---

## License

This project is released under the [MIT License](LICENSE).
