# MetaME

R pipeline for GWAS meta-analysis of Myalgic Encephalomyelitis / Chronic Fatigue Syndrome (ME/CFS), using three European-ancestry cohorts: DecodeME, MVP, and UK Biobank (EBI).

This repository is the companion code to:

> *[citation placeholder]*

---

## Cohorts

| Cohort | Cases (EUR) | Controls (EUR) | Regression | Assembly | Reference |
|--------|-------------|----------------|------------|----------|-----------|
| DecodeME | 15,579 | 259,909 | Logistic | GRCh38 | [Rayner et al.](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-/) |
| MVP | 3,891 | 439,202 | Logistic Mixed Model (SAIGE) | GRCh38 | [Verma et al. 2024](https://pubmed.ncbi.nlm.nih.gov/39024449/) |
| UK Biobank EBI | 2,092 | 482,506 | Linear (BOLT-LMM) | GRCh37 | [Boyle et al.](https://europepmc.org/article/MED/33959723) |

Summary statistics are downloaded automatically by the pipeline on first run.

---

## Dependencies

### R packages

The following packages are installed automatically on first run if absent:

- `BiocManager`
- `MungeSumstats`
- `BSgenome.Hsapiens.NCBI.GRCh38`
- `SNPlocs.Hsapiens.dbSNP155.GRCh38`
- `SNPlocs.Hsapiens.dbSNP155.GRCh37`
- `BSgenome.Hsapiens.1000genomes.hs37d5`

The following CRAN packages must be installed manually:

```r
install.packages(c("httr", "R.utils", "data.table", "yaml", "ggplot2", "patchwork"))
```

### METAL

The pipeline uses METAL for meta-analysis, called via WSL from RStudio on Windows.

1. Download the METAL source from the [2020-05-05 release](https://github.com/statgen/METAL/releases/tag/2020-05-05)
2. Compile within WSL2:
   ```bash
   mkdir build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make
   ```
3. Set the path to the compiled binary in `MetaME_config.yml`:
   ```yaml
   METAL:
     path_metal_exe: "/mnt/c/your/path/to/metal"
   ```

---

## Configuration

All pipeline settings are controlled via `MetaME_config.yml`:

```yaml
filters:
  maf_uncommon: 0.01   # minimum MAF threshold
  info_cutoff: 0.9     # minimum imputation quality (INFO)

samples:
  DME: 1               # 1 = include, 0 = exclude
  MVP: 1
  UKBEIB: 1

METAL:
  path_metal_exe: "/mnt/c/your/path/to/metal"
```

---

## Usage

Open `MetaME_main.R` in RStudio and run it from top to bottom. The script will:

1. Install missing R packages
2. Source `MetaME_func.R`, which reads the config and downloads all summary statistics
3. Filter, munge, and liftover each cohort's summary statistics
4. Generate and run the METAL script
5. Post-process the METAL output and produce GRCh38 and GRCh37 versions
6. Generate Z-score comparison plots for the top associations

Set the R working directory to the repository root before running.

---

## Outputs

All outputs are written to subdirectories created automatically:

| Directory | Contents |
|-----------|----------|
| `Data/` | Downloaded raw summary statistics |
| `Munged/` | Filtered and munged per-cohort summary statistics (GRCh38 and GRCh37) |
| `Output/` | METAL meta-analysis results and final munged sumstats |

Key output files:

| File | Description |
|------|-------------|
| `Output/GWAS_METAL_DME_MVP_UKBEIB_GRCh38.tsv.gz` | Meta-analysis summary statistics, GRCh38 |
| `Output/GWAS_METAL_DME_MVP_UKBEIB_GRCh37.tsv.gz` | Meta-analysis summary statistics, GRCh37 (recommended for web-FUMA) |
| `Output/GWAS_METAL_DME_MVP_UKBEIB_GRCh38.pdf` | Z-score plots for top associations |

Effect sizes (BETA) and standard errors (SE) in the output are reconstructed from METAL Z-scores using the formula of Vukcevic et al. (2011, *Genetic Epidemiology*). Z-scores and p-values are exact as computed by METAL.

---

## FUMA

Functional mapping and annotation of the meta-analysis results was performed using [FUMA v1.5.2](https://fuma.ctglab.nl) (SNP2GENE). The input file was `GWAS_METAL_DME_MVP_UKBEIB_GRCh37.tsv.gz` (GRCh37). Key settings:

| Parameter | Value |
|-----------|-------|
| Genome build | GRCh37 |
| Lead SNP threshold | P < 5×10⁻⁸ |
| GWAS-wide SNP inclusion | P < 0.05 |
| LD r² (lead SNPs) | 0.6 |
| LD r² (locus boundaries) | 0.1 |
| Locus merge distance | 250 kb |
| MHC region | Excluded |
| Reference panel | 1000 Genomes Phase 3, EUR |
| Gene type | Protein-coding |
| Ensembl | v102 |
| eQTL mapping | GTEx v8, all 49 tissues |
| MAGMA | Enabled, 0 kb window |
| MAGMA expression | GTEx v8 tissue-specific and general avg log2TPM |

The full FUMA parameter file is available in [params.config](FUMA/params.config).

---

## Environment

Developed and tested on:
- Windows 11 24H2 (build 10.0.26200.8037)
- RStudio 2026.1.1.403
- R 4.4.1 (2024-06-14 ucrt)
- WSL 2.3.24.0, kernel 5.15.153.1-2, Ubuntu 24.04.3 LTS

---

## License

[MIT License](LICENSE)
