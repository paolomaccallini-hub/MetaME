# META_GWAS

META_GWAS is an R-based pipeline for running a genome-wide association study (GWAS) meta-analysis of chronic fatigue syndrome (CFS) summary statistics. The workflow automates downloading public summary statistics, munging them into a common format, filtering variants, and exporting harmonised datasets that can be combined downstream.

## Repository contents

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

## Data sources

The pipeline retrieves five summary statistics:

| Database                 | Symbol     | cases       | controls    | Trait            | Regression  | Ancestry  | Assembly  | Reference     | Summary Statistics | 
| :----------------------- | :--------- | :---------- | :---------- | :--------------- | :---------- | :-------- | :-------- | :------------ | :----------------- |
| **DecodeME** | DME | 15579 | 259909 | ME/CFS (CCC/IOM)  | Logistic | EUR | GRCh38 | ([Preprint 2025](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-)) | ([GWAS-1](https://osf.io/rgqs3/files/osfstorage))
| **Million Veteran Project** | MVP | 3891 | 443093 | CFS (PheCode 798.1) | Logistic | EUR | GRCh38 | ([Verma 2024](https://pubmed.ncbi.nlm.nih.gov/39024449/)) | ([GCST90479178](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479178/)) |
| **UKBiobank (NealeLab)** | UKBNL | 1659 | 359482 | CFS (self-reported) | Linear | EUR | GRCh37 | ([NealeLab](https://www.nealelab.is/uk-biobank)) | ([20002_1482](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=178908679#gid=178908679))
| **UKBiobank (EIB)**      | UKBEIB| 2092 | 482506 | CFS (self-reported) | Linear | EUR | GRCh37 | ([Dönertaş 2021](https://europepmc.org/article/MED/33959723)) | ([GCST90038694](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038694)) |
| **FinnGen** | FG | 283 | 663029 | Post-viral fatigue | Logistic | FIN | GRCh38 | ([Kurki 2023](https://pubmed.ncbi.nlm.nih.gov/36653562/)) | ([R12_G6_POSTVIRFAT](https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/)) |

Consult the respective project documentation for licensing or data use agreements that may apply to the downloaded data.

## License

This project is released under the terms described in the included `LICENSE` file.
