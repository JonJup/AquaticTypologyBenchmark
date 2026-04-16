**Pan-European freshwater ecotypologies and ecological benchmarks from data-driven typologies and joint species distribution models**

<!-- Optional badges -->
<!-- [![DOI](https://zenodo.org/badge/DOI/xx.xxxx/zenodo.xxxxxxx.svg)](https://doi.org/xx.xxxx/zenodo.xxxxxxx) -->
<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) -->
<!-- [![Shiny app](https://img.shields.io/badge/Shiny-app-blue)](https://<your-shinyapps-url>) -->

This repository accompanies the manuscript:    

Jupke et al (*in preparation*) Grading on a Curve: Simulation-Based Benchmarks for the Biological Evaluation of Freshwater Typology Performance


<!-- > **[Full title]** (202x). [Authors]. *[Journal]*. DOI: [xx.xxxx/xxxxxx]. Preprint: [link]. --> 

It contains the code, configuration, and (where licensing permits) derived data needed to reproduce the ecotypologies, HMSC joint species distribution models, QRF benchmarks, and fuzzy-clustering analyses reported in the paper, together with the Shiny application used to explore the results.

---

## Overview

PULSE develops **pan-European ecological benchmarks** for freshwater biomonitoring by combining:

- **Data-driven ecotypologies** built from catchment- and reach-scale environmental descriptors, using fuzzy spatial clustering (SKATER-CON).
- **Joint species distribution models** (HMSC) fitted separately for diatoms, fish, macroinvertebrates, and macrophytes, used to partition community variation into environmental, spatial, and residual components and to generate counterfactual simulations.
- **Context-dependent benchmarks** derived from Quantile Regression Forests (QRF), producing reference-condition expectations conditioned on the typological and environmental context of each site.
- **Typology evaluation** against crisp and fuzzy alternatives using ANOSIM, PERMANOVA, classification strength, and coherence diagnostics.

The pipeline is deliberately modular so that each component — typology construction, HMSC fitting, QRF benchmarking, and evaluation — can be re-run independently against updated inputs.

## Scope of the paper

- **Taxonomic groups.** Diatoms, fish, macroinvertebrates, macrophytes.
- **Geographic scope.** Continental Europe.
- **Dataset.** ~400,000 biological samples harmonised from ~90 national and regional datasets, paired with catchment-derived environmental descriptors.
- **Key outputs.** Ecotypology assignments, HMSC posterior summaries, QRF benchmark surfaces, and a set of reproducible figures and tables.

## Repository structure

```
.
├── R/                    # Reusable functions sourced by the pipeline
├── scripts/              # Top-level, numbered scripts (one per pipeline stage)
│   ├── 01_prepare_data.R
│   ├── 02_build_typology.R
│   ├── 03_fit_hmsc.R
│   ├── 04_simulate_counterfactuals.R
│   ├── 05_fit_qrf_benchmarks.R
│   └── 06_evaluate.R
├── hpc/                  # SLURM / Singularity templates for HMSC fitting
├── config/               # YAML configuration for each taxonomic group
├── shiny/                # PULSE Shiny application
├── data/
│   ├── raw/              # (Gitignored) source data — see Data section
│   └── processed/        # Harmonised inputs produced by 01_prepare_data.R
├── outputs/              # Fitted models, predictions, figures, tables
├── renv.lock             # Pinned R package versions (renv)
├── DESCRIPTION           # Optional — if the repo is structured as an R package
└── README.md
```

*(Adjust the tree to match the actual layout before publishing.)*

## Data

Raw biological and environmental data are **not redistributed** in this repository because several source datasets are governed by institutional or national licences. Each script in `scripts/` documents the expected input files and their provenance; the harmonisation code in `01_prepare_data.R` shows exactly how the processed inputs are derived.

Where licensing allows, harmonised inputs and fitted model objects are archived on Zenodo:

- **Processed inputs:** [DOI placeholder]
- **Fitted HMSC models and QRF benchmarks:** [DOI placeholder]

The Shiny app loads fitted models remotely from the Zenodo archive by default, so end users do not need to refit the models to explore the benchmarks.

## Installation

The analyses were run with **R ≥ 4.3**. Package versions are pinned with `renv`.

```r
# Clone the repository, then from the project root:
install.packages("renv")
renv::restore()
```

System-level dependencies:

- GDAL, GEOS, PROJ (for the spatial pipeline)
- A C++17 toolchain (for `Hmsc` and `randomForest`)
- Optional: Singularity/Apptainer for reproducing the HPC runs

## Reproducing the analyses

The full pipeline can be reproduced by running the numbered scripts in order. Intermediate outputs are cached under `outputs/`, so individual stages can be re-run without repeating the whole workflow.

```bash
Rscript scripts/01_prepare_data.R        config/diatoms.yml
Rscript scripts/02_build_typology.R      config/diatoms.yml
Rscript scripts/03_fit_hmsc.R            config/diatoms.yml
Rscript scripts/04_simulate_counterfactuals.R config/diatoms.yml
Rscript scripts/05_fit_qrf_benchmarks.R  config/diatoms.yml
Rscript scripts/06_evaluate.R            config/diatoms.yml
```

Repeat with the configuration file for each taxonomic group.

### HPC workflow

HMSC fitting for the full dataset is computationally intensive. Templates in `hpc/` are provided for SLURM clusters running Singularity. A typical invocation:

```bash
sbatch hpc/fit_hmsc.slurm config/macroinvertebrates.yml
```

The template exposes knobs for the number of chains, thinning, transient period, and memory per task, matching the settings reported in the manuscript.

## Shiny application

The PULSE Shiny app provides an interactive view of the benchmarks:

- QRF prediction and density visualisation for a user-supplied site or catchment
- Batch prediction from uploaded CSV/Excel tables
- KNN-based imputation of missing environmental predictors
- Diagnostic overlays against the reference distribution of each ecotype

Run locally with:

```r
shiny::runApp("shiny")
```

A hosted version is available at: **[URL placeholder]**.

## Related software

- **SKATER-CON** — R package for fuzzy spatial clustering via coassociation matrices, used to build the ecotypology: [link placeholder].
- **DiaThor** — diatom ecological indicator package, used for cross-checking diatom-based inference: [link placeholder].

## Citation

If you use PULSE, please cite both the paper and the archived code release:

```bibtex
@article{pulse_paper_202x,
  title   = {[Full title]},
  author  = {[Authors]},
  journal = {[Journal]},
  year    = {202x},
  doi     = {xx.xxxx/xxxxxx}
}

@software{pulse_code_202x,
  title   = {PULSE: code and fitted models},
  author  = {[Authors]},
  year    = {202x},
  doi     = {xx.xxxx/zenodo.xxxxxxx},
  version = {vX.Y.Z}
}
```

## Funding and acknowledgments

This work was supported by the **DFG Walter Benjamin Fellowship** (grant no. [xxx]) and benefited from collaboration with **JRC ECOSTAT** and the wider European freshwater monitoring community. We thank the data providers listed in the manuscript supplement for sharing harmonised biological and environmental records.

## License

- **Code:** [MIT / GPL-3.0 / Apache-2.0 — choose one] — see `LICENSE`.
- **Derived data and figures:** CC BY 4.0, unless otherwise stated in the relevant Zenodo record.
- Source datasets retain their original licences; contact the listed providers for redistribution terms.

## Contact

Questions, bug reports, and collaboration enquiries are welcome via GitHub issues, or by email to **[corresponding author]**.
