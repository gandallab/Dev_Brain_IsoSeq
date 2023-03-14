<<<<<<< HEAD
# devBrainIsoSeq
=======
# Fetal_bulk_IsoSeq

<https://gandallab.github.io/devBrainIsoSeq/>

## How to use this repository

All ipynb, Rmd, and qmd notebooks should run using data files
included in the repository or easily accessible elsewhere.

Paths are always relative to the root of the project directory, e.g.
`analysis/BulkTxomeAnalysis.qmd`, `data/cp_vz_0.75_min_7_recovery_talon.gtf.gz`.

Render the Quarto site (or individual qmd files) locally with
`quarto render`. html output will appear in `_site`.

Publish the site with `quarto publish gh-pages`. Changes will
be pushed to GitHub on the gh-pages branch.

The entire site must be published at once. You can use
`quarto publish gh-pages --no-render` to only update files that have been
manually rendered first.

See `_quarto.yml` for details.

To add new analysis notebooks to the website, please add their path to
`_quarto.yml` under both `render:` and `contents:`. WARNING: Even when trying
to render a single file, `quarto render [path]` will fail unless the path is
listed under `render:` in `_quarto.yml` due to a glitch.

## Directory structure

- analysis - fully reproducible notebooks that will be knitted into the site
- code - other code that may take more effort to run/won't generate output
  for site
- data - unprocessed data files specific to the project
- ouput - processed data files and other outputs generated from analysis or
  code
- local - scratch directory specific to the user's local environment
- ref - unprocessed data files that come from an external source (if files
  are too large to fit in the GitHub, there will be clear directions to obtain
  them in a README)
- workflow - Snakemake workflow files
- _quarto.yml - Quarto config
- styles.css - stylesheet for the Quarto website because there's nowhere else
  obvious to put it

## Hoffman2 directories

Rough mirror of this repository + Snakemake workflow results:

```
/u/project/gandalm/jops/ucdavis_cpvz_isoseq
```

Raw reads as downloaded from UC Davis:

```
/u/project/gandalm/shared/isoSeq/UCDavis/
```
>>>>>>> 8f3d6ce (add code, data, analyses for Figures 1 to 3)
