# WTx-CosMx SMI on TVA/CRC

## setup

- this repository is a file size-limit version; any large data are omitted
- let `<x>` denote a slide identifier (run 1: 11,12; run 2: 21,22,23,24)
  - Zarr stores of IF stains should be located at `imgs/<x>/images/`
  - flat files should be located at `data/raw/<x>/`
- in addition, snPATHO-seq un/filtered barcodes should be  
  located at `data/ref/raw` and `data/ref/fil`, respectively
- Gut Cell Atlas reference data should be at `data/gca.rds`
  (retrievable from https://www.gutcellatlas.org/)

## notes

- R version and library have are specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)
- `.Rprofile` is used for handling command line arguments
- `logs/` capture `.Rout` files from `R CMD BATCH` executions
- intermediate results are written to `outs/` (as *.rds*)
- visualizations are written to `plts/` (as *.pdf*)

## steps

**data setup**

- `01-raw.R`
  - read flat files as `SingleCellExperiment`
  - stash non-RNA targets as `altExps`
  - write out as *.h5*-backed object

**quality control**

- `02-fil.R`
  - exclude cells too close to any FOV border
  - exclude cells with low counts, low counts per area,  
  and high negative probe or false code count

**metadata**

- `03-pol.R`
  - read polygon data from *.csv*
  - filter for cells passing QC 
  - write to *.parquet*
- `03-roi.R`
  - retrieve Napari-based  shape annotations
    - histopathological regions (REV/TVA/CRC)
    - invasion-free blood vessels (BV)
    - lymphovascular invarsions (LI)
  - align coordinates with CosMx data
  - stash annotations as cell metadata

**clustering**  

- `03-ist.R`:  
`InSituType` clustering supervised by reference profiles  
extracted from (pooled) snPATHO-seq data on adjacent sections

- `04-sub.R`: subset cells into epi(thelia), imm(unne), (str)mal

- `05-jst.R`:  
`InSituType` subclustering (analogous to above),  
using distinct references profiles for each subset
  
- `00-lab.R`: used to relabel  
`03-ist.R` outputs (automated; using `meta/lab/lv1.json`), and  
`05-jst.R` outputs (manual; using `meta/lab/lv2,<sub>.json`)
  
- downstream analyses
    - `03-pro.R` and `05-rep.R`: PCA using feature subset  
    according to `03-ist.R` and `05-jst.R`, respectively
    - `03-sig.R`: gene set signature scoring
    - `04-ccc.R`: cell cell-communication analysis
    - `06-trj.R`: epithelial trajectory inference
    - `06-ctx.R`: spatial context/niche analysis
    - `06-cty.R`: niches blinded to epithelia

- `10-plt__<by1>__<by2>-<out1>,<out2>,<plt>.R`
  - pool outputs `<out1/2>` according to `<by1/2>`  
  - generates `plts/<out1>,<out2>,<plt>.pdf`  
  (depending on `<by1/2>`, name may also include  
  subset (`sub`) or section (`sid`) identifier)
    - `sid` = one section
    - `all_sid` = all sections
    - `all_sip` = also include LN (241)
    - `all_sid_all_sub` = all sections and subsets
    