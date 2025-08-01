# Code to: Tracing colorectal <br/> malignancy transformation <br/> from cell to tissue scale

Helena L. Crowell, Irene Ruano, Zedong Hu, Yourae Hong, Gin Caratù, Hubert Piessevaux,  
Ashley Heck, Rachel Liu, Max Walter, Megan Vandenberg, Kim Young, Dan McGuire,  
Evelyn Metzger, Margaret L. Hoang4, Joseph M. Beechem, Sabine Tejpar,  
Anna Pascual-Reguant* & Holger Heyn*

[PREPRINT](https://doi.org/10.1101/2025.06.23.660674)

## setup

- this repository is file size-limit; any data beyond KBs are omitted
- let `<x>` denote a slide identifier (run 1: 11,12; run 2: 21,22,23,24)
  - Zarr stores of IF stains should be at `imgs/<x>/images/`
  - corresponding flat files should be at `data/raw/<x>/`
- in addition, snPATHO-seq un- and filtered barcodes  
  should be at `data/ref/raw` and `-/fil`, respectively
- Gut Cell Atlas reference data should be at `data/gca.rds`  
  (it can be retrieved from https://www.gutcellatlas.org/)

## notes

- Software  versions used throughout this study are  
  captured in the session information [here](inf.txt)
- A complete list of package dependencies (as well as  
  installation commands) are provided in `code/09-inf.R`
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
  
**downstream analyses**

- `03-pro.R` and `05-rep.R`: PCA using feature subset  
  according to `03-ist.R` and `05-jst.R`, respectively
- `03-sig.R`: gene set signature scoring
- `04-ccc.R`: cell cell-communication analysis
- `06-trj.R`: epithelial trajectory inference
- `06-ctx.R`: spatial context/niche analysis
- `06-cty.R`: niches blinded to epithelia

**visualization**

- majority of scripts serve the purpose of collecting  
  specific subsets of results, and visualizing them
- `10-plt__<by1>__<by2>-<out1>,<out2>,<plt>.R`
  - pool outputs `<out1/2>` according to `<by1/2>`  
  - generates `plts/<out1>,<out2>,<plt>.pdf`  
  (depending on `<by1/2>`, name may also include  
  subset (`sub`) or section (`sid`) identifier)
    - `sid` = one section
    - `all_sid` = all colon sections
    - `all_sip` = also include LN (241)
    - `all_sid_all_sub` = all sections and subsets
    