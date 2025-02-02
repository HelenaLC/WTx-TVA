# WTx-CosMx SMI on TVA/CRC

## setup

- R version and library have are specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)
- `.Rprofile` is used for handling command line arguments
- `logs/` capture `.Rout` files from `R CMD BATCH` executions
- intermediate results are written to `outs/` (as *.rds*)
- visualizations are written to `plts/` (as *.pdf*)

## steps

- `01-raw.R`
  - read flat files as `SingleCellExperiment`
  - stash non-RNA targets as `altExps`
  - write out as *.h5*-backed object
  
- `02-fil.R`
  - exclude cells too close to any FOV border
  - exclude cells with low counts, low counts per area,  
  and high negative probe or false code count
  
- `03-ist.R`
  - `InSituType` clustering supervised by reference profiles  
  extracted from snPATHO-seq data on adjacent sections
  
- `04-sub.R`
  - split and subset cells into  
  epi(thelia), imm(unne), (str)mal
  
- `05-jst.R`
  - `InSituType` subclustering (analogous to above);  
  using distinct references profiles for each subset