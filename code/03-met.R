# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpmisc)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])

# wrangling
df <- data.frame(colData(sce))
xs <- grep("^Mean", names(df), value=TRUE)
df <- mutate(df, across(all_of(xs), ~asinh(./200)))

# PanCK
d <- density(x <- df$Mean.PanCK)
y <- d$x[find_peaks(-d$y)]
(th_pck <- y+diff(range(x))*0.05)

# CD45
d <- density(x <- df$Mean.CD45)
y <- d$x[find_peaks(-d$y)]
(th_cd45 <- y-diff(range(x))*0.05)

# filtering
table(epi <- 
    df$Mean.PanCK > th_pck & 
    df$Mean.CD45 < th_cd45)

# saving
kid <- rep("epi", sum(epi))
names(kid) <- colnames(sce)[epi]
ist <- list(clust=kid)
saveRDS(ist, args[[2]])
