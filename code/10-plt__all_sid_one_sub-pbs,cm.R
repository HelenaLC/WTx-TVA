# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(pheatmap)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
pbs <- lapply(args[[1]], readRDS)

# wrangling
df <- lapply(pbs, \(x) {
    z <- t(assay(x[rowData(x)$sel, ]))
    data.frame(
        s=factor(x$sid), 
        k=factor(colnames(x)), 
        z, check.names=FALSE)
}) |> bind_rows()

# analysis
fd <- select(df, where(is.numeric))
es <- as.matrix(fd)
es[is.na(es)] <- 0
cm <- cor(t(es))

# plot & save
pheatmap(cm,
    filename=args[[2]], width=15/2.54, height=12/2.54,
    # aesthetics
    cellwidth=2.5, cellheight=2.5, treeheight_row=5, treeheight_col=5,
    fontsize=6, border_color=NA, show_rownames=FALSE, show_colnames=FALSE,
    # legend
    breaks=seq(0, 1, length=n <- 100),
    color=grDevices::colorRampPalette(pals::coolwarm())(n),
    # annotations
    annotation_col=select(df, !where(is.numeric)),
    annotation_colors=list(
        s=setNames(pals::jet(nlevels(df$s)), levels(df$s)), 
        k=setNames(pals::trubetskoy(nlevels(df$k)), levels(df$k))))
