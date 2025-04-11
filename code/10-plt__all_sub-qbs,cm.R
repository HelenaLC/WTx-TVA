# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(pheatmap)
    library(SingleCellExperiment)
})

# loading
args=list(list.files("outs", "qbs", full.names = TRUE), "plts/qbs,cm.pdf")
pbs <- lapply(args[[1]], readRDS)

# wrangling
df <- lapply(pbs, \(x) {
    x <- x[, x$sid != "all"]
    z <- t(assay(x[rowData(x)$sel, ]))
    df <- data.frame(
        s=factor(x$sid), 
        k=factor(x$kid), 
        z, check.names=FALSE)
})

# plot & save
pdf(args[[2]], onefile=TRUE, width=15/2.54, height=12/2.54)
ps <- lapply(df, \(fd) {
    df_es <- select(fd, where(is.numeric))
    df_id <- select(fd, !where(is.numeric))
    pal_s <- setNames(pals::jet(length(.)), . <- levels(fd$s))
    pal_k <- setNames(pals::trubetskoy(length(.)), . <- levels(fd$k))
    es <- as.matrix(df_es)
    es[is.na(es)] <- 0
    cm <- cor(t(es))
    nm <- rownames(df_id)
    dimnames(cm) <- list(nm, nm)
    pheatmap(cm,
        # legend
        breaks=seq(0, 1, length=n <- 99),
        color=rev(pals::brewer.rdylbu(n)),
        # annotations
        annotation_col=df_id, 
        annotation_colors=list(s=pal_s, k=pal_k),
        # aesthetics
        cellwidth=2, cellheight=2, treeheight_row=8, treeheight_col=8,
        fontsize=4, border_color=NA, show_rownames=FALSE, show_colnames=FALSE)
})
dev.off()
