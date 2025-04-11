# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
df <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        sce <- readRDS(x)
        ist <- readRDS(y)$clust
        sub <- rowData(sce)$sub
        ns <- sapply(sub, length)
        gs <- rep.int(names(sub), ns)
        gs <- split(gs, unlist(sub))
        cs <- split(names(ist), ist)
        # analysis
        lapply(names(gs), \(sub) {
            se <- sce[gs[[sub]], cs[[sub]]]
            cm <- cor(t(assay(se))); diag(cm) <- NA
            df <- data.frame(i=rownames(cm), cm)
            fd <- pivot_longer(df, -i, names_to="j")
            data.frame(sid=sce$sid[1], sub, fd)
        }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

# plotting
ps <- by(df, list(df$sid, df$sub), \(fd) {
    mx <- pivot_wider(fd, id_cols="i", names_from="j")
    my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
    ggplot(fd, aes(i, j, fill=value)) + 
        scale_x_discrete(limits=.xo(my)) +
        scale_y_discrete(limits=.xo(my)) +
        scale_fill_gradient2(
            "corr.", breaks=seq(-1, 1, 0.2), 
            low="blue", high="red", na.value="grey") +
        ggtitle(paste0(fd$sid[1], "-", fd$sub[1])) +
        geom_tile() + coord_equal(expand=FALSE) +
        .thm_fig_c("minimal") + theme(
            axis.title=element_blank(),
            axis.text=element_text(size=3),
            axis.text.x=element_text(angle=45, hjust=1))
})

# saving
pdf(args[[3]], onefile=TRUE, width=10/2.54, height=8/2.54)
for (p in ps) print(p); dev.off()
