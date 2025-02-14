# color palette ----
.cmy <- c("cyan2", "magenta2", "gold")
.pal <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# prettified title ----
.lab <- \(id, n) bquote(bold(.(id))~"(N ="~.(format(n, big.mark=","))*")")

# hierarchical cluster -----
.hc <- \(x, dim=c("x", "y")) {
    dim <- match.arg(dim)
    if (!is.matrix(x)) {
        y <- as.matrix(x[, -1])
        rownames(y) <- x[[1]]
        x <- y
    }
    switch(dim, 
        x=rownames(x)[hclust(dist(x))$order],
        y=colnames(x)[hclust(dist(t(x)))$order])
}

# thresholded z-normalization ----
.z <- \(x, th=2.5) {
    if (is.null(dim(x))) {
        x[x < 0] <- 0
        sd <- sd(x, na.rm=TRUE)
        x <- x-mean(x, na.rm=TRUE)
        if (sd != 0) x <- x/sd
    } else {
        mus <- colMeans(x)
        sds <- colSds(x)
        x <- sweep(x, 2, mus, `-`)
        x <- sweep(x, 2, sds, `/`)
    }
    x[x > +th] <- +th
    x[x < -th] <- -th
    return(x)
}

# compositional barplot ----

.plt_fq <- \(sce, x, y, id="", pal=.pal) {
    suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(SummarizedExperiment)
    })
    df <- data.frame(colData(sce))
    fd <- tally(group_by(df, !!sym(y), !!sym(x)))
    ggplot(fd, aes(n, fd[[y]], fill=fd[[x]])) +
        geom_col(position="fill", key_glyph="point", width=1, col="black", linewidth=0.1) +
        guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, size=2))) +
        scale_fill_manual(x, values=pal) +
        scale_x_continuous(n.breaks=2) +
        coord_cartesian(expand=FALSE) +
        ggtitle(.lab(id, nrow(df))) +
        theme_linedraw(6) + theme(
            aspect.ratio=2/3,
            axis.title=element_blank(),
            legend.key.size=unit(0.2, "lines"))
}

# dimensionality reduction ----

.plt_dr <- \(sce, ids, dr="UMAP") {
    suppressPackageStartupMessages({
        library(ggplot2)
        library(patchwork)
        library(SingleCellExperiment)
    })
    dr <- reducedDim(sce, dr)
    colnames(dr) <- c("x", "y")
    df <- data.frame(dr, colData(sce)[ids])
    ps <- lapply(ids, \(.) {
        n <- length(unique(df[[.]]))
        ggplot(df, aes(x, y, col=.data[[.]])) +
        scale_color_manual(values=if (n == 3) .cmy else .pal)
    })
    wrap_plots(ps) +
        plot_layout(nrow=1) &
        geom_point(shape=16, stroke=0, size=0.2) &
        guides(col=guide_legend(override.aes=list(size=2, alpha=1))) &
        theme_void() & theme(aspect.ratio=1, legend.key.size=unit(0, "lines"))
}

# principal component regression ----

.pcr <- \(sce, ids) {
    suppressPackageStartupMessages({
        library(dplyr)
        library(SummarizedExperiment)
    })
    pcs <- reducedDim(sce, "PCA")
    pcr <- lapply(ids, \(id) {
        fit <- summary(lm(pcs ~ sce[[id]]))
        r2 <- sapply(fit, \(.) .$adj.r.squared)
        data.frame(id, pc=seq_along(r2), r2)
    }) |> 
        do.call(what=rbind) |>
        mutate(id=factor(id, ids))
}

.plt_pcr <- \(df, n=c(1, 30), x="pc", y="r2", c="id") {
    suppressPackageStartupMessages(library(ggplot2))
    df[[y]][df[[y]] < 0] <- NA
    ggplot(df, 
        aes(.data[[x]], .data[[y]], col=.data$id)) +
        geom_point(na.rm=TRUE) +
        geom_line(show.legend=FALSE) +
        coord_cartesian(xlim=n, ylim=c(0, 1)) +
        scale_x_continuous(breaks=c(1, seq(5, max(n), 5))) +
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
        guides(col=guide_legend("predictor", override.aes=list(size=2))) +
        labs(x="principal component", y="coeff. of determination") +
        theme_minimal() + theme(
            panel.grid.minor=element_blank(),
            legend.key.size=unit(0.5, "lines"))
}

# signatures ----

.plt_hm <- \(sce, mtx, by, id="") {
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(AUCell)
        library(ggplot2)
        library(SingleCellExperiment)
    })
    if (is(mtx, "SummarizedExperiment")) mtx <- assay(mtx)
    if (!is.matrix(mtx)) mtx <- as.matrix(mtx)
    cd <- colData(sce)[by]
    df <- data.frame(t(mtx), cd, check.names=FALSE)
    fd <- group_by(df, across(all_of(by)))
    for (. in seq_along(by))
        fd <- summarize(fd, .groups="drop_last",
            across(all_of(rownames(mtx)), mean))
    fd <- mutate(fd, across(rownames(mtx), .z))
    ggplot(
        pivot_longer(fd, -!!sym(by[1])), 
        aes(.data[[by[1]]], name, fill=value)) +
        geom_tile() + 
        coord_equal(expand=FALSE) + 
        ggtitle(.lab(id, nrow(cd))) +
        scale_x_discrete(limits=.hc(fd, "x")) +
        scale_y_discrete(limits=.hc(fd, "y")) +
        labs(fill="z-scaled\nsig. score") +
        scale_fill_gradientn(
            colors=rev(hcl.colors(9, "BrBG")),
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) +
        theme_bw(9) + theme(
            axis.title=element_blank(), 
            plot.title=element_text(hjust=0.5), 
            legend.key.width=unit(0.4, "lines"),
            legend.key.height=unit(0.8, "lines"))
}
