# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(AUCell)
    library(scuttle)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
.sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)
auc <- mapply(
    SIMPLIFY=FALSE, x=args[[1]],
    y=split(args[[2]], .sid(args[[2]])), \(x, y) {
        sce <- readRDS(x)
        ist <- lapply(y, readRDS)
        kid <- lapply(ist, \(.) .$clust)
        sub <- rep.int(.sub(y), sapply(kid, length))
        kid <- unlist(kid); names(sub) <- names(kid)
        idx <- intersect(colnames(sce), names(kid))
        sid <- .sid(x); sub <- sub[idx]; kid <- kid[idx]
        out <- SummarizedExperiment(
            list(counts=assay(sce[, idx])),
            colData=data.frame(sid, sub, kid))
        cs <- split(seq_len(ncol(out)), out$sub)
        names(sub) <- sub <- names(cs)
        lapply(sub, \(.) {
            rd <- rowData(sce)$sub
            gs <- sapply(rd, match, x=.)
            gs <- !is.na(gs)
            out[gs, cs[[.]]]
        }) 
    }) 

# wrangling
names(sub) <- sub <- c("epi", "imm", "str")
auc <- lapply(sub, \(sub) do.call(cbind, lapply(auc, \(.) .[[sub]])))
#auc <- lapply(auc, \(.) { assay(.) <- scale(assay(.)); . })

# plotting
.a <- \(.) gsub("^epi\\.", "", .)
.b <- \(.) {
    . <- gsub("HALLMARK_", "H_", .)
    . <- gsub("GAVISH_3CA_MALIGNANT_METAPROGRAM_[0-9]+_", "MP_", .)
    .
}
.p <- \(x, i, xo=TRUE, yo=TRUE) {
    n <- length(unique(x$kid))
    x <- select(x, kid, name, value)
    y <- pivot_wider(x, names_from="kid")
    z <- as.matrix(y[, -1]); rownames(z) <- y[[1]]
    ggplot(x, aes(kid, name, fill=value)) +
        (if (xo) scale_x_discrete(limits=.yo(z), position="top", labels=.a)) +
        (if (yo) scale_y_discrete(limits=.xo(z), labels=.b)) +
        scale_fill_gradient2(
            "z-scaled\nmean AUCell",
            low="turquoise", high="purple",
            limits=c(-2.5, 2.5), n.breaks=5) +
        ggtitle(bquote(bold(.(i))~"(N ="~.(n)*")")) +
        geom_tile() + coord_equal(3/4, expand=FALSE) +
        .thm_fig_c("bw") + theme(
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.text.y=element_text(size=3),
            axis.text.x=element_text(size=4, angle=45, vjust=0.5, hjust=0))
}

# joint
ps <- lapply(auc, \(se) {
    ex <- c("goblet", "EE", "entero")
    ex <- paste0("epi.", ex)
    sf <- se[, !se$kid %in% ex]
    mu <- aggregateAcrossCells(sf, 
        colData(sf)[ids <- c("kid")],
        statistics="mean")
    df <- t(assay(mu)) |>
        data.frame(colData(mu)[ids]) |>
        pivot_longer(all_of(rownames(mu))) |>
        group_by(name) |>
        mutate(value=.z(value)) |>
        ungroup()
    .p(df, se$sub[1])
})

# split
qs <- lapply(sub, \(sub) {
    se <- auc[[sub]]
    mu <- aggregateAcrossCells(se,
        colData(se)[ids <- c("sid", "kid")],
        statistics="mean")
    df <- t(assay(mu)) |>
        data.frame(colData(mu)[ids]) |>
        pivot_longer(all_of(rownames(mu))) |>
        group_by(sid, name) |>
        mutate(value=.z(value)) |>
        ungroup()
    by(df, df$sid, \(fd) .p(fd, paste0(sub, "-", fd$sid[1])))
}) |> Reduce(f=c)

# saving
pdf(args[[3]], width=15/2.54, height=12/2.54, onefile=TRUE)
for (p in c(ps, qs)) print(p); dev.off()

# saving
gg <- c(ps, qs)
tf <- replicate(length(gg), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(gg)) {
    df <- gg[[.]]$data
    nx <- length(unique(df$kid))/6+max(nchar(paste(df$name)))/12
    ny <- length(unique(df$name))/8+max(nchar(paste(df$kid)))/9/2
    pdf(tf[[.]], 
        width=(2+nx)/2.54, 
        height=(0.5+ny)/2.54)
    print(gg[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
