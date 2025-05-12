# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scuttle)
    library(SingleCellExperiment)
})

.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
.sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)
.foo <- \(., by) split(., get(paste0(".", by))(.))
se <- mapply(
    x=.foo(args[[1]], "sid"),
    y=.foo(args[[2]], "sid"),
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        lys <- readRDS(x)
        ist <- lapply(.foo(y, "sub"), readRDS)
        kid <- lapply(ist, \(.) data.frame(kid=.$clust))
        kid <- bind_rows(kid, .id="sub")
        # construct 'SummarizedExperiment'
        cd <- DataFrame(kid, sid=.sid(x))
        cs <- intersect(rownames(lys[[1]]), rownames(cd))
        es <- lapply(lys, \(.)
            `rownames<-`(mx <- t(.[cs, ]),
                gsub("^(s|r)-", "", rownames(mx))))
        # characterize sender/receiver as
        # total, ligand-receptor, pathway
        rd <- data.frame(
            row.names=nm <- rownames(es[[1]]),
            typ=ifelse(grepl("total", nm), "tl",
                ifelse(grepl("-", nm), "lr", "pw")))
        se <- SummarizedExperiment(es, rowData=rd, colData=cd[cs, ])
    }) |> do.call(what=cbind)

# average by cluster
id <- colData(se)["kid"]
names(sr) <- sr <- c("s", "r")
mu <- lapply(sr, \(.) {
    sf <- summarizeAssayByGroup(se, id, assay.type=., statistics="mean")
    df <- data.frame(t(assay(sf)), colData(sf), check.names=FALSE)
    fd <- pivot_longer(df, all_of(rownames(sf)), names_to="lr")
    fd |>
        mutate(typ=rowData(se)[fd$lr, "typ"]) |>
        mutate(sub=colData(se)[match(fd$kid, se$kid), "sub"])
}) |> bind_rows(.id="sr")

.hm <- \(df, k1, k2, n1, n2, ht=FALSE) {
    df <- mutate(df, i=row_number())
    # exclude homotypic
    fd <- filter(df, kid %in% k1)
    idx <- if (ht) fd else 
        fd |>
        group_by(lr, sr) |>
        slice_max(value, n=1) |>
        group_by(lr) |>
        mutate(n=length(unique(kid))) |>
        filter(n > 1)
    # get top-N sender/receiver
    one <- fd |>
        filter(i %in% idx$i) |>
        group_by(kid) |>
        slice_max(value, n=n1)
    # get top-M receiver/sender
    two <- mapply(
        .sr=one$sr, .lr=one$lr, 
        SIMPLIFY=FALSE, \(.sr, .lr) {
            df |>
                filter(kid %in% k2) |>
                filter(lr == .lr) |>
                filter(sr == setdiff(c("s", "r"), .sr)) |>
                slice_max(value, n=n2)
        }) |> do.call(what=rbind)
    gg <- df |>
        filter(lr %in% one$lr) |>
        filter(kid %in% rbind(one, two)$kid) |>
        filter(
            (sr == "s" & kid %in% k1) |
            (sr == "r" & kid %in% k2)) |>
        group_by(sr, lr) |> mutate_at("value", .z) |> ungroup() |>
        mutate_at("sr", factor, c("s", "r"), c("sender", "receiver"))
    mx <- if (ht) mx <- filter(gg, sr == "sender") else gg
    mx <- mx |>
        select(c(kid, lr, value)) |>
        pivot_wider(names_from="lr")
    my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
    asp <- length(unique(gg$kid))/length(unique(gg$lr))
    ggplot(gg, aes(kid, lr, fill=value)) + geom_tile() +
        facet_grid(~sr, space="free", scales="free_x") +
        scale_fill_gradient2(
            "z-scaled\nmean CCC",
            limits=c(-2.5, 2.5), n.breaks=6,
            low="cadetblue", mid="ivory", high="firebrick") +
        scale_x_discrete(limits=\(.) intersect(.xo(my), .)) +
        scale_y_discrete(limits=.yo(my)) +
        .thm_fig_c("minimal") + theme(
            legend.position="none",
            panel.grid=element_blank(),
            axis.title=element_blank(),
            axis.text.y=element_text(size=2),
            axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5))
}

# plotting
ij <- unique(mu$sub)
ij <- expand.grid(ij, ij)
ps <- mapply(
    i=ij[, 1], j=ij[, 2], 
    SIMPLIFY=FALSE, \(i, j) {
        k1 <- unique(mu$kid[mu$sub == i])
        k2 <- unique(mu$kid[mu$sub == j])
        hm <- .hm(filter(mu, typ == "lr"), k1, k2, 5, 1, i == j)
        hm + ggtitle(paste(i, j, sep="-"))
    })

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    n <- length(unique(unlist(ij[., ])))
    n <- ifelse(n == 1, 0.5, 1)
    x <- ps[[.]]$data
    h <- length(unique(x$lr))/12+max(nchar(x$kid))/12
    w <- length(unique(x$kid))/8/n+max(nchar(x$lr))/12
    pdf(tf[[.]], 
        width=(w)/2.54, 
        height=(1+h)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
