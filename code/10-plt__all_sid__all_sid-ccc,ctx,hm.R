# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scuttle)
    library(SingleCellExperiment)
})

se <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        ccc <- readRDS(x)
        ctx <- readRDS(y)
        # construct 'SummarizedExperiment'
        cs <- intersect(rownames(ccc[[1]]), ctx$cid)
        cd <- ctx[match(cs, ctx$cid), ]; row.names(cd) <- cs
        es <- lapply(ccc, \(.)
            `rownames<-`(mx <- t(.[cs, ]),
            gsub("^(s|r)-", "", rownames(mx))))
        # characterize sender/receiver as
        # total, ligand-receptor, pathway
        rd <- data.frame(
            row.names=nm <- rownames(es[[1]]),
            typ=ifelse(grepl("total", nm), "tl",
                ifelse(grepl("-", nm), "lr", "pw")))
        se <- SummarizedExperiment(es, rowData=rd, colData=cd)
    }) |> do.call(what=cbind)

# average by context
id <- colData(se)["ctx"]
names(sr) <- sr <- c("s", "r")
mu <- lapply(sr, \(.) {
    sf <- summarizeAssayByGroup(se, id, assay.type=., statistics="mean")
    df <- data.frame(t(assay(sf)), colData(sf), check.names=FALSE)
    fd <- pivot_longer(df, all_of(rownames(sf)), names_to="lr")
    fd |> mutate(typ=rowData(se)[fd$lr, "typ"])
}) |> bind_rows(.id="sr") |>
    # average sender/receiver
    group_by(ctx, typ, lr) |>
    summarise_at("value", mean) |>
    ungroup()

# add inter-context FCs
cs <- levels(mu$ctx)
mx <- pivot_wider(select(mu, -typ), names_from="lr")
my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
fc <- sapply(rownames(my), \(.) {
    mz <- my[-match(., rownames(my)), ]
    my[., ]/colMeans(mz)
})
fc <- data.frame(lr=rownames(fc), fc, row.names=NULL)
fc <- pivot_longer(fc, -lr, names_to="ctx", values_to="fc")
mu <- left_join(mu, fc, by=c("lr", "ctx"))
mu <- mutate(mu, ctx=factor(ctx, unique(ctx)))

ps <- lapply(c("pw", "lr", cs), \(id) {
    if (id == "pw") {
        n <- 10 
        mu <- mu |> 
            filter(typ == "pw") |> 
            select(-typ)
    } else if (id == "lr") { 
        n <- 10
        mu <- mu |> 
            filter(typ == "lr") |> 
            select(-typ)
    } else { 
        n <- 30 
        cs <- id
        mu <- mu |> 
            filter(typ == "lr") |> 
            select(-typ)
    }
    # selection
    nm <- mu |>
        filter(ctx %in% cs) |>
        group_by(ctx) |>
        slice_max(fc, n=n) |>
        pull(lr)
    mv <- mu |>
        select(-fc) |>
        filter(lr %in% nm) |>
        group_by(lr) |>
        mutate_at("value", .z)
    # hierarchical clustering
    mx <- pivot_wider(mv, names_from="lr")
    my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
    # plotting
    ggplot(mv, aes(ctx, lr, fill=value)) + 
        geom_tile() + ggtitle(id) +
        scale_fill_gradient2(
            "z-scaled\nmean CCC",
            limits=c(-2.5, 2.5), n.breaks=6,
            low="cadetblue", mid="ivory", high="firebrick") +
        scale_x_discrete(limits=\(.) intersect(.xo(my), .)) +
        scale_y_discrete(limits=.yo(my)) +
        coord_equal(2/3, expand=FALSE) +
        .thm_fig_c("minimal") + theme(
            legend.position="none",
            panel.grid=element_blank(),
            axis.title=element_blank(),
            axis.text.y=element_text(size=2),
            axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5)) 
})

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    x <- ps[[.]]$data
    w <- nlevels(x$ctx)/8+max(nchar(x$lr))/12
    h <- length(unique(x$lr))/12+max(nchar(levels(x$ctx)))/12
    pdf(tf[[.]], 
        width=(w)/2.54, 
        height=(0.5+h)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
