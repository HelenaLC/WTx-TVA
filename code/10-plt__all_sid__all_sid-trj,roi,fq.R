# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tidytext)
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# args <- list(
#     list.files("outs", "trj", full.names=TRUE),
#     list.files("outs", "roi", full.names=TRUE),
#     "plts/trj,roi,fq.pdf")
# args <- lapply(args, \(.) grep("241", ., value=TRUE, invert=TRUE))

df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        trj <- readRDS(x)
        sce <- readRDS(y)
        # wrangling
        pt <- .q(averagePseudotime(trj$slingshot))
        cs <- colnames(sce)[!is.na(sce$roi) & !grepl("(BV|LI)", sce$roi)]
        df <- data.frame(colData(sce)[cs, ], pt=pt[cs])
    }) |> do.call(what=rbind)

# average by binned pseudotime
xs <- seq(-(dx <- 0.005), 1.005, 0.01)
mu <- df |>
    filter(!is.na(pt), !is.na(typ)) |>
    mutate(tyq=gsub(".*(REF|TVA|CRC).*", "\\1", typ)) |>
    mutate(tyq=factor(tyq, rev(names(.pal_roi)))) |>
    mutate(
        t=cut(pt, breaks=xs),
        t=xs[as.integer(t)]+dx,
        t=factor(t, paste((xs+dx)[-length(xs)]))) |>
    group_by(sid, t, tyq) |>
    dplyr::count(tyq) |>
    group_by(sid, t) |>
    mutate(p=n/sum(n))

# plotting
ps <- by(mu, mu$sid, \(mv) {
    ggplot(mv, aes(t, p, fill=tyq)) +
        geom_col(width=1, key_glyph="point", position="fill") +
        labs(x="pseudotime", y="frequency", title=mv$sid[1]) +
        scale_fill_manual(values=.pal_roi) +
        geom_col(position="fill") +
        .thm_fig_d("minimal") + theme(
            legend.position="none",
            axis.text=element_blank(),
            panel.grid=element_blank(),
            strip.text.x=element_blank())
})

# saving
pdf(args[[3]], onefile=TRUE, width=8/2.54, height=3/2.54)
for (p in ps) print(p); dev.off()
