# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tidytext)
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

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

# grouping
fd <- df |>
    filter(!is.na(pt)) |>
    arrange(pt) |>
    mutate(grp=case_when(
        grepl("RT[0-9]?$", roi) ~ "RT",
        grepl("RC[0-9]?$", roi) ~ "RC",
        grepl("T[0-9]?C$", roi) ~ "TC",
        TRUE ~ "TT"), grp=factor(grp, c("RT", "TT", "RC", "TC")))

# reorder
yo <- fd |>
    group_by(roi) |>
    summarise_at("pt", mean) |>
    arrange(-pt) |> pull(roi)

# plotting
gg <- ggplot(fd, aes(pt, roi, fill=pt)) +
    facet_grid(grp~1, space="free", scales="free_y") +
    labs(x="pseudotime", y=NULL) +
    scale_y_discrete(limits=\(.) intersect(yo, .)) +
    geom_col(position="fill") +
    scale_fill_gradientn(
        "q-scaled\npseudotime", n.breaks=6,
        limits=c(0, 1), colors=pals::jet()) +
    .thm_fig_c("minimal") + theme(
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x=element_blank())

# saving
ggsave(args[[3]], gg, units="cm", width=8, height=8)
