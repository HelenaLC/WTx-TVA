# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# wrangling
.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
.sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)
ist_by_sid <- split(args[[2]], .sid(args[[2]]))
ist_by_sub <- split(args[[2]], .sub(args[[2]]))
kid_by_sub <- lapply(ist_by_sub, lapply, \(.) readRDS(.)$clust)
kid_by_sub <- lapply(kid_by_sub, \(.) sort(unique(unlist(.))))

ps <- mapply(
    x=args[[1]], y=ist_by_sid,
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- lapply(y, readRDS)
        # wrangling
        names(ist) <- names(sub) <- sub <- .sub(y)
        df <- lapply(sub, \(.) {
            ks <- kid_by_sub[[.]]
            .ks <- gsub("^epi\\.", "", ks)
            kid <- ist[[.]]$clust
            kid <- factor(kid, ks, .ks)
            cd <- colData(sce)[names(kid), ]
            data.frame(cd, sub=., kid) |>
                mutate(roi=gsub("^.*_", "", typ)) |>
                mutate(roi=factor(roi, names(.pal_roj))) |>
                mutate(roi=droplevels(roi))
        })
        # plotting
        ps <- lapply(df, \(fd) {
            gg <- .plt_fq(fd, 
                x="roi", y="kid", hc=FALSE) +
                coord_equal(15, expand=FALSE) + 
                .thm_fig_d("minimal", "f") + theme(
                    legend.margin=margin(),
                    legend.title=element_blank(),
                    axis.title=element_blank(),
                    axis.ticks=element_blank(),
                    axis.text.y=element_blank(),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
            gg$guides$guides$fill$params$override.aes$size <- 1
            sid <- fd$sid[1]; sub <- fd$sub[1]
            id <- ifelse(sub == "epi", paste0(sid, "-", sub), sub)
            gg + ggtitle(.lab(as.character(id), nrow(fd)))
        })
        wrap_plots(ps, nrow=1)
    })

# saving
pdf(args[[3]], onefile=TRUE, width=10/2.54, height=4/2.54)
for (p in ps) print(p); dev.off()
