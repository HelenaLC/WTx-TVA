# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(slingshot)
    library(SingleCellExperiment)
})

.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
.sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)
.split <- \(., by) split(., get(paste0(".", by))(.))

# for each compartment, get unique 
# subpopulations across sections
kid <- lapply(.split(args[[3]], "sub"), 
    \(.) unique(unlist(lapply(., 
        \(.) readRDS(.)$clust))))

ps <- mapply(
    x=args[[1]], y=args[[2]], 
    z=.split(args[[3]], "sid"), \(x, y, z) {
        # loading
        res <- readRDS(x)
        sce <- readRDS(y)
        ist <- lapply(z, readRDS)
        # wrangling
        xy <- "Center(X|Y)_global_mm"
        xy <- grep(xy, names(colData(sce)))
        xy <- as.matrix(colData(sce)[xy])
        t <- .q(slingAvgPseudotime(res$slingshot))
        # analysis
        names(sub) <- names(ist) <- sub <- .sub(z)
        idx <- lapply(ist, \(lys) intersect(colnames(sce), names(lys$clust)))
        df <- lapply(sub, \(.) {
            k <- 300; i <- TRUE; if (. == "epi") { k <- k+1; i <- -1 }
            nn <- nn2(xy[idx[[.]], ], xy[idx$epi, ], searchtype="radius", r=0.05, k=k)
            is <- nn$nn.idx[, i]; is[is == 0] <- NA
            ks <- ist[[.]]$clust; ks <- factor(ks, sort(kid[[.]]))
            ns <- lapply(seq(nrow(is)), \(.) table(ks[is[., ]]))
            fq <- prop.table(do.call(rbind, ns), 1)
            df <- data.frame(check.names=FALSE, t=t[idx$epi], fq)
            fd <- pivot_longer(df, -t, names_to="k", values_to="p")
        })
        # plotting
        xs <- seq(-(dx <- 0.005), 1.005, 0.01)
        ps <- lapply(sub, \(.) {
            fd <- df[[.]] |>
                mutate(
                    t=cut(t, breaks=xs),
                    t=xs[as.integer(t)]+dx,
                    t=factor(t, paste((xs+dx)[-length(xs)]))) |>
                group_by(k, t) |>
                summarize_at("p", mean, na.rm=TRUE)
            if (. == "epi") fd$k <- gsub("^epi\\.", "", fd$k)
            ggplot(fd, aes(t, p, fill=k)) +
                (if (. == "epi") ggtitle(.lab(paste(sce$sid[1]), ncol(res)))) +
                geom_col(width=1, key_glyph="point", position="fill") +
                labs(x=if (. == "str") "pseudotime", y="frequency") +
                scale_fill_manual(., values=.pal_kid) +
                coord_cartesian(expand=FALSE) 
        }) 
        wrap_plots(ps, ncol=1, guides="collect") &
            .thm_fig_d("minimal", "f") & theme(
                aspect.ratio=1/2,
                plot.margin=margin(),
                legend.margin=margin(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                legend.spacing=unit(0.2, "lines"),
                legend.key.size=unit(0, "lines")) &
            guides(fill=guide_legend(override.aes=list(
                shape=21, stroke=0, size=1)))
    })

# saving
pdf(args[[4]], onefile=TRUE, width=9/2.54, height=9/2.54)
for (p in ps) print(p); dev.off()
