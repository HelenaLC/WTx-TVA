# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

.sid_one <- \(.) gsub(".*-(.*)\\.rds", "\\1", .)
.sid_two <- \(.) gsub(".*-(.*),.*\\.rds", "\\1", .)

# wrangling
names(sce) <- .sid_one(basename(args[[1]]))
df <- lapply(names(sce), \(sid) {
    i <- .sid_two(basename(args[[2]])) == sid
    k <- unlist(lapply(ist[i], \(lys) lys$clust))
    k <- k[match(colnames(sce[[sid]]), names(k))]
    xy <- grep("global_mm", names(colData(sce[[sid]])))
    xy <- setNames(colData(sce[[sid]])[xy], c("x", "y"))
    data.frame(xy, k, sid)
}) |> do.call(what=rbind)

# analysis
i <- split(seq_len(nrow(df)), df$sid)
fq <- lapply(i, \(.) {
    # cellular neighborhoods
    yx <- as.matrix(df[., c("x", "y")])
    nn <- nn2(yx, r=0.05, k=(k <- 300)+1, searchtype="radius")
    range(rowSums((is <- nn$nn.idx[, -1]) > 0))
    is[is == 0] <- NA
    # subpopulation frequencies
    ks <- factor(df$k[.][c(is)], sort(unique(df$k[.])))
    ns <- array2DF(by(ks, rep(seq(nrow(yx)), k), table))
    fq <- prop.table(as.matrix(ns[, -1]), 1)
    df <- data.frame(fq, i=df$k[.], check.names=FALSE)
    fd <- pivot_longer(df, -i, names_to="j", values_to="p") 
}) |> bind_rows(.id="sid")
head(fq)

mu <- fq |>
    group_by(i, j, sid) |>
    # average across cells
    summarize_at("p", mean, na.rm=TRUE)
um <- mu |>
    # average across section
    summarize_at("p", mean, na.rm=TRUE) 

# aesthetics
thm <- list(
    geom_tile(),
    coord_equal(expand=FALSE),
    scale_fill_gradient2(
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2),
        low="royalblue3", mid="ivory", high="tomato3"),
    theme_minimal(6), theme(
        panel.grid=element_blank(),
        axis.title=element_blank(),
        plot.title=element_text(hjust=0.5),
        legend.key.size=unit(0.4, "lines"),
        panel.background=element_rect(fill=NA),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)))

# plotting
gg <- um |>
    group_by(j) |>
    mutate(p=.z(p)) |>
    filter(!is.na(i), !is.na(j))
y <- pivot_wider(gg, names_from="j", values_from="p")
z <- as.matrix(y[, -1]); rownames(z) <- y[[1]]
p0 <- ggplot(gg, aes(j, i,fill=p)) + thm +
    scale_x_discrete(limits=.xo(z)) +
    scale_y_discrete(limits=.xo(z)) 

gg <- mu |>
    group_by(sid, j) |>
    mutate(p=.z(p)) |>
    filter(!is.na(i), !is.na(j))
ps <- lapply(split(seq(nrow(gg)), gg$sid), \(.) {
    y <- pivot_wider(gg[., ], names_from="j", values_from="p")
    z <- as.matrix(y[, -1]); rownames(z) <- y[[1]]
    nk <- nrow(z); id <- gg$sid[.][1]
    ggplot(gg[., ], aes(j, i,fill=p)) + thm +
        scale_x_discrete(limits=.xo(z)) +
        scale_y_discrete(limits=.xo(z)) +
        ggtitle(bquote(bold(.(id))~"(N ="~.(nk)*")")) +
        theme(legend.position="none", axis.text=element_text(size=3))
}) |> wrap_plots(nrow=2) 

pdf(args[[3]], width=18/2.54, height=10/2.54, onefile=TRUE)
for (p in list(p0, ps)) print(p); dev.off()
