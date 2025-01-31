# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(AUCell)
    library(scuttle)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# wrangling
.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
names(sce) <- sid <- .sid(args[[1]])

df <- lapply(sid, \(sid) {
    .sce <- sce[[which(.sid(args[[1]]) == sid)]]
    .ist <- ist[which(.sid(args[[2]]) == sid)]
    kid <- unlist(lapply(.ist, \(lys) lys$clust))
    kid <- kid[match(colnames(.sce), names(kid))]
    mtx <- cor(t(assay(.sce)))
    data.frame(sid, p=c(mtx),
        i=rownames(mtx)[c(row(mtx))],
        j=rownames(mtx)[c(col(mtx))])
}) |> do.call(what=rbind)

# plotting
.p <- \(df) {
    y <- pivot_wider(df, names_from="j", values_from="p")
    z <- as.matrix(y[, -1]); rownames(z) <- y[[1]]
    ggplot(df, aes(i, j, fill=p)) +
        scale_x_discrete(limits=.yo(z)) +
        scale_y_discrete(limits=rev(.yo(z))) +
        scale_fill_gradient2("corr.", 
            low="royalblue3", mid="ivory", high="tomato3",
            na.value="lightgrey", limits=c(-1, 1), n.breaks=3) + 
        coord_equal(expand=FALSE) +
        geom_tile() +
        theme_bw(6) + theme(
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            panel.grid=element_blank(),
            legend.key.size=unit(0.4, "lines"),
            legend.title=element_text(vjust=1),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
}
    
fd <- df |>
    group_by(i, j) |>
    summarise_at("p", mean) |>
    mutate(p=case_when(i == j ~ NA, TRUE ~ p)) 

p0 <- .p(fd)

ps <- lapply(sid, \(.) {
    fd <- df |>
        filter(sid == .) |> 
        select(-sid) |>
        mutate(p=case_when(i == j ~ NA, TRUE ~ p)) 
    nc <- format(ncol(sce[[.]]), big.mark=",")
    .p(fd) + ggtitle(bquote(bold(.(.))~"(N ="~.(nc)*")"))
})

# saving
pdf(args[[3]], width=8/2.54, height=8/2.54, onefile=TRUE)
for (p in c(list(p0), ps)) print(p); dev.off()
