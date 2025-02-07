# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(ggnewscale)
    library(RColorBrewer)
    library(SingleCellExperiment)
})

# loading
df <- mapply(x=rep(args[[1]], each=3), y=args[[2]], \(x, y) {
    sce <- readRDS(x); ist <- readRDS(y)
    sub <- gsub(".*(epi|imm|str).*", "\\1", y)
    kid <- (k <- ist$clust)[match(colnames(sce), names(k))]
    data.frame(colData(sce)[c("sid", "typ")], sub, kid)
}, SIMPLIFY=FALSE) |> 
    do.call(what=rbind) |>
    mutate(roi=typ) |>
    mutate(sid=gsub("^([0-9]{3}).*", "\\1", roi)) |>
    mutate(typ=gsub(".*(REF|TVA|CRC).*", "\\1", typ)) |>
    mutate(typ=factor(typ, c("REF", "TVA", "CRC"))) |>
    filter(!is.na(typ), !is.na(kid))

# aesthetics
aes <- list(
    scale_fill_manual(values=.pal_roi),
    facet_grid(rows="sub", scales="free_y", space="free"),
    geom_col(alpha=2/3, position="fill", width=1, key_glyph="point"),
    scale_x_continuous(
        "proportion of cells", n.breaks=6, 
        labels=scales::label_percent()),
    coord_cartesian(expand=FALSE),
    .thm_fig_d("minimal"), theme(
        legend.position="none",
        panel.grid=element_blank(),
        axis.title.y=element_blank()))

# reorder according to REF proportion
foo <- \(df) {
    o <- df |>
        filter(typ == "REF") |> 
        arrange(p) |> pull(kid) |> unique()
    mutate(df, kid=factor(kid, o))
}

# joint
fd <- df |>
    group_by(sub, kid, typ) |> tally() |>
    mutate(p=n/sum(n)) |> ungroup() 
p0 <- ggplot(foo(fd), aes(p, kid, fill=typ)) + 
    aes + ggtitle(.lab("", sum(fd$n)))

# split
fd <- df |>
    group_by(sub, kid, sid, typ) |> tally() |>
    mutate(p=n/sum(n)) |> ungroup() 
ps <- by(fd, fd$sid, \(.) {
    ggplot(foo(.), aes(p, kid, fill=typ)) + 
    aes + ggtitle(.lab(.$sid[1], sum(.$n)))
})

# saving
gs <- gridExtra::marrangeGrob(
    grobs=c(list(p0), ps), 
    nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=6, height=8)
