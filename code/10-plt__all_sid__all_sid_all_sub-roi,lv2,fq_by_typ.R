# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(ggnewscale)
    library(SingleCellExperiment)
})

# loading
lys <- mapply(x=rep(args[[1]], each=3), y=args[[2]], \(x, y) {
    sce <- readRDS(x); ist <- readRDS(y)
    sub <- gsub(".*(epi|imm|str).*", "\\1", y)
    kid <- (k <- ist$clust)[match(colnames(sce), names(k))]
    data.frame(colData(sce)[c("sid", "typ")], sub, kid)
}, SIMPLIFY=FALSE)

# wrangling
df <- lys |>
    do.call(what=rbind) |>
    filter(!is.na(typ), !is.na(kid)) |>
    mutate(sid=gsub(".*([0-9]{3}).*", "\\1", typ)) |>
    mutate(typ=gsub("^.*_", "\\1", typ)) |>
    mutate(tyq=gsub("[0-9]$", "", typ))

# aesthetics
pal <- c(
    "REF"="seagreen",
    "TVA"="royalblue",
    "TVA1"="dodgerblue",
    "TVA2"="royalblue",
    "TVA3"="slateblue",
    "CRC"="tomato",
    "CRC1"="tomato",
    "CRC2"="indianred")
aes <- list(
    geom_col(
        position="fill", key_glyph="point",
        col="white", alpha=4/5, linewidth=0.1, width=1),
    facet_grid(rows="sub", scales="free_y", space="free"),
    coord_cartesian(expand=FALSE),
    scale_fill_manual(values=pal),
    labs(x="frequency"),
    .thm_fig_d("minimal"), theme(
        legend.position="none",
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank()))

# reorder according to REF proportion
foo <- \(df) {
    ty <- grep("typ|q", names(df), value=TRUE)
    fd <- filter(df, !!sym(ty) == "REF")
    ks <- unique(arrange(fd, p)$kid)
    mutate(df, kid=factor(kid, ks))
}

# joint
fd <- df |>
    group_by(sub, kid, tyq) |> tally() |>
    mutate(p=n/sum(n)) |> ungroup() 
p0 <- ggplot(foo(fd), aes(p, kid, fill=tyq)) + 
    aes + ggtitle(.lab("all", sum(fd$n)))

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
