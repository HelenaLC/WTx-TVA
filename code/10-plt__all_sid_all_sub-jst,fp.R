# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
    library(InSituType)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
pat <- ".*-(.*),(.*)\\.rds"
sid <- gsub(pat, "\\1", basename(args[[1]]))
sub <- gsub(pat, "\\2", basename(args[[1]]))
idx <- split(seq_along(ist), sub)

# plotting by subset
ps <- lapply(idx, \(.) {
    ll <- lapply(ist[.], \(jst) data.frame(
        jst$logliks, check.names=FALSE))
    ps <- lapply(ist[.], \(jst) data.frame(
        jst$profiles, check.names=FALSE))
    ll <- as.matrix(bind_rows(ll))
    ps <- as.matrix(bind_rows(ps))
    ll[is.na(ll)] <- -Inf
    ps[is.na(ps)] <- 0
    y <- flightpath_layout(logliks=ll, profiles=ps)
    k <- unlist(lapply(ist[.], \(jst) jst$clust))
    df <- data.frame(y$cellpos, k)
    nc <- format(nrow(df), big.mark=",")
    cs <- lapply(split(seq(nrow(df)), df$k), 
        \(.) sample(., min(1e4, length(.))))
    df <- df[sample(unlist(cs)), ]
    fd <- summarize(
        group_by(df, k), 
        across(c("x", "y"), median))
    ggplot(df, aes(x, y, col=k)) + 
        ggtitle(bquote(bold(.(sub[.][1]))~"(N ="~.(nc)*")")) +
        geom_point_rast(shape=16, stroke=0, size=0.1, alpha=0.2) +
        geom_text(data=fd, aes(label=k), size=1.5, col="black")
})

# plotting by section
qs <- lapply(idx, \(i) {
    df <- lapply(seq_along(sid[i]), \(j) {
        jst <- ist[i][[j]]
        y <- flightpath_layout(
            logliks=jst$logliks, 
            profiles=jst$profiles)
        df <- data.frame(y$cellpos, k=jst$clust)
        cs <- lapply(split(seq(nrow(df)), df$k), 
            \(.) sample(., min(1e4, length(.))))
        data.frame(
            df[sample(unlist(cs)), ], 
            sid=sid[i][j], sub=sub[i][j])
    }) |> do.call(what=rbind)
    fd <- summarize(
        group_by(df, sid, sub, k), 
        across(c("x", "y"), median))
    ggplot(df, aes(x, y, col=k)) + 
        facet_wrap(~sid, nrow=2, scales="free") +
        geom_point_rast(shape=16, stroke=0, size=0.1, alpha=0.2) +
        geom_text(data=fd, aes(label=k), size=1.5, col="black")
})

gg <- lapply(c(ps, qs), \(p)
    p + theme_void(6) + theme(
    panel.grid=element_blank(),
    plot.title=element_text(hjust=0.5),
    strip.text=element_text(face="bold", margin=margin(b=2)),
    aspect.ratio=1, legend.position="none",
    panel.background=element_rect(color="grey", fill=NA)) +
    scale_color_manual(NULL, values=.pal))
gg <- c(list(wrap_plots(gg[seq(3)], nrow=1)), gg[-seq(3)])

# saving
pdf(args[[2]], width=12/2.54, height=6.5/2.54, onefile=TRUE)
for (g in gg) print(g); dev.off()
