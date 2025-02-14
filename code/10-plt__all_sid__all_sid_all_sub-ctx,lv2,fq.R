# loading
dfs <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# wrangling
.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
.sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)
fd <- mapply(lys=ist, 
    sid=.sid(args[[2]]), 
    sub=.sub(args[[2]]), 
    \(lys, sid, sub) {
        data.frame(sid, 
            kid=lys$clust, 
            cid=names(lys$clust))
    }, SIMPLIFY=FALSE) |>
    do.call(what=rbind)
df <- do.call(rbind, dfs); df$kid <- NULL
df <- dplyr::left_join(df, fd, by=c("sid", "cid"))

# plotting
ps <- by(df, df$sub, \(.) .plt_fq(., "ctx", "kid", .$sub[1], hc=TRUE, h=TRUE))
xo <- ps[[1]]$scales$scales[[2]]$limits
ps <- lapply(ps, `+`, scale_x_discrete(limits=xo))
gg <- wrap_plots(ps, nrow=1) & theme(axis.title.y=element_blank())

# saving
ggsave(args[[3]], gg, units="cm", width=15, height=5)
