# loading
require(dplyr, quietly=TRUE)
roi <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# wrangling
pat <- ".*-(.*),(.*)\\.rds"
sid <- gsub(pat, "\\1", args[[2]])
sub <- gsub(pat, "\\2", args[[2]])
idx <- split(seq_along(ist), sub)

args <- lapply(args, \(.) grep("241", ., invert=TRUE, value=TRUE))

# plotting by subset
ps <- lapply(idx, \(i) {
    ks <- unlist(lapply(ist[i], \(jst) jst$clust))
    ll <- lapply(ist[i], \(jst) data.frame(jst$logliks))
    ps <- lapply(ist[i], \(jst) data.frame(jst$profiles))
    ll <- as.matrix(bind_rows(ll)); ll[is.na(ll)] <- -Inf
    ps <- as.matrix(bind_rows(ps)); ps[is.na(ps)] <- 0
    x <- list(clust=ks, logliks=ll, profiles=ps)
    .plt_fp(x, sub[i][1])
}) |>
    wrap_plots(nrow=1) & 
    theme(plot.margin=margin(r=2))

# plotting by section
qs <- lapply(idx, \(i) { 
    ps <- lapply(seq_along(sid[i]), \(j)
        .plt_fp(ist[i][[j]], sid[i][j]))
    qs <- lapply(seq_along(ps), \(i) {
        p <- ps[[i]]
        roj <- roi[[i]]
        cs <- match(rownames(p$data), colnames(roj))
        typ <- gsub("^.*_", "", roj$typ)[cs]
        p$data$typ <- typ
        p$mapping$colour <- quo(typ)
        p$data <- p$data[!is.na(p$data$typ), ]
        p + scale_color_manual(values=.pal_roj)
    })
    lapply(list(ps, qs), \(p) 
        p |> wrap_plots(nrow=2) & 
        theme(plot.margin=margin(r=2)))
}) |> Reduce(f=c)

# saving
pdf(args[[3]], width=12/2.54, height=6.5/2.54, onefile=TRUE)
for (p in c(list(ps), qs)) print(p); dev.off()
