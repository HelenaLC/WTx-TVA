# loading
require(dplyr, quietly=TRUE)
ist <- lapply(args[[1]], readRDS)

# wrangling
pat <- ".*-(.*),(.*)\\.rds"
sid <- gsub(pat, "\\1", args[[1]])
sub <- gsub(pat, "\\2", args[[1]])
idx <- split(seq_along(ist), sub)

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
    lapply(seq_along(sid[i]), \(j) {
        .plt_fp(ist[i][[j]], sid[i][j])
    }) |>
        wrap_plots(nrow=2) & 
        theme(plot.margin=margin(r=2))
})

# saving
pdf(args[[2]], width=12/2.54, height=6.5/2.54, onefile=TRUE)
for (p in c(list(ps), qs)) print(p); dev.off()
