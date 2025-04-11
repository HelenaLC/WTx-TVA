ps <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        ist <- readRDS(x)
        sce <- readRDS(y)
        # plotting
        sid <- gsub(".*([0-9]{3}).*", "\\1", x)
        sub <- gsub(".*(epi|imm|str).*", "\\1", x)
        .plt_xy(sce, ist$clust, paste0(sid, "-", sub))
    }) |> Reduce(f=c)

# saving
df <- ps[[1]]$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=4+dx/2, height=0.5+dy/2)
