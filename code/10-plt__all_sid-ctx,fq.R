#  loading
dfs <- lapply(args[[1]], readRDS)
df <- do.call(rbind, dfs)

# plotting
p1 <- .plt_fq(df, "sid", "ctx", "", hc=FALSE)
p2 <- .plt_fq(df, "ctx", "sid", "", hc=TRUE)
ws <- c(
    length(unique(df$sid)),
    length(unique(df$ctx)))
gg <- wrap_plots(p1, p2, widths=ws) & 
    theme(plot.title=element_blank(), 
        axis.title.x=element_blank())

# saving
ggsave(args[[2]], gg, units="cm", width=8, height=6)
