# loading
sce <- readRDS(args[[1]])

# wrangling
typ <- gsub(".*(REF|TVA|CRC).*", "\\1", sce$typ)
typ <- factor(typ, c("REF", "TVA", "CRC"))
typ <- setNames(typ, colnames(sce))

# plotting
gg <- .plt_xy(sce, typ, wcs$sid, split=FALSE) + 
    scale_color_manual(values=.pal_roi) +
    theme(legend.position="none")

# saving
df <- gg$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
ggsave(args[[3]], gg, units="cm", width=2+dx/2, height=0.5+dy/2)
