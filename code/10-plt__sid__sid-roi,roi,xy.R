# loading
sce <- readRDS(args[[1]])

# wrangling
pal <- c("seagreen", "royalblue", "tomato")
typ <- gsub(".*(REF|TVA|CRC).*", "\\1", sce$typ)
typ <- factor(typ, c("REF", "TVA", "CRC"))
names(typ) <- colnames(sce)

# plotting
gg <- .plt_xy(sce, typ, wcs$sid)[[1]] + 
    scale_color_manual(values=pal) +
    theme(legend.position="none")

# saving
df <- gg$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
ggsave(args[[3]], gg, units="cm", width=2+dx/2, height=0.5+dy/2)
