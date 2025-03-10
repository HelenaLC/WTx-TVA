# loading
sce <- readRDS(args[[1]])

# wrangling
typ <- gsub(".*(REF|TVA|CRC).*", "\\1", sce$typ)
typ <- factor(typ, c("REF", "TVA", "CRC"))
names(typ) <- colnames(sce)

# plotting
gg <- .plt_xy(sce, typ, wcs$sid, split=FALSE, na=TRUE) + 
    scale_color_manual(values=.pal_roi, na.value="grey") +
    theme(legend.position="none")

# saving
df <- gg$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=list(gg), nrow=1, ncol=1, top=FALSE)
ggsave(args[[2]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
