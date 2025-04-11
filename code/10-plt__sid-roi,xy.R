# loading
sce <- readRDS(args[[1]])

# wrangling
typ <- gsub("^.*_", "", sce$typ)
typ <- factor(typ, names(.pal_roj))
names(typ) <- colnames(sce)

# plotting
gg <- .plt_xy(sce, typ, wcs$x, split=FALSE, na=TRUE) + 
    scale_color_manual(values=.pal_roj, na.value="grey") +
    theme(legend.position="none")

# saving
df <- gg$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
ggsave(args[[2]], gg, unit="cm", width=2+dx/2, height=0.5+dy/2)
