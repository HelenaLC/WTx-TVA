# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(slingshot)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])

# wrangling
ex <- grep("sling", names(cd <- colData(sce)))
pt <- slingAvgPseudotime(sce$slingshot)
df <- data.frame(cd[-ex], pt)

# plotting
p <- .plt_xy(df, unname(.q(pt)), wcs$x) +
    scale_color_gradientn(
        "q-scaled\npseudotime",
        colors=pals::jet(), limits=c(0, 1),
        n.breaks=2, labels=c("early", "late"))

# saving
df <- p$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=list(p), nrow=1, ncol=1, top=FALSE)
ggsave(args[[2]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
