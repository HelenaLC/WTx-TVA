# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])

# wrangling
cd <- colData(sce)
ts <- grep("seudo", names(cd))
pt <- unname(.q(cd[[ts]]))
ex <- grep("sling", names(cd))
df <- data.frame(cd[-ex])

# plotting
t <- unname(.q(pt))
p <- .plt_xy(df, t, wcs$sid)

# saving
df <- p$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=list(p), nrow=1, ncol=1, top=FALSE)
ggsave(args[[2]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
