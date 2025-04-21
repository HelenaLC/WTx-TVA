# dependencies
suppressPackageStartupMessages({
    library(zoo)
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(SingleCellExperiment)
})

# loading
raw <- readRDS(args[[1]])
fil <- readRDS(args[[2]])

# wrangling
ex <- !colnames(raw) %in% colnames(fil)
ex <- setNames(ex, colnames(raw))
df <- data.frame(colData(raw), ex)
df <- df[order(df$ex), ]

# plotting
p <- round(100*mean(!ex), 2)
n <- format(ncol(raw), big.mark=",")
m <- format(ncol(fil), big.mark=",")
gg <- .plt_xy(df, df$ex, split=FALSE) +
    theme(legend.position="none") +
    scale_color_manual(values=c("lavender", "purple")) +
    ggtitle(bquote(bold(.(wcs$sid))~"(N ="~.(n)*";"~.(m)~"-"~.(p)*"%)"))

# saving
df <- gg$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
ggsave(args[[3]], gg, unit="cm", width=2+dx/2, height=0.5+dy/2)
