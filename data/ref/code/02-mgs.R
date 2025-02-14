# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(scran)
})

# loading
sce <- readRDS(args[[1]])

# analysis
mgs <- findMarkers(sce, groups=sce[[wcs$kid]], add.summary=TRUE)

# wrangling
df <- lapply(names(mgs), \(kid) {
    gene <- rownames(df <- mgs[[kid]])
    df <- data.frame(gene, kid, df)
    select(df, !starts_with("logFC"))
}) |> bind_rows()
rownames(df) <- NULL

# filtering
df <- df[df$summary.logFC > 0.5, ]
df <- df[df$self.detected > 0.1, ]
table(df$kid)

# saving
write.csv(df, args[[2]])