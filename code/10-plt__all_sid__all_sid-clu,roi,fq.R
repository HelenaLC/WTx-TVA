# loading
library(SingleCellExperiment)
ps <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        roi <- readRDS(y)
        # wrangling
        roi <- roi[, !is.na(roi$typ)]
        idx <- intersect(colnames(sce), colnames(roi))
        sce <- sce[, idx]; roi <- roi[, idx]
        sid <- gsub(".*([0-9]{3}).*", "\\1", x)
        roi <- gsub("^.*_", "", roi$typ)
        roi <- factor(roi, names(.pal_roj))
        df <- data.frame(
            clu=sce[, idx]$clu, 
            roi=droplevels(roi))
        # plotting
        n <- length(unique(df$clu))
        m <- length(unique(df$roi))
        p <- .plt_fq(df, "clu", "roi", sid, h=TRUE) + 
            scale_fill_manual("roi", values=.pal_roj) +
            labs(y=NULL) 
        q <- .plt_fq(df, "roi", "clu", hc=FALSE, h=TRUE) + 
            scale_fill_manual("cnv", values=pals::kelly()) +
            theme(plot.title=element_blank())
        aes <- list(shape=21, stroke=0, size=1)
        wrap_plots(p, q, ncol=1, heights=c(n, m), guides="collect") & 
            guides(fill=guide_legend(override.aes=aes)) &
            theme(
                legend.margin=margin(),
                legend.spacing=unit(0.2, "lines"),
                legend.justification=c(0, 0.5),
                axis.text.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.y=element_text(hjust=1))
    })

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    h <- sum(sapply(ps[[.]], \(.) length(unique(.$data[[2]]))))
    pdf(tf[[.]], width=6/2.54, height=(0.5+h/5)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
