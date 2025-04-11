# dependencies
suppressPackageStartupMessages({
    library(scater)
    library(scuttle)
    library(circlize)
    library(ComplexHeatmap)
})

ps <- mapply(x=args[[1]], y=args[[2]], SIMPLIFY=FALSE, \(x, y) {
    # loading
    sce <- readRDS(x)
    roi <- readRDS(y)
    # aggregation
    pbs <- aggregateAcrossCells(sce, sce$clu, 
        use.assay.type=assayNames(sce), 
        statistics="mean")
    # quantify by ROI
    cs <- intersect(colnames(sce), colnames(roi))
    fq <- prop.table(table(clu=sce[, cs]$clu, roi=roi[, cs]$typ), 1)
    colnames(fq) <- gsub("^.*_", "", colnames(fq))
    fq <- fq[, match(names(.pal_roj), colnames(fq), nomatch=0)]
    fq <- unclass(fq[colnames(pbs), ])
    colData(pbs) <- cbind(colData(pbs), fq)
    # aesthetics
    sid <- gsub(".*([0-9]{3}).*", "\\1", x)
    bin <- tail(names(sort(rowMeans(abs(assay(pbs))))), 200)
    bin <- intersect(rownames(sce), bin)
    col <- colorRampPalette(c("navy" ,"ivory", "maroon"))(n <- 101)
    col <- circlize::colorRamp2(seq(-4, 4, l=n), col)
    names(roi) <- roi <- colnames(fq)
    pal <- lapply(roi, \(.) colorRamp2(c(0, 1), c("grey95", .pal_roj[.])))
    chr <- rowData(pbs[bin, ])$chr; names(chr) <- bin
    chr <- factor(chr, unique(chr[order(as.integer(gsub("chr", "", chr)))]))
    qal <- list(chr=setNames(pals::stepped(nlevels(chr)), levels(chr)))
    # plotting
    mtx <- t(scale(assay(pbs[bin, ])))
    ht_opt$simple_anno_size <- unit(1, "mm")
    lgd <- list(
        title_gp=gpar(fontsize=4),
        labels_gp=gpar(fontsize=3),
        grid_width=unit(1, "mm"),
        grid_height=unit(0.2, "mm"),
        legend_height=unit(1, "cm"))
    Heatmap(mtx, col=col, 
        name="z-scaled\nCNV score",
        heatmap_legend_param=lgd,
        row_names_side="left",
        row_names_gp=gpar(fontsize=4),
        show_row_names=TRUE,
        show_column_names=FALSE,
        cluster_columns=FALSE,
        left_annotation=rowAnnotation(
            df=data.frame(fq), col=pal, 
            annotation_name_gp=gpar(fontsize=4),
            annotation_legend_param=c(lgd, list(at=c(0, 1)))),
        top_annotation=columnAnnotation(
            df=data.frame(chr), col=qal,
            annotation_name_gp=gpar(fontsize=4),
            annotation_legend_param=lgd))
})

# saving
pdf(args[[3]], onefile=TRUE, width=12/2.54, height=5/2.54)
for (p in ps) print(p); dev.off()
