# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(scater)
    library(scuttle)
    library(circlize)
    library(ComplexHeatmap)
})

# args <- list(
#     list.files("outs", "clu", full.names=TRUE),
#     list.files("outs", "roi", full.names=TRUE),
#     "plts/clu,roi,hm.pdf")
#     source("code/_utils.R")

.pal_roj <- c(
    .pal_roj, 
    MET1="gold",
    MET2="orange",
    MET3="brown2",
    MET4="brown4")

lys <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
    # loading
    sce <- readRDS(x)
    roi <- readRDS(y)
    # aggregation
    pbs <- aggregateAcrossCells(sce, sce$clu,
        use.assay.type=assayNames(sce),
        statistics="mean")
    # quantify by ROI
    cs <- intersect(colnames(sce), colnames(roi))
    if (sce$sid[1] == "241") {
        fq <- prop.table(table(clu=sce[, cs]$clu, roi=roi[, cs]$roi), 1)
        colnames(fq) <- gsub("^.*_ROI([0-9])_(.*)", "\\2\\1", colnames(fq))
    } else {
        fq <- prop.table(table(clu=sce[, cs]$clu, roi=roi[, cs]$typ), 1)
        colnames(fq) <- gsub("^.*_", "", colnames(fq))
    }
    rs <- match(names(.pal_roj), colnames(fq), nomatch=0)
    fq <- as.matrix(unclass(fq[colnames(pbs), rs]))
    fqs <- sweep(fq, 2, colMaxs(fq), `/`)
    chr <- rowData(pbs)$chr
    list(
        id=gsub(".*([0-9]{3}).*", "\\1", x),
        mx=mx <- t(assay(pbs)),
        rd=data.frame(fqs,
            clu=colnames(pbs), 
            tot=.z(rowMeans(abs(mx)))),
        cd=data.frame(row.names=rownames(pbs), chr))
})

# collect results across sections
names(lys) <- sapply(lys, `[[`, "id")
mx <- do.call(rbind, lapply(lys, `[[`, "mx"))
rd <- bind_rows(lapply(lys, `[[`, "rd"), .id="sid")
rd$sid <- factor(rd$sid)
id <- with(rd, paste(sid, clu, sep="."))
rownames(rd) <- rownames(mx) <- id
cd <- lys[[1]]$cd

ps <- lapply(c("all", "240", levels(rd$sid)), \(sid) {
    # subset section(s) of interest
    if (sid == "240") {
        idx <- rownames(rd)[grep("^24", rd$sid)]
        mx <- mx[idx, ]; rd <- rd[idx, ]
    } else if (sid != "all") {
        idx <- rownames(rd)[rd$sid == sid]
        mx <- mx[idx, ]; rd <- rd[idx, ]
        rownames(mx) <- rownames(rd) <- rd$clu
    }
    # selection
    bin <- names(sort(colMeans(abs(mx))))
    bin <- intersect(colnames(mx), tail(bin, 400))
    mx <- mx[, bin]; cd <- cd[bin, , drop=FALSE]
    # palettes
    col <- colorRampPalette(c("navy" ,"ivory", "maroon"))(n <- 101)
    col <- circlize::colorRamp2(seq(-4, 4, l=n), col)
    names(rs) <- rs <- intersect(names(.pal_roj), names(rd))
    pal <- lapply(rs, \(.) colorRamp2(c(0, 1), c("grey95", .pal_roj[.])))
    pal$tot <- colorRamp2(c(-2.5, 0, 2.5), c("blue2", "grey95", "red2"))
    cd$chr <- factor(cd$chr, chr <- paste0("chr", seq_len(22)))
    qal <- list(chr=setNames(pals::stepped(22), chr))
    # aesthetics
    ht_opt$simple_anno_size <- unit(1, "mm")
    lgd <- list(
        title_gp=gpar(fontsize=4),
        labels_gp=gpar(fontsize=3),
        grid_width=unit(1, "mm"),
        grid_height=unit(0.2, "mm"),
        legend_height=unit(1, "cm"))
    # plotting
    Heatmap(t(scale(t(mx))),
        use_raster=TRUE, raster_quality=10,
        col=col, name="z-scaled\nCNV score",
        heatmap_legend_param=lgd,
        row_names_gp=gpar(fontsize=3),
        show_row_names=sid != "all",
        show_column_names=FALSE,
        row_names_side="left",
        row_dend_gp=gpar(lwd=0.2),
        right_annotation=rowAnnotation(
            df=select(rd, tot),
            show_annotation_name=FALSE,
            col=pal, show_legend=FALSE,
            annotation_name_gp=gpar(fontsize=4)),
        left_annotation=rowAnnotation(
            df=select(
                rd[, !colAlls(is.na(rd))], 
                any_of(names(.pal_roj))), 
            col=pal, show_legend=FALSE,
            annotation_name_gp=gpar(fontsize=4)),
        top_annotation=columnAnnotation(
            show_annotation_name=FALSE,
            df=cd, col=qal, show_legend=FALSE,
            annotation_name_gp=gpar(fontsize=4)),
        row_split=if (sid == "all") rd$sid, 
        row_title=if (sid != "all") sid else character(),
        column_split=cd$chr,
        cluster_columns=FALSE,
        column_title_rot=45,
        row_gap=unit(2, "pt"),
        column_gap=unit(2, "pt"),
        row_title_gp=gpar(fontsize=4),
        column_title_gp=gpar(fontsize=4))
})

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    if (. == 1) {
        w <- 15; h <- 9
    } else {
        w <- 12; h <- 5
    }
    pdf(tf[[.]], 
        width=w/2.54, 
        height=h/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
