# pal ----
.pal_sub <- c(imm="cyan", epi="gold", str="magenta")
.pal_kid <- unname(pals::trubetskoy())
.pal <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# thresholded z-normalization
.z <- \(x, th=2.5) {
    if (is.null(dim(x))) {
        x[x < 0] <- 0
        sd <- sd(x, na.rm=TRUE)
        if (is.na(sd)) sd <- 1
        x <- x-mean(x, na.rm=TRUE)
        if (sd != 0) x <- x/sd
    } else {
        mus <- colMeans(x)
        sds <- colSds(x)
        x <- sweep(x, 2, mus, `-`)
        x <- sweep(x, 2, sds, `/`)
    }
    x[x > +th] <- +th
    x[x < -th] <- -th
    return(x)
}

# upper/lower quantile scaling
.q <- \(x, margin=1, q=0.01) {
    if (length(q) == 1)
        q <- c(q, 1-q)
    if (!is.matrix(x)) {
        qs <- quantile(x, q, na.rm=TRUE)
        x <- (x-qs[1])/diff(qs)
    } else {
        qs <- c(rowQuantiles, colQuantiles)[[margin]]
        qs <- matrix(qs(x, probs=q), ncol=2)
        x <- switch(margin, 
            `1`=(x-qs[, 1])/(qs[, 2]-qs[, 1]), 
            `2`=t((t(x)-qs[, 1])/(qs[, 2]-qs[, 1])))
    }
    x[x < 0] <- 0
    x[x > 1] <- 1
    return(x)
}

# prettified plot title in the style of
# 'title (N = count)' with bold 'title'
.lab <- \(x, n=NULL) {
    if (is.null(n)) {
        if (is.null(x)) "" else # blanc
            bquote(bold(.(x)))  # 'x' only
    } else {
        n <- format(n, big.mark=",")
        if (is.null(x)) bquote("N ="~.(n)) else # 'n' only
            bquote(bold(.(x))~"(N ="~.(n)*")")  # both
    }
}

# run PCA using 'RSpectra'
.pca <- \(sce, gs=TRUE, k=30) {
    library(RSpectra)
    library(SingleCellExperiment)
    mtx <- as(t(logcounts(sce[gs, ])), "dgCMatrix")
    mtx <- svds(mtx, k=k, opts=list(center=TRUE, scale=FALSE))
    mtx$d <- mtx$d * sqrt(nrow(mtx$u)-1); pca <- mtx$u %*% diag(mtx$d)
    dimnames(pca) <- list(colnames(sce), paste0("PC", seq(ncol(pca))))
    reducedDim(sce, "PCA") <- pca
    return(sce)
}

# principal component regression
.pcr <- \(sce, x) {
    y <- reducedDim(sce, "PCA")
    lapply(x, \(.) {
        z <- summary(lm(y ~ sce[[.]]))
        r2 <- sapply(z, \(.) .$adj.r.squared)
        data.frame(x=., pc=seq_along(r2), r2)
    }) |> do.call(what=rbind)
}

# run 'InSituType'
# gs = features to use
# nk = number of clusters
.ist <- \(sce, nk, gs=TRUE, pbs=NULL, bkg=TRUE, ns=c(1e4, 2e4, 1e5)) {
    # dependencies
    library(InSituType)
    library(SingleCellExperiment)
    # load counts
    mtx <- counts(sce[gs, ])
    mtx <- as(t(mtx), "dgCMatrix")
    # cohorting based on IF data
    j <- names(cd <- colData(sce))
    i <- grep("^Mean", j, value=TRUE)
    i <- setdiff(i, "Mean.G")
    i <- c("Area", "AspectRatio", i)
    coh <- fastCohorting(as.matrix(cd[i]))
    # background estimation
    neg <- grep("^neg", altExpNames(sce), ignore.case=TRUE, value=TRUE)
    neg <- sce$nCount_negprobes/nrow(altExp(sce, neg))
    # update reference profiles
    pbs <- if (!is.null(pbs)) {
        bkg <- if (bkg) {
            rna <- sce$nCount_RNA
            rna*mean(neg)/mean(rna) 
        }
        updateReferenceProfiles(
            reference_profiles=pbs, reference_sds=NULL,
            counts=mtx, neg=neg, bg=bkg)$updated_profiles
    }
    # clustering
    insitutype(mtx, 
        reference_profiles=pbs,
        update_reference_profiles=FALSE,
        neg=neg, cohort=coh, n_clusts=nk,
        n_chooseclusternumber=ns[1],
        n_benchmark_cells=ns[1],
        n_phase1=ns[1],
        n_phase2=ns[2],
        n_phase3=ns[3])
}

# subset 'InSituType' clustering results
.sub_ist <- \(x, y) {
    if (all(y %in% x$clust)) {
        i <- x$clust %in% y
        x$clust <- x$clust[i]
        j <- y
    } else {
        m <- match(y, names(x$clust), nomatch=0)
        x$clust <- x$clust[m]
        i <- names(x$clust)
        j <- unique(x$clust)
    }
    x$prob <- x$prob[i]
    x$profiles <- x$profiles[, j]
    x$logliks <- x$logliks[i, j]
    return(x)
}

# relabel 'InSituType' clustering results
.lab_ist <- \(x, df) {
    if (!is.data.frame(df))
        df <- data.frame(names(df), df, check.names=FALSE)
    i <- x$clust
    j <- df[[2]][match(i, df[[1]])]
    j[.] <- i[. <- is.na(j)]
    names(j) <- names(i)
    x$clust <- j
    
    i <- colnames(x$logliks)
    j <- df[[2]][match(i, df[[1]])]
    j[.] <- i[. <- is.na(j)]
    colnames(x$logliks) <- j
    colnames(x$profiles) <- j
    
    i <- split(seq_along(j), j)
    x$logliks <- sapply(i, \(.) rowMeans(x$logliks[, ., drop=FALSE]))
    x$profiles <- sapply(i, \(.) rowMeans(x$profiles[, ., drop=FALSE]))
    
    return(x)
}

# compute each cell's distance (in px) 
# to t(op), b(ottom), r(ight), l(eft) 
# field of view (FOV) borders 
.d2b <- \(obj, xy=c("CenterX_local_px", "CenterY_local_px")) {
    df <- if (is(obj, "SingleCellExperiment"))
        data.frame(colData(obj)) else obj
    is <- split(seq(nrow(df)), df$fov)
    xy <- match(xy, names(df))
    yx <- c(".x", ".y")
    names(df)[xy] <- yx
    ds <- lapply(is, \(.) 
        with(df[., ], cbind(
            .x-min(.x), max(.x)-.x, 
            .y-min(.y), max(.y)-.y)))
    ds <- as.matrix(do.call(rbind, ds))
    ds <- ds[order(unlist(is)), ]
    bs <- c("l", "r", "b", "t")
    colnames(ds) <- bs; df[, bs] <- ds
    df <- df[, setdiff(names(df), yx)]
    if (!is(obj, "SingleCellExperiment")) return(df) 
    colData(obj) <- cbind(colData(obj), ds); obj
}

.plt_qc_ns <- \(x, i, id="") {
    # dependencies
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
    # wrangling
    df <- data.frame(colData(x))
    df <- list(raw=df, fil=df[i, ])
    df <- bind_rows(df, .id="did")
    ys <- c("RNA", "negprobes", "falsecode")
    ps <- lapply(ys, \(y) {
        switch(y, 
            RNA={bw=NULL;nb=30;w=10}, 
            {bw=1;nb=NULL;w=0.5})
        z <- paste0("(^nC|nF).*", y, "$")
        z <- grep(z, names(df), value=TRUE)
        fd <- pivot_longer(df, all_of(z))
        ms <- fd |>
            group_by(did, name) |>
            summarise_at("value", mean) |>
            mutate(y=ifelse(did == "raw", 1, 0.9))
        # upper quantile filtering
        if (y != "negprobes")
            fd <- mutate(group_by(fd, did, name),
                th=quantile(value, 0.99), ol=value > th,
                value=case_when(ol ~ th, TRUE ~ value))
        ggplot(fd, 
            aes(value, after_stat(ncount), fill=did)) +
            facet_wrap(~name, ncol=1, scales="free_x") +
            geom_histogram(alpha=0.5,
                position=position_dodge(w),
                col=NA, bins=nb, binwidth=bw) +
            geom_vline(
                linewidth=0.2, show.legend=FALSE,
                data=ms, aes(xintercept=value, col=did)) +
            geom_text(
                hjust=1, vjust=1, size=2, x=Inf, show.legend=FALSE,
                data=ms, aes(label=round(value), y=y, col=did)) +
            guides(fill=guide_legend(override.aes=list(alpha=1))) +
            scale_color_manual(NULL, values=c("orange", "blue")) + 
            scale_fill_manual(NULL, values=c("orange", "blue")) + 
            scale_y_continuous(n.breaks=2) + if (y == "RNA")
            labs(y="normalized count") else ylab(NULL) 
    })
    n <- format(ncol(x), big.mark=",")
    m <- format(length(i), big.mark=",")
    t <- bquote(bold(.(id))~"(N ="~.(n)*"; "*.(m)*")")
    wrap_plots(ps, nrow=1) +
        plot_layout(guides="collect") &
        plot_annotation(title=t) &
        .theme_black & theme(
            aspect.ratio=2/3,
            legend.position="bottom",
            axis.title.x=element_blank(),
            panel.grid.major=element_line(color="grey")) 
}

.plt_qc_bs <- \(x, md, id="") {
    # dependencies
    library(zoo)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggrastr)
    # setup
    xy <- data.frame(
        x=x$CenterX_local_px,
        y=x$CenterY_local_px)
    ds <- with(xy, cbind(
        top=max(y)-y, bottom=y-min(y),
        right=max(x)-x, left=x-min(x)))
    n <- format(ncol(x), big.mark=",")
    th_d <- md[[grep("th_d", names(md))]]
    th_n <- md[[grep("th_n", names(md))]]
    m <- format(sum(apply(ds > th_d, 1, all)), big.mark=",")
    t <- bquote(bold(.(id))~"(N ="~.(n)*"; "*.(m)*")")
    # wrangling
    cd <- pivot_longer(data.frame(colData(x), ds), all_of(colnames(ds)))
    dc <- cd[with(cd, value > th_d & nFeature_RNA > th_n), ]
    df <- bind_rows(list(raw=cd, fil=dc), .id="did")
    df$did <- factor(df$did, c("raw", "fil"))
    # averaging
    mu <- summarize(group_by(df, did), y=mean(nFeature_RNA))
    fd <- group_by(df, did, name, value) |>
        summarise_at("nFeature_RNA", mean) |>
        arrange(value) |> group_by(name) |> mutate(
            mu_x=rollmean(value, 10, mean, align="right", fill=NA),
            mu_y=rollmean(nFeature_RNA, 10, mean, align="right", fill=NA))
    # plotting
    ggplot(df, 
        aes(value, nFeature_RNA)) + 
        facet_grid(did~name) +
        geom_point_rast(
            shape=16, stroke=0, col="blue",
            size=1e5/ncol(x)*0.1, alpha=0.05+1e4/ncol(x)) +
        labs(x="distance to FOV border (px)") +
        coord_cartesian(xlim=c(0, 400), ylim=c(0, 2*max(mu$y))) +
        geom_line(data=fd, aes(mu_x, mu_y), col="navy", linewidth=0.2) +
        geom_hline(data=mu, aes(yintercept=y), col="orange", linewidth=0.2) +
        geom_hline(yintercept=md$th_ndet, col="orange", linewidth=0.2, lty=2) +
        geom_vline(xintercept=md$th_dist, col="orange", linewidth=0.2, lty=2) +
        .theme_black + theme(aspect.ratio=2/3) + ggtitle(t)
}

# plot dimensionality reduction
.plt_dr <- \(x, c, id="", dr="UMAP", s=0.2, 
    t=c("n", "z", "q"), th=2.5, qs=0.01) {
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
    if (is.data.frame(x)) {
        df <- x
    } else {
        f <- \(.) inherits(tryCatch(as.vector(.), error=\(.) .), "error")
        i <- !sapply(as.list(colData(x)), f)
        df <- reducedDim(x, dr)[, c(1, 2)]
        if (dr == "PCA") colnames(df) <- paste0("PC", c(1, 2))
        df <- data.frame(colData(x)[i], df, check.names=FALSE)
        if (c %in% rownames(x)) df[[c]] <- logcounts(x)[c, ]
    }
    df <- df[sample(nrow(df)), ]
    if (dr == "PCA") dr <- "PC"
    xy <- paste0(dr, c(1, 2))
    if (is.numeric(df[[c]])) {
        aes <- switch(match.arg(t), 
            n={ # no transformation
                scale_color_gradientn(colors=c("navy", "red", "gold", "ivory"))
            },
            z={ # thresholded z-normalization
                df[[c]] <- .z(df[[c]], th)
                scale_color_gradient2(low="blue", mid="ivory", high="red")
            },
            q={ # upper/lower quantile scaling
                df[[c]] <- .q(df[[c]], qs)
                scale_color_gradientn(
                    limits=c(0, 1), breaks=c(0, 1),
                    colors=c("navy", "red", "gold", "ivory"))
            })
        df <- df[order(abs(df[[c]]), na.last=FALSE), ]
    } else {
        aes <- list(
            scale_color_manual(values=.pal),
            guides(col=guide_legend(ncol=1, override.aes=list(size=2))))
    }
    ggplot(df, aes(.data[[xy[1]]], .data[[xy[2]]], col=.data[[c]])) +
        geom_vline(xintercept=0, col="grey", linewidth=0.1) +
        geom_hline(yintercept=0, col="grey", linewidth=0.1) +
        geom_point_rast(shape=16, stroke=0, size=s) +
        theme(legend.key.size=unit(0.5, "lines")) +
        ggtitle(.lab(id, nrow(df))) + 
        aes + theme_void(6) + theme(
            aspect.ratio=1, 
            panel.grid=element_blank(),
            legend.key.size=unit(0.5, "lines")) 
}

# plot embedding of 'InSituType' likelihoods
.plt_pb <- \(x, id="", s=0.1, a=0.5) {
    set.seed(112358)
    library(ggplot2)
    library(ggrastr)
    library(InSituType)
    # wrangling
    y <- flightpath_layout(x$logliks, profiles=x$profiles)
    k <- names(p <- y$meanconfidence)
    p <- sprintf("bold(`%s`)(%s)", k, round(p, 2))
    fd <- data.frame(y$clustpos, k, p)
    df <- data.frame(y$cellpos, k=x$clust)
    df <- df[sample(nrow(df)), ]
    # plotting
    ggplot(df, aes(x, y, col=k)) + 
        ggtitle(.lab(id, nrow(df))) +
        scale_color_manual(NULL, values=.pal) +
        geom_vline(xintercept=0, col="grey", linewidth=0.1) +
        geom_hline(yintercept=0, col="grey", linewidth=0.1) +
        geom_point_rast(shape=16, stroke=0, size=s, alpha=a) +
        geom_text(data=fd, aes(label=p), parse=TRUE, size=2, col="white") +
        .theme_black_void + theme(
            axis.title=element_text(hjust=0.5),
            aspect.ratio=1, legend.position="none",
            panel.background=element_rect(fill="black"))
}

# 'flightpath plot', i.e., embedding of 
# 'InSituType' assignment likelihoods
.plt_fp <- \(x, id="") {
    library(dplyr)
    library(ggplot2)
    library(InSituType)
    ks <- x$clust
    ll <- x$logliks
    ps <- x$profiles
    set.seed(194849)
    f <- \(.) sample(., min(1e4, length(.)))
    i <- unlist(lapply(split(seq(nrow(ll)), ks), f))
    y <- flightpath_layout(logliks=ll[i, ], profiles=ps)
    df <- data.frame(y$cellpos, k=ks[i])[sample(i), ]
    fd <- summarize(group_by(df, k), across(c("x", "y"), median))
    ggplot(df, aes(x, y, col=k)) + 
        .thm_xy_d(0.1) + theme(aspect.ratio=1, legend.position="none") +
        geom_text(data=fd, aes(label=k), size=1.5, col="black") +
        scale_color_manual(values=.pal_kid) +
        ggtitle(.lab(id, sum(!is.na(ks))))
}

# plot cluster composition by field of view
.plt_fq <- \(z, x, y, id="", hc=TRUE, h=FALSE) {
    library(ggplot2)
    library(SingleCellExperiment)
    # tabulate cell counts
    if (is(z, "SingleCellExperiment"))
        z <- data.frame(colData(z))
    ns <- table(z[[x]], z[[y]])
    ys <- sort(unique(z[[y]]))
    df <- as.data.frame(ns)
    # hierarchical clustering
    xo <- if (hc) {
        hc <- hclust(dist(prop.table(ns, 1)))
        hc$labels[hc$order]
    } else rownames(ns)
    # plotting
    aes <- if (h) {
        coord_flip(expand=FALSE)
    } else {
        list(coord_cartesian(expand=FALSE), theme(
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)))
    }
    ggplot(df, aes(Var1, Freq, fill=Var2)) +
        geom_col(
            position="fill", col="white", 
            linewidth=0.1, width=1, key_glyph="point") +
        scale_fill_manual(NULL, values=.pal_kid) +
        labs(x=x, y=NULL) +
        ggtitle(.lab(id, nrow(z))) +
        scale_x_discrete(limits=xo) +
        scale_y_continuous(n.breaks=2) +
        .thm_fig_d("minimal", "f") + aes + 
        theme(axis.ticks=element_blank()) 
}

# gene x cluster heatmaps including look-up, joint & split markers
.plt_de <- \(x, id="", z=TRUE, gs=NULL, n1=30, n2=80) {
    library(dplyr)
    library(tidyr)
    library(scran)
    library(scuttle)
    library(ggplot2)
    # DGE analysis
    ids <- factor(. <- x$clust, sort(unique(.)))
    nk <- length(names(ks) <- ks <- names(ns <- table(ids)))
    # selection 
    es <- x$profiles
    es <- normalizeCounts(es)
    fcs <- lapply(ks, \(i) {
        j <- setdiff(ks, i)
        a <- es[, i]
        b <- rowMeans(es[, j])
        c <- (c <- a/b)[!is.na(c) & is.finite(c)]
        c*rowDiffs(rowRanges(es))[names(c), ]
    })
    gs <- intersect(gs, rownames(es))
    sel <- \(n) lapply(ks, \(k) names(tail(sort(fcs[[k]]), n)))
    top <- unique(c(gs, unlist(sel(n1)))); all <- sel(n2)
    # aesthetics
    pal <- if (z) {
        es <- t(apply(es, 1, .z))
        scale_fill_gradientn(
            "z-scaled\nmean expr.", 
            colors=hcl.colors(9, "Blue-Red 3"),
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) 
    } else {
        scale_fill_gradientn(
            "mean\nexpr.", limits=c(0, NA), 
            colors=c("navy", "red", "gold"))
    }
    aes <- list(
        coord_equal(4/3, expand=FALSE),
        geom_tile(), theme_bw(6), theme(
            legend.position="bottom",
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            plot.background=element_blank(),
            legend.key=element_blank(),
            legend.background=element_blank(),
            legend.title=element_text(vjust=1.5),
            legend.key.width=unit(0.8, "lines"),
            legend.key.height=unit(0.4, "lines"),
            axis.text.y=element_text(size=4),
            axis.text.x=element_text(size=3, 
                angle=90, hjust=1, vjust=0.5)))
    # hierarchical clustering
    .x <- \(.) rownames(.)[hclust(dist(.))$order]
    .y <- \(.) colnames(.)[hclust(dist(t(.)))$order]
    df <- data.frame(gene=rownames(es), es, check.names=FALSE)
    df <- pivot_longer(df, -gene, names_to="k", values_to="y") 
    # plotting
    p1 <- ggplot(
        df[df$gene %in% top, ],
        aes(gene, k, fill=y)) + 
        scale_x_discrete(limits=.x(es[top, ])) + 
        scale_y_discrete(limits=.y(es[top, ])) +
        ggtitle(.lab(id, length(x$clust))) +
        aes + pal 
    ps <- lapply(ks, \(k) ggplot(
        df[df$gene %in% (gs <- all[[k]]), ], 
        aes(gene, k, fill=y)) + 
        scale_x_discrete(limits=.x(es[gs, ])) + 
        scale_y_discrete(limits=.y(es[gs, ])) +
        ggtitle(.lab(k, ns[[k]])) +
        aes + pal)
    c(list(p1), ps)
}

.plt_es <- \(x, gs, by="ist", xo=NULL, yo=NULL, th=2.5, asp=1, flip=FALSE, n=TRUE) {
    # dependencies
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(scuttle)
    })
    pbs <- aggregateAcrossCells(x, 
        x[[by]], subset.row=unique(unlist(gs)),
        statistics="mean", use.assay.type="logcounts")
    mtx <- t(assay(pbs))
    df <- data.frame(
        colData(pbs), mtx,
        check.names=FALSE) |>
        pivot_longer(rownames(pbs)) |>
        group_by(name, !!sym(by)) |>
        summarise_at("value", mean) |>
        mutate_at("value", .z, th=th)
    y <- pivot_wider(df)
    z <- as.matrix(y[, -1])
    rownames(z) <- y[[1]]
    xo <- if (is.null(xo) & is.list(gs)) names(gs) else .xo(z)
    yo <- if (is.null(yo) & is.list(gs)) rownames(pbs) else .yo(z)
    ns <- data.frame(n=format(pbs$ncells, big.mark=",")); ns[[by]] <- colnames(pbs)
    ns <- format(pbs$ncells, big.mark=",")
    
    gg <- if (flip) {
        ggplot(df, aes(name, .data[[by]], fill=value)) +
            scale_x_discrete(limits=yo) +
            scale_y_discrete(limits=xo) +
            geom_tile() + if (n) annotate("text", 
                nrow(pbs), colnames(pbs), 
                hjust=-0.5, size=1, label=ns)
    } else {
        ggplot(df, aes(.data[[by]], name, fill=value)) +
            scale_x_discrete(limits=xo) +
            scale_y_discrete(limits=yo) + 
            geom_tile() + if (n) annotate("text", 
                colnames(pbs), nrow(pbs)+1,
                hjust=-0.5, size=1, label=ns)
    }
    gg + scale_fill_gradient2(
        "z-scaled\nmean expr.", 
        low="blue", mid="ivory", high="red",
        limits=c(-th, th), breaks=c(-th+0.5, 0, th-0.5)) +
        coord_equal(asp, expand=FALSE, clip=ifelse(n, "off", "on")) +
        .theme_white_void + theme(
            legend.position="bottom",
            legend.text=element_text(size=3),
            legend.title=element_text(size=3, vjust=1),
            legend.key.width=unit(0.4, "lines"),
            legend.key.height=unit(0.2, "lines"),
            axis.text.y=element_text(size=3, hjust=1),
            axis.text.x=element_text(size=3/asp, angle=90, hjust=1, vjust=0.5)) +
        if (n) theme(plot.margin=margin(r=1, unit="lines")) 
}

# spatial plot
.plt_xy <- \(x, k, id="", s=NULL, split=TRUE) {
    # dependencies
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
    # wrangling
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(x)))
    names(colData(x))[xy] <- c("x", "y")
    if (length(names(k))) {
        cs <- match(colnames(x), names(k))
        nk <- length(ks <- levels(k <- as.factor(k[cs])))
        pal <- if (nk == 3) .pal_sub else .pal_kid
        pal <- setNames(pal[seq_len(nk)], ks)
    }
    # aesthetics
    df <- data.frame(colData(x), k)
    dx <- diff(range(df$x))
    dy <- diff(range(df$y))
    pt <- if (is.null(s)) min(dx, dy)/100/4 else s
    # plotting
    if (is.factor(df$k)) {
        fd <- df[!is.na(df$k), ]
        p0 <- ggplot(fd, aes(x, y, col=k)) + .thm_xy_d(pt) +
            scale_color_manual(NULL, drop=FALSE, values=pal) +
            ggtitle(.lab(id, nrow(fd)))
        ps <- if (split) lapply(c(ks, NA), \(k) {
            df$. <- if (is.na(k)) is.na(df$k) else grepl(sprintf("^%s$", k), df$k)
            ggplot(df[order(df$.), ], aes(x, y, col=.)) + 
                .thm_xy_d(pt) + theme(legend.position="none") +
                scale_color_manual(NULL, values=c("lavender", "purple")) +
                ggtitle(.lab(k, sum(df$.)))
        })
        c(list(p0), ps)
    } else {
        ggplot(df, aes(x, y, col=k)) + .thm_xy_c(pt) +
            scale_color_gradientn(colors=pals::jet()) +
            ggtitle(.lab(id, nrow(df)))
    }
}

.plt_rgb <- \(x, id) {
    # dependencies
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
    # wrangling
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(x)))
    names(colData(x))[xy] <- c("x", "y")   
    y <- reducedDim(x, "PCA")[, seq_len(3)]
    z <- sweep(y, 1, rowMins(y), `-`)
    z <- sweep(z, 1, rowMaxs(z), `/`)
    z <- apply(z, 1, \(.) rgb(.[1], .[2], .[3]))
    df <- data.frame(colData(x), y, z)
    # aesthetics
    dx <- diff(range(df$x))
    dy <- diff(range(df$y))
    pt <- min(dx, dy)/100/4
    # plotting
    p0 <- ggplot(df, aes(x, y, col=z)) + 
        ggtitle(.lab(id, nrow(df))) +
        scale_color_identity() + 
        .thm_xy_d(pt)
    ps <- lapply(colnames(y), \(.) {
        ggplot(df, aes(x, y, col=.q(.data[[.]]))) + 
            .thm_xy_c(pt) +
            scale_color_gradientn(
                paste0("q-scaled\n", ., " value"), 
                n.breaks=6, colors=pals::jet()) 
    })
    c(list(p0), ps)
}

# df = `arrow::Table` containing cell boundaries
# sce = corresponding 'SingleCellExperiment'
# c = character string; feature name or 'colData' to color points by
# t = "n"(o transformation), "z"(-normalization), or "q"(uantile) scaling
# th = scalar numeric; threshold to use when 't == "z"'
# qs = scalar or length-2 numeric; quantiles to use when 't == "q"'
# hl = logical/character vector; cells to highlight (others are 'blacked out')
.plt_ps <- \(df, sce=NULL, c="white", 
    t=c("n", "z", "q"), th=2.5, qs=0.01, 
    hl=NULL, lw=0.1, lc="black", id="") {
    library(dplyr)
    library(ggplot2)
    # filtering
    df <- df |>
        filter(cell %in% colnames(sce)) |>
        mutate(
            x=.px2mm(x_global_px), 
            y=.px2mm(y_global_px))
    i <- pull(df, "cell", as_vector=TRUE)
    i <- match(i, colnames(sce))
    j <- setdiff(names(colData(sce)), names(df))
    df <- cbind(as.data.frame(df), colData(sce)[i, j])
    if (c %in% rownames(sce)) {
        df[[c]] <- logcounts(sce)[c, i]
        # continuous coloring
        pal <- switch(match.arg(t), 
            n={ # no transformation
                scale_fill_gradientn(colors=c("navy", "red", "gold", "ivory"))
            },
            z={ # thresholded z-normalization
                df[[c]] <- .z(df[[c]], th)
                scale_fill_gradient2(low="blue", mid="ivory", high="red")
            },
            q={ # lower/upper quantile scaling
                df[[c]] <- .q(df[[c]], qs)
                scale_fill_gradientn(
                    limits=c(0, 1), breaks=c(0, 1),
                    colors=c("navy", "red", "gold", "ivory"))
            })
            thm <- list(pal, theme(            
                legend.key.width=unit(1, "lines"),
                legend.key.height=unit(0.25, "lines")))
    } else if (c %in% names(df)) {
        if (!is.numeric(df[[c]])) {
        # discrete coloring
        df[[c]] <- droplevels(factor(df[[c]]))
        pal <- scale_fill_manual(
            breaks=levels(df[[c]]),
            values=.pal(levels(df[[c]])))
        thm <- list(pal, 
            theme(legend.key.size=unit(0.5, "lines")), 
            guides(fill=guide_legend(override.aes=list(
                size=3, stroke=0.2, shape=21))))
        } else thm <- NULL
    } else {
        df[[c <- "foo"]] <- c
        thm <- list(scale_fill_identity(NULL))
    }
    # highlighting
    if (!is.null(hl)) {
        if (is.logical(hl)) 
            hl <- colnames(sce)[hl]
        df[[c]][!df$cell %in% hl] <- NA
    }
    # plotting
    ggplot(df, aes(x, y, fill=.data[[c]], group=cell)) + 
        geom_polygon(col=lc, linewidth=lw, key_glyph="point") +
        ggtitle(.lab(id, length(unique(df$cell)))) +
        coord_equal(expand=FALSE) + 
        .theme_black_void + 
        thm + theme(
            legend.position="bottom",
            legend.title=element_blank(),
            panel.background=element_rect(fill="black"))
}

# px to mm conversion
.px2mm <- \(.) .*0.00012028

# align shape layer exported from napari with
# local/global (FOV/tissue) cell coordinates
# sce = 'SingleCellExperiment'
# id = character string; shape identifier
# dir = character string; directory housing shape files
# xy_global/local = pattern matching tissue/FoV coordinates
.align_shape <- \(sce, dir=".", 
    xy_global="Center._global_px", 
    xy_local="Center._local_px") {
    # dependencies
    library(dplyr)
    library(SingleCellExperiment)
    . <- c("md", "px", "mm")
    . <- paste0(., ".csv")
    . <- file.path(dir, .)
    # get constants
    df <- read.csv(.[1])
    y_delta_mm <- df$y_delta_mm
    x_delta_mm <- df$x_delta_mm
    y_px_mm_ratio <- df$y_px_mm_ratio
    x_px_mm_ratio <- df$x_px_mm_ratio
    fov_id <- df$chosen_FOV
    fov_px <- 4256
    # get coordinates
    xy_px <- read.csv(.[2])
    xy_mm <- read.csv(.[3])
    # select cells in reference FOV
    cd <- data.frame(colData(sce)) |> 
        dplyr::filter(fov == fov_id) |> 
        select(cell_id, 
            matches(xy_global), 
            matches(xy_local)) |> 
        slice_head(n=1)
    # cd <- mutate(cd,
    #     x_slide_px=x_slide_mm/x_px_mm_ratio,
    #     y_slide_px=y_slide_mm/x_px_mm_ratio)
    xy_global <- grep(xy_global, names(cd))
    xy_local <- grep(xy_local, names(cd))
    # compute location of the FOV on the slide in px via difference 
    # between location of cell in FOV and slide; subtract FOV height
    # from y-coordinates as origin is in the top-left corner
    fov_px_x <- cd[[xy_global[1]]]-cd[[xy_local[1]]]
    fov_px_y <- cd[[xy_global[2]]]-(fov_px-cd[[xy_local[2]]])
    # convert position of the FOV's top left 
    # corner within the slide from px to mm
    fov_mm_x <- fov_px_x*x_px_mm_ratio
    fov_mm_y <- fov_px_y*x_px_mm_ratio+fov_px*x_px_mm_ratio
    # convert shape coordinates to fit objects frame of reference
    xy <- dplyr::rename(
        mutate(xy_mm, y_vals=-y_vals), 
        y_global_mm="y_vals", x_global_mm="x_vals")
    # move shape into right location using the 
    # reference FOV's top-left corner as anchor
    xy$y_global_mm <- xy$y_global_mm-xy$y_global_mm[1]+fov_mm_y-y_delta_mm 
    xy$x_global_mm <- xy$x_global_mm-xy$x_global_mm[1]+fov_mm_x+x_delta_mm 
    xy$x_global_px <- xy$x_global_mm/0.0001203
    xy$y_global_px <- xy$y_global_mm/0.0001203
    return(xy)
}

# subset 'SingleCellExperiment' by ROI's xy-coordinates
.subset_shape <- \(se, df, xy="Center._global_px") {
    require(SummarizedExperiment, quiet=TRUE)
    xy <- grep(xy, names(colData(se)))
    se[, sp::point.in.polygon(
        se[[xy[1]]], se[[xy[2]]], 
        df$x_global_px, df$y_global_px) == 1]
}

# thm ----
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

# base figure theme
.thm_fig <- \(.="minimal") list(
    get(paste0("theme_", .))(6),
    theme(
        legend.key=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank(),
        legend.background=element_blank(),
        plot.title=element_text(hjust=0.5)))

# discrete coloring
.thm_fig_d <- \(., l=c("c", "f")) {
    aes <- switch(match.arg(l),
        c=list(alpha=2, shape=19, size=2),
        f=list(alpha=1, shape=21, stroke=0, col=NA, size=2))
    thm <- list(theme(
        legend.key.size=unit(0, "lines")),
        guides(col=guide_legend(override.aes=aes)),
        guides(fill=guide_legend(override.aes=aes)))
    c(.thm_fig(.), list(thm))
}

# continuous coloring
.thm_fig_c <- \(.) {
    thm <- theme(
        legend.key.width=unit(0.4, "lines"),
        legend.key.height=unit(0.8, "lines"))
    c(.thm_fig(.), list(thm))
}

# theme for spatial plots
.thm_xy <- \(s) list(
    ggrastr::geom_point_rast(shape=16, stroke=0, size=s, raster.dpi=600),
    scale_x_continuous(expand=expansion(0, 0.1)),
    scale_y_continuous(expand=expansion(0, 0.1)),
    coord_equal(), theme(
        plot.margin=margin(),
        plot.title=element_text(hjust=0.5),
        panel.background=element_rect(color="grey", fill=NA)))
.thm_xy_d <- \(s) c(.thm_fig_d("void"), .thm_xy(s))
.thm_xy_c <- \(s) c(.thm_fig_c("void"), .thm_xy(s))

.theme_w <- theme(
    panel.grid=element_blank(),
    strip.background=element_blank(),
    legend.key.size=unit(0.5, "lines"))

.theme_b <- .theme_w + theme(
    legend.key=element_rect(fill="black"),
    strip.text=element_text(color="white"),
    plot.title=element_text(color="white"),
    legend.text=element_text(color="white"),
    legend.title=element_text(color="white"),
    plot.background=element_rect(fill="black", color="black"),
    panel.background=element_rect(fill="white", color="black"),
    legend.background=element_rect(fill="black", color="black"))
  
.theme_white <- 
    theme_linedraw(6) + .theme_w + theme(
        strip.text=element_text(color="black"))

.theme_black <- 
    theme_linedraw(6) + .theme_b + theme(
        axis.text=element_text(color="white"),
        axis.ticks=element_line(color="white"),
        axis.title=element_text(color="white"))

.theme_white_void <- theme_void(6) + .theme_w
.theme_black_void <- theme_void(6) + .theme_b

.xo <- \(.) rownames(.)[hclust(dist(.))$order]
.yo <- \(.) colnames(.)[hclust(dist(t(.)))$order]

.theme_fig_void <- theme(
    panel.grid=element_blank(),
    plot.margin=margin(2,2,2,2),
    plot.background=element_blank(),
    panel.background=element_blank(),
    legend.key=element_blank(),
    legend.background=element_blank(),
    legend.key.size=unit(0.2, "lines"),
    plot.title=element_text(size=4, color="black"),
    legend.text=element_text(size=3, color="black"),
    legend.title=element_text(size=3, color="black"),
    strip.text=element_text(size=4, color="black", margin=margin(1,1,1,1)))

.theme_fig <- .theme_fig_void + theme(
    axis.text=element_text(size=3, color="black"),
    axis.title=element_text(size=4, color="black"))
