# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# restrict to epithelia
ist <- grep("epi", args[[2]], value=TRUE)

df <- mapply(
    x=args[[1]], y=ist,
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- readRDS(y)$clust
        sid <- paste(sce$sid[1])
        # wrangling
        idx <- intersect(colnames(sce), names(ist))
        typ <- setNames(sce$typ, colnames(sce))[idx]
        ns <- table(typ=typ[idx], kid=ist[idx])
        ns <- as.data.frame(ns, responseName="n")
        data.frame(sid, ns)
    }) |> do.call(what=rbind)

# more wrangling
fd <- df |>
    group_by(typ) |>
    mutate(p=n/sum(n)) |>
    filter(grepl("CRC", typ)) |>
    mutate(kid=gsub("^epi\\.", "", kid)) 

# specify grading
gs <- list(G1=c(120, 210, 221), G2=110, G3=c(231, 232, 242))
gt <- setNames(rep.int(names(gs), sapply(gs, length)), unlist(gs))
table(fd$G <- factor(gt[match(fd$sid, names(gt))], names(gs)))

# spot check
mu <- fd |>
    group_by(G, kid) |>
    summarize_at("p", mean)

# get coefficients
fdd <- fd |>
    group_by(kid) |>
    mutate(p=p-mean(p))
fit <- lm(p*as.integer(G)~-1+kid, fdd)
round(sort(coef(fit)), 2)

# wrangling again
res <- data.frame(
    est=100*coef(fit),
    kid=sort(unique(fd$kid))) |>
    arrange(est) |>
    mutate(kid=factor(kid, unique(kid)))

# plotting
gg <- ggplot(res, aes(est, kid, fill=sign(est)*sqrt(abs(est)))) + geom_col() +
    geom_vline(xintercept=0, linewidth=0.2) +
    scale_fill_gradient2(low="navy", high="firebrick") +
    labs(x=bquote(Delta*"(%Y) per"~Delta*"G"), y=NULL) +
    .thm_fig_d("bw") + theme(
        aspect.ratio=1,
        legend.position="none", 
        axis.ticks.y=element_blank(),
        panel.grid.minor=element_blank()) 

# saving
ggsave(args[[3]], gg, units="cm", width=5, height=4)
