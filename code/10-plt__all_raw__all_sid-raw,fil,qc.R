# dependencies
suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
fd <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        raw <- readRDS(x)
        fil <- readRDS(y)
        assays(raw) <- assays(fil) <- list()
        fil <- colnames(raw) %in% colnames(fil)
        data.frame(colData(raw), fil)
    }) |> do.call(what=rbind)

# wrangling
lt <- c("counts", "features", "area")
xs <- c(lt, "counts/um2", "features/um2")
df <- fd |>
    mutate(area=Area, 
        counts=nCount_RNA, 
        features=nFeature_RNA) |>
    mutate(`counts/um2`=counts/area) |>
    mutate(`features/um2`=features/area) |>
    mutate(sid=factor(sid, unique(sid))) |>
    mutate(fil=factor(fil,
        levels=c(TRUE, FALSE),
        labels=c("kept", "removed"))) |>
    mutate(across(all_of(lt), log10)) |>
    pivot_longer(all_of(xs)) |>
    mutate(name=factor(name, xs))

# annotations
mu <- df |>
    group_by(sid, fil, name) |>
    summarise_at("value", median) |>
    mutate(value=round(value, 2)) |>
    mutate(value=format(value, big.mark=","))

# plotting
gg <- ggplot(df, aes(sid, value, fill=sid)) + 
    geom_boxplot(key_glyph="point", alpha=2/3, linewidth=0.2,
        outlier.shape=16, outlier.size=0.2, outlier.stroke=0) +
    scale_fill_manual(NULL, values=.pal_sid) +
    facet_grid(name~fil, scales="free_y") +
    geom_text(
        data=mu, aes(y=Inf, label=value),
        angle=90, size=1.8, hjust=1.2, col="darkgrey") +
    scale_y_continuous(limits=\(.) c(min(.), max(.)*1.2)) +
    .thm_fig_d("bw", "f") + theme(
        aspect.ratio=1,
        legend.position="none",
        axis.title=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))

# saving
ggsave(args[[3]], gg, unit="cm", width=8, height=12)
