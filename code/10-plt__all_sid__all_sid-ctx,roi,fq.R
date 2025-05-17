args <- list(
   list.files("outs", "ctx-", full.names=TRUE),
   list.files("outs", "roi-", full.names=TRUE),
   "plts/ctx,roi,fq.pdf")
args <- lapply(args, \(.) grep("241", ., value=TRUE, invert=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(SingleCellExperiment)
})

# loading
df <- mapply(x=args[[1]], y=args[[2]], \(x, y) {
    df <- readRDS(x)
    se <- readRDS(y)
    cs <- match(colnames(se), df$cid)
    data.frame(df[cs, ], colData(se))
}, SIMPLIFY=FALSE) |> do.call(what=rbind)

# wrangling
df <- df |>
    mutate(tyq=gsub(".*(REF|TVA|CRC).*", "\\1", typ)) |>
    mutate(tyq=factor(tyq, c("REF", "TVA", "CRC"))) |>
    filter(!is.na(typ), !is.na(ctx))

# plotting
p0 <- .plt_fq(df, "ctx", "tyq", h=TRUE) + labs(x="niche") +
    scale_fill_manual("domain", values=.pal_roi) +
    theme(aspect.ratio=length(unique(df$ctx))/20)

ps <- by(df, df$tyq, \(fd) {
    .plt_fq(fd, "typ", "ctx", h=TRUE) + labs(x="region") +
        scale_fill_manual("niche", values=.pal_ctx) +
        theme(aspect.ratio=length(unique(fd$typ))/20)
})

gg <- wrap_plots(c(list(p0), ps), ncol=1, guides="collect") & 
    theme(
        plot.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(hjust=0))

# saving
ggsave(args[[3]], gg, units="cm", width=8, height=10)
