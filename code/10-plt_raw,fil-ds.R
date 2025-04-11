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

# FOVs are 4256 x 4256 px; position
# coordinates  are top-left corners
pos <- "data/raw/%s/positions.csv.gz"
pos <- read.csv(sprintf(pos, raw$did[1]))
w <- h <- 4256
x <- raw$CenterX_global_px
y <- raw$CenterY_global_px
i <- match(raw$fov, pos$FOV)
d <- cbind(
    l=x-pos$x_global_px[i],
    t=pos$y_global_px[i]-y,
    r=pos$x_global_px[i]+w-x,
    b=y-pos$y_global_px[i]+h)
colData(raw)[colnames(d)] <- d

# wrangling
df <- data.frame(colData(raw)) |>
    pivot_longer(all_of(colnames(d))) |>
    mutate(name=factor(name, 
        c("t", "r", "b", "l"),
        c("top", "right", "bottom", "left"))) |>
    group_by(name, value) |>
    summarise_at("nCount_RNA", mean) |>
    arrange(value) |> group_by(name) |> mutate(
        mu_x=rollmean(value, 10, mean, align="left", fill=NA),
        mu_y=rollmean(nCount_RNA, 10, mean, align="left", fill=NA))
mu <- data.frame(name=unique(df["name"]), y=mean(raw$nCount_RNA))

# aesthetics
dy <- ceiling(max(mu$y)/500)*500
n <- ncol(raw); m <- ncol(fil)
n <- format(n, big.mark=",")
m <- format(m, big.mark=",")
md <- metadata(fil)

# plotting
gg <- ggplot(df, aes(mu_x, mu_y)) + 
    facet_grid(~name) + geom_line(col="navy", linewidth=0.2) +
    geom_vline(xintercept=md$th_d, col="red", linewidth=0.2, lty=2) +
    geom_hline(yintercept=md$th_n, col="red", linewidth=0.2, lty=2) +
    geom_hline(data=mu, aes(yintercept=y), col="blue", linewidth=0.2) +
    ggtitle(bquote(bold(.(wcs$sid))~"(N ="~.(n)*";"~.(m)*")")) +
    labs(x="distance to FOV border (px)", y="counts (RNA)") +
    coord_cartesian(xlim=c(0, 400), ylim=c(0, dy)) +
    theme_bw(6) + theme(
        aspect.ratio=2/3,
        panel.grid=element_blank(),
        plot.background=element_blank())

# saving
ggsave(args[[3]], gg, units="cm", width=8, height=2.75)
