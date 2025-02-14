p <- .plt_dr(readRDS(args[[1]]), wcs$kid, id=wcs$sid)
ggsave(args[[2]], p, width=8, height=7, units="cm")
