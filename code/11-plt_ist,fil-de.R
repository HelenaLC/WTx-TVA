# loading
ist <- readRDS(args[[1]])

# plotting
ps <- .plt_de(ist, id=wcs$sid, z=TRUE)

# saving
h <- (2+length(unique(ist$clust))/5)
w <- length(unique(ps[[1]]$data$gene))/8
pdf(args[[3]], onefile=TRUE, width=w/2.54, height=h/2.54)
for (p in ps) print(p); dev.off()
