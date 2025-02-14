# loading
sce <- readRDS(args[[1]])
lab <- read.csv(args[[2]])

# labelling
i <- match(sce$kid_lv1, lab[[1]])
j <- factor(lab[[2]][i], unique(lab[[2]]))
table(sce$lab_lv1 <- j)

# saving
saveRDS(sce, args[[3]])
