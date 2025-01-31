# dependencies
suppressPackageStartupMessages({
    
})

h5ls(h5 <- "imgs/11/targets.hdf5")
fs <- h5read(h5, "table/columns/fov")
length(cs <- which(fs$data == 1))

gs <- h5read(h5, "/table/columns/target")
i <- gs$dictionary$indices+1
i <- mapply(i=i[-length(i)], j=i[-1]-1, \(i, j)
    paste(gs$dictionary$data[seq(i, j)], collapse=""))
j <- i[match(1+gs$indices$data, seq_along(i))]

df <- data.frame(
    gene_id=j[cs],
    cell_id=h5read(h5, "/table/columns/CellId/data", index=list(cs)),
    x=h5read(h5, "/table/columns/x/data", index=list(cs)),
    y=h5read(h5, "/table/columns/y/data", index=list(cs)),
    z=h5read(h5, "/table/columns/z/data", index=list(cs)))

fd <- df[!grepl("Negative|System", df$gene_id), ]

at <- arrow::arrow_table(df)
write_parquet(at, "transcripts.parquet")