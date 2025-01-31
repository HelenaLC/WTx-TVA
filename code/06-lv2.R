ist <- lapply(args[[1]], readRDS)
sub <- gsub(".*(epi|imm|str).*", "\\1", basename(args[[1]]))

# jst <- list(
#     prob=unlist(lapply(ist, \(.) .$clust)),
#     clust=unlist(lapply(ist, \(.) .$clust)),
#     profiles=setNames(lapply(ist, \(.) .$profiles), sub),
#     sub=rep.int(sub, vapply(ist, \(.) length(.$clust), numeric(1))))
# 
# table(jst$clust, jst$sub)
# table(jst$sub)

ist <- setNames(ist, sub)
saveRDS(ist, args[[2]])