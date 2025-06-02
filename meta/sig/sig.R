setwd("~/projects/cosmx-wtx/meta/sig")
url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=%s&fileType=json"
.msig <- \(set) {
    names(set) <- set
    lapply(set, \(.) {
        url <- sprintf(url, .)
        tf <- tempfile(fileext=".json")
        download.file(url, tf, quiet=TRUE)
        jsonlite::fromJSON(tf)[[1]]$geneSymbols
    })
}
.save <- \(new, out) {
    if (!file.exists(out)) {
        jsonlite::write_json(new, out)
    } else {
        old <- jsonlite::fromJSON(out)
        if (!identical(old, new))
            jsonlite::write_json(new, out)
    }
}

# epi ----
gs <- .msig(read.delim("sig-epi.txt", header=FALSE)[[1]])
df <- readxl::read_xlsx("Socias2022.xlsx", skip=3)
gs <- c(gs, list(HRC_CORE=df$gene[df$core_HRC == 1]))
.save(gs, "sig-epi.json")

# imm ----
gs <- .msig(read.delim("sig-imm.txt", header=FALSE)[[1]])
.save(gs, "sig-imm.json")

# str ----
gs <- .msig(read.delim("sig-str.txt", header=FALSE)[[1]])
.save(gs, "sig-str.json")
