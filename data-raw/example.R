# Builds the small synthetic example data shipped in inst/extdata.
# Run once from the package root:  source("example.R")
# Nothing here is real data. Counts are made up so the signals are known.
#
# Tree (written the way IQ-TREE writes it, outgroup in a basal split):
#   ((sp1,sp2),sp3) and (sp4,sp5) are the two main groups
# Planted gene flow (extra shared sites):
#   sp2 <-> sp5  (+250)     sp3 <-> sp4  (+150)

library(ape)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

tree_text <- "(((sp1,sp2),sp3),(sp4,sp5),outgroup);"
writeLines(tree_text, "inst/extdata/example_tree.treefile")

# the tree with unit branches, used only to find each trio's closest pair
tr <- compute.brlen(root(read.tree(text = tree_text), "outgroup", resolve.root = TRUE), 1)
dep <- node.depth.edgelength(tr)
M   <- mrca(tr)

key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))
extra <- c(key("sp2", "sp5"), key("sp3", "sp4"))
extra_n <- c(250, 150)

# shared-site count for a pair inside a given trio
pair_count <- function(a, b, trio) {
  others <- setdiff(trio, c(a, b))
  closest <- combn(trio, 2, function(p) dep[M[p[1], p[2]]])
  pairs   <- combn(trio, 2, function(p) key(p[1], p[2]))
  base <- if (pairs[which.max(closest)] == key(a, b)) 1000 else 300
  base + sum(extra_n[extra == key(a, b)])
}

rows <- list()
for (trio in combn(paste0("sp", 1:5), 3, simplify = FALSE)) {
  for (h3 in trio) {
    h12 <- sort(setdiff(trio, h3))
    abba <- pair_count(h12[2], h3, trio)   # H2 and H3 share
    baba <- pair_count(h12[1], h3, trio)   # H1 and H3 share
    D  <- (abba - baba) / (abba + baba)
    se <- sqrt((1 - D^2) / (abba + baba))
    rows[[length(rows) + 1]] <- data.frame(
      H1 = h12[1], H2 = h12[2], H3 = h3, nABBA = abba, nBABA = baba,
      Dstat = signif(D, 7), jackEst = signif(D, 7), SE = signif(se, 7),
      Z = signif(D / se, 7))
  }
}
tab <- do.call(rbind, rows)
write.table(tab, "inst/extdata/example_dstats.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("wrote", nrow(tab), "rows\n")