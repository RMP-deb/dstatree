# Basic checks using plain base R (no extra packages needed)
library(dstatree)

tree_file  <- tempfile(fileext = ".tre")
table_file <- tempfile(fileext = ".txt")
writeLines("(((A,B),C),D,O);", tree_file)
write.table(
  data.frame(H1 = c("A", "A", "A", "A", "C"), H2 = c("B", "C", "B", "B", "D"),
             H3 = c("C", "B", "D", "O", "A"), Z = c(4.2, 3.5, -3.8, 9, 2)),
  table_file, sep = "\t", quote = FALSE, row.names = FALSE)

f <- DStatFilter$new()$loadTree(tree_file, "O")$loadTable(table_file)$filter()

# only (A,B,C) and (A,B,D) agree with the tree; the outgroup row is dropped
stopifnot(nrow(f$result) == 2)
stopifnot(all(f$result$H1 == "A" & f$result$H2 == "B"))

# sister pair signals: positive Z -> H2 (B), negative Z -> H1 (A)
s <- f$sisterAdmixture()
stopifnot(nrow(s) == 2)
stopifnot(identical(s$recipient[s$H3 == "C"], "B"))
stopifnot(identical(s$recipient[s$H3 == "D"], "A"))

# essential error checks fire
stopifnot(inherits(try(DStatFilter$new()$loadTree(tree_file, "nope"), silent = TRUE),
                   "try-error"))
stopifnot(inherits(try(DStatFilter$new()$filter(), silent = TRUE), "try-error"))


# --- clade level and per-split summary ---------------------------------
tree2  <- tempfile(fileext = ".tre")
table2 <- tempfile(fileext = ".txt")
writeLines("(((A,B),C),D,O);", tree2)
write.table(
  data.frame(H1 = c("A", "B", "A"), H2 = c("C", "C", "B"),
             H3 = c("D", "D", "C"), Z = c(-5, -4, 4.2)),
  table2, sep = "\t", quote = FALSE, row.names = FALSE)
g <- DStatFilter$new()$loadTree(tree2, "O")$loadTable(table2)$filter()

# cherry level only sees (A,B;C); clade level sees all three congruent rows
stopifnot(nrow(g$sisterAdmixture()) == 1)
stopifnot(nrow(g$sisterAdmixture(level = "clade")) == 3)

# clade A,B,C with H3 = D: two pairs tested, both significant,
# both pointing at the same branch (the A,B side)
s <- g$sisterAdmixture(level = "clade", summary = TRUE)
row <- s[s$H3 == "D", ]
stopifnot(row$n_tests == 2, row$n_significant == 2)
stopifnot(setequal(strsplit(row$toward_side, ",")[[1]], c("A", "B")))