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
