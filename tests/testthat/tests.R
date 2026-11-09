# Core tests for dstatree.


# Writes a tiny tree and table to temporary files and returns their paths.
# Tree: ((((A,B),C),D),O) once rooted on O. A+B is a cherry.
make_tiny <- function(rows = NULL, tree_text = "(((A,B),C),D,O);") {
  if (is.null(rows)) {
    rows <- data.frame(
      H1 = c("A", "A", "A", "A", "C"),
      H2 = c("B", "C", "B", "B", "D"),
      H3 = c("C", "B", "D", "O", "A"),
      Z  = c(4.2, 3.5, -3.8, 9, 2),
      stringsAsFactors = FALSE)
  }
  tree_file  <- tempfile(fileext = ".tre")
  table_file <- tempfile(fileext = ".txt")
  writeLines(tree_text, tree_file)
  write.table(rows, table_file, sep = "\t", quote = FALSE, row.names = FALSE)
  list(tree = tree_file, table = table_file)
}

# Three rows inside the clade A,B,C, used by the sister-pair tests.
sister_rows <- data.frame(
  H1 = c("A", "B", "A"),
  H2 = c("C", "C", "B"),
  H3 = c("D", "D", "C"),
  Z  = c(-5, -4, 4.2),
  stringsAsFactors = FALSE)

filtered <- function(rows = sister_rows) {
  p <- make_tiny(rows)
  DStatFilter$new()$loadTree(p$tree, "O")$loadTable(p$table)$filter()
}


test_that("only tree-congruent rows are kept and outgroup rows are dropped", {
  p <- make_tiny()
  f <- DStatFilter$new()$loadTree(p$tree, "O")$loadTable(p$table)$filter()
  # (A,B;C) and (A,B;D) fit the tree. (A,C;B) and (C,D;A) do not. The O row is dropped.
  expect_equal(nrow(f$result), 2)
  expect_true(all(f$result$H1 == "A" & f$result$H2 == "B"))
  expect_true(all(f$result$congruent))
})

test_that("clear errors for the common mistakes", {
  p <- make_tiny()
  expect_error(DStatFilter$new()$loadTree(p$tree, "nope"), "not a tip")

  bad <- make_tiny(data.frame(H1 = "A", H2 = "B", Z = 1))
  expect_error(DStatFilter$new()$loadTable(bad$table), "needs columns")

  wrong <- make_tiny(data.frame(H1 = "A", H2 = "B", H3 = "zzz", Z = 1,
                                stringsAsFactors = FALSE))
  f <- DStatFilter$new()$loadTree(wrong$tree, "O")$loadTable(wrong$table)
  expect_error(f$filter(), "Names not in tree")

  expect_error(DStatFilter$new()$filter(), "loadTree")
  g <- DStatFilter$new()$loadTree(p$tree, "O")$loadTable(p$table)
  expect_error(g$sisterAdmixture(), "filter")
})


test_that("cherry level only sees true sister pairs", {
  s <- filtered()$sisterAdmixture()
  expect_equal(nrow(s), 1)
  expect_equal(c(s$H1, s$H2, s$H3), c("A", "B", "C"))
  expect_equal(s$clade_size, 2)
})

test_that("clade level sees every congruent pair, with clade information", {
  s <- filtered()$sisterAdmixture(level = "clade")
  expect_equal(nrow(s), 3)
  expect_true(all(c("clade", "clade_size", "recipient", "toward_side") %in% names(s)))
})

test_that("positive Z points at H2 and negative Z at H1", {
  s <- filtered()$sisterAdmixture(level = "clade")
  expect_equal(s$recipient[s$Z > 0], s$H2[s$Z > 0])
  expect_equal(s$recipient[s$Z < 0], s$H1[s$Z < 0])
})

test_that("summary groups pairs that point at the same side", {
  s <- filtered()$sisterAdmixture(level = "clade", summary = TRUE)
  row <- s[s$H3 == "D", ]
  expect_equal(row$n_tests, 2)
  expect_equal(row$n_significant, 2)
  expect_setequal(strsplit(row$toward_side, ",")[[1]], c("A", "B"))
})

test_that("summary reports 'mixed' when pairs point at different sides", {
  mixed <- data.frame(H1 = c("A", "B"), H2 = c("C", "C"), H3 = c("D", "D"),
                      Z = c(-5, 5), stringsAsFactors = FALSE)
  s <- filtered(mixed)$sisterAdmixture(level = "clade", summary = TRUE)
  expect_equal(s$toward_side, "mixed")
})