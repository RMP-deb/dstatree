# Copyright 2026 Ruben Muller-Perol 
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# D-statistic filter, object oriented version (one file, base R classes)
# One class (DStatFilter) holds the tree and the table, and the methods
# do the steps in order: loadTree -> loadTable -> filter -> sisterAdmixture

# ---------------------------------------------------------------
# Helper function (not part of the class)
# Makes the same label for a pair whichever way round it is given
# e.g. key("B", "A") and key("A", "B") both give "A B"
# ---------------------------------------------------------------
pairKey <- function(a, b) {
  ifelse(a < b, paste(a, b), paste(b, a))
}


# ---------------------------------------------------------------
# The class
# ---------------------------------------------------------------
DStatFilter <- setRefClass("DStatFilter",

  # Fields = the data each object stores
  fields = list(
    tree     = "ANY",         # rooted tree (ape phylo object)
    outgroup = "character",   # name of the outgroup tip
    depth    = "matrix",      # depth of the common ancestor for every pair of tips
    data     = "data.frame",  # the D-statistic table as loaded
    result   = "data.frame"   # the table after filtering
  ),

  # Methods = the things an object can do
  methods = list(

    # Read the tree, root it, and work out the depth lookup table
    loadTree = function(file, outgroup_name) {
      tr <- read.tree(file)

      # check the outgroup name is really in the tree
      if (!(outgroup_name %in% tr$tip.label)) {
        stop("Outgroup '", outgroup_name, "' is not a tip in the tree")
      }

      # root on the outgroup
      tr <- root(tr, outgroup = outgroup_name, resolve.root = TRUE)

      # give every branch length 1, so depth = number of steps from the root
      # (bigger depth = more recent)
      tr <- compute.brlen(tr, 1)

      # matrix of common ancestors for every pair of tips, worked out once
      m <- mrca(tr)
      node_depth <- node.depth.edgelength(tr)

      # save into the object's fields (<<- changes the field, not a local copy)
      tree     <<- tr
      outgroup <<- outgroup_name
      depth    <<- matrix(node_depth[m], nrow(m), dimnames = dimnames(m))

      invisible(.self)   # returning the object lets you chain methods
    },

    # Read the D-statistic table
    loadTable = function(file) {
      d <- read.table(file, header = TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE)

      # rename Dsuite style column names to H1 H2 H3 Z
      names(d)[names(d) == "P1"]      <- "H1"
      names(d)[names(d) == "P2"]      <- "H2"
      names(d)[names(d) == "P3"]      <- "H3"
      names(d)[names(d) == "Z-score"] <- "Z"

      # check the columns we need are there
      if (!all(c("H1", "H2", "H3", "Z") %in% names(d))) {
        stop("Table needs columns H1, H2, H3 and Z")
      }

      data <<- d
      invisible(.self)
    },

    # Keep only rows where H1 and H2 are closer to each other than to H3
    filter = function(z_threshold = 3) {
      if (nrow(depth) == 0) stop("Run loadTree() first")
      d <- data

      # drop any rows that contain the outgroup
      has_outgroup <- d$H1 == outgroup | d$H2 == outgroup | d$H3 == outgroup
      d <- d[!has_outgroup, ]

      # check every species name in the table is a tip in the tree
      table_names <- unique(c(d$H1, d$H2, d$H3))
      if (!all(table_names %in% rownames(depth))) {
        bad <- table_names[!(table_names %in% rownames(depth))]
        stop("Names not in tree: ", paste(head(bad, 5), collapse = ", "))
      }

      # depth of the common ancestor for each pair in every row
      d12 <- depth[cbind(d$H1, d$H2)]
      d13 <- depth[cbind(d$H1, d$H3)]
      d23 <- depth[cbind(d$H2, d$H3)]

      # congruent = H1 and H2 share a more recent ancestor than either does with H3
      # (H1 != H2 stops a species being compared with itself)
      d$congruent <- d12 > d13 & d12 > d23 & d$H1 != d$H2
      d <- d[d$congruent, ]

      # significance: p value from Z, then correct for many tests (BH)
      d$abs_Z       <- abs(d$Z)
      d$p           <- 2 * pnorm(-d$abs_Z)
      d$p_adj       <- p.adjust(d$p, method = "BH")
      d$significant <- d$abs_Z > z_threshold

      # strongest signals first
      d <- d[order(-d$abs_Z), ]
      rownames(d) <- NULL

      result <<- d
      invisible(.self)
    },

     # Look at significant rows and which split of the tree H1 and H2 sit on
    #   level = "cherry": only true sister pairs (clades of exactly 2 species)
    #   level = "clade":  any tree-congruent pair, annotated with its clade
    #   summary = TRUE:   one row per clade + H3 instead of one row per trio
    sisterAdmixture = function(z_threshold = 3, level = c("cherry", "clade"),
                               summary = FALSE) {
      level <- match.arg(level)
      if (nrow(result) == 0) stop("Run filter() first")
      r <- result

      ntip   <- length(tree$tip.label)
      tip_id <- seq_len(ntip)
      names(tip_id) <- tree$tip.label
      below  <- ape::prop.part(tree)   # species below each internal node

      # species names below any node (a tip is just itself)
      members <- function(n) {
        if (n <= ntip) tree$tip.label[n] else tree$tip.label[below[[n - ntip]]]
      }

      # the node where H1 and H2 meet = the split the pair sits on
      M    <- ape::mrca(tree)
      node <- M[cbind(r$H1, r$H2)]
      key  <- as.character(node)

      # for each split: its species, and the species under each of its branches
      split_nodes <- unique(key)
      clade_txt <- character(0)
      side_map  <- list()   # which branch each tip sits under
      side_txt  <- list()   # species list for each branch
      for (k in split_nodes) {
        n    <- as.integer(k)
        kids <- tree$edge[tree$edge[, 1] == n, 2]
        m    <- rep(NA_integer_, ntip)
        txt  <- character(length(kids))
        for (i in seq_along(kids)) {
          sp <- members(kids[i])
          m[tip_id[sp]] <- i
          txt[i] <- paste(sp, collapse = ",")
        }
        clade_txt[k]  <- paste(members(n), collapse = ",")
        side_map[[k]] <- m
        side_txt[[k]] <- txt
      }
      r$clade      <- unname(clade_txt[key])
      r$clade_size <- vapply(key, function(k) sum(!is.na(side_map[[k]])), numeric(1))
      r$clade_size <- unname(r$clade_size)

      # positive Z = H3 shares more with H2, negative Z = more with H1
      r$recipient <- ifelse(r$Z > 0, r$H2, r$H1)

      # which branch of the split the recipient is on
      r$toward_side <- NA_character_
      for (k in split_nodes) {
        i <- key == k
        r$toward_side[i] <- side_txt[[k]][side_map[[k]][tip_id[r$recipient[i]]]]
      }

      if (level == "cherry") r <- r[r$clade_size == 2, ]
      r$hit <- r$abs_Z > z_threshold

      # one row per significant trio
      if (!summary) {
        out <- r[r$hit, ]
        out$hit <- NULL
        out <- out[order(-out$abs_Z), ]
        rownames(out) <- NULL
        return(out)
      }

      # one row per clade + H3: how many pairs tested, how many significant,
      # and whether they all point at the same branch
      groups <- split(seq_len(nrow(r)), paste(r$clade, r$H3, sep = " || "))
      out <- do.call(rbind, lapply(groups, function(i) {
        g <- r[i, ]
        h <- g[g$hit, ]
        sides <- unique(h$toward_side)
        toward <- if (length(sides) == 0) NA_character_
                  else if (length(sides) == 1) sides
                  else "mixed"
        data.frame(clade = g$clade[1], H3 = g$H3[1], n_tests = nrow(g),
                   n_significant = nrow(h), toward_side = toward,
                   max_abs_Z = max(g$abs_Z), stringsAsFactors = FALSE)
      }))
      out <- out[order(-out$n_significant, -out$max_abs_Z), ]
      rownames(out) <- NULL
      return(out)
    },


    # Write the filtered table to a file
    saveResult = function(file) {
      write.csv(result, file, row.names = FALSE)
      cat("Saved", nrow(result), "rows to", file, "\n")
      invisible(.self)
    },

    # What you see when you type the object's name
    show = function() {
      cat("DStatFilter object\n")
      cat("  tips in tree:   ", nrow(depth), "\n")
      cat("  rows loaded:    ", nrow(data), "\n")
      cat("  rows after filter:", nrow(result), "\n")
    }
  )
)
