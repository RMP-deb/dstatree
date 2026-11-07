# dstatree

Tree-aware filtering of D-statistic (ABBA-BABA) result tables. Keeps only the
trios whose topology agrees with a rooted species tree, adjusts for multiple
testing, and reports significant signals involving true sister pairs.

## Install

```r
# install.packages("remotes")
remotes::install_github("RMP-deb/dstatree")

# or from a downloaded source archive
install.packages("dstatree_0.0.1.tar.gz", repos = NULL, type = "source")
```

## Use

```r
library(dstatree)

f <- DStatFilter$new()
f$loadTree("alignment_chr6.fa.treefile", "macrovipera_chr6")
f$loadTable("vipera_dstat_jackknife.txt")
f$filter()
f                                  # summary
f$sisterAdmixture()                # significant signals with true sister pairs as H1/H2
f$saveResult("dstat_clean.csv")
```

## Limits

A D-statistic cannot test for gene flow between the two members of a sister
pair themselves, and it cannot give the direction of gene flow.

## Licence

Apache License 2.0. See `LICENSE`.
