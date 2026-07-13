# spDNS

spDNS is an extension of scDNS for single-cell splicing-regulatory network
analysis.

The original scDNS framework quantifies condition-associated changes in
pairwise network relationships using joint and conditional probability
distributions, Jensen-Shannon divergence, network-based null models,
node-level Z-scores, and cell-specific point-wise divergence.

**spDNS retains the original scDNS network-divergence framework and adapts it
to integrated splicing-factor expression and splicing-event PSI data connected
by a predefined directed bipartite SF-to-event network.**

spDNS introduces three primary adaptations:

1. **Heterogeneous node integration:** source nodes contain splicing-factor
   expression and target nodes contain splicing-event PSI.
2. **Scale harmonization:** SF expression is transformed to a bounded 0-1
   scale, matching the numerical range of PSI.
3. **Directed bipartite topology:** observed and randomized networks preserve
   the SF-to-event direction, SF out-degree, and event in-degree.

spDNS does not reconstruct the SF-event network de novo. A predefined
SF-to-event network must be supplied. It does not perform PSI imputation,
event-coordinate matching, disease-specific preprocessing, GAT or Node2Vec
network reconstruction, graph contrastive learning, or automatic expansion
of an SF-to-gene network. It does not define its result as absolute SF
activity.

spDNS is derived from scDNS and retains its original network-divergence and
cell-contribution framework. The S4 class remains named `scDNS`, and the five
core function names are retained for compatibility.

## Installation

```r
devtools::install_github("ztian608-oss/spDNS")
library(spDNS)
```

## Input contract

The `data` matrix is node by cell:

```r
data_mat <- rbind(sf_data, event_psi)
```

- `sf_data`: SF by cell, with values in 0-1.
- `event_psi`: event by the same cells, with PSI retained in 0-1.
- SF expression and event PSI are harmonized to a common bounded numerical
  range of 0-1. This does not mean that they have identical statistical
  distributions.
- The network has exactly the directed interpretation `source SF -> target
  splicing event`. Source and target node sets must not overlap.

The inherited `counts` slot is assembled as:

```r
counts_mat <- rbind(sf_counts, event_support)
```

`sf_counts` contains raw SF expression counts. `event_support` should
preferably contain the real event denominator/support. If support is not
available, a matrix of ones may be used only as a structural compatibility
placeholder.

**A constant event-count matrix is only a compatibility placeholder for the
inherited scDNS object structure and does not represent actual event read
support.** Consequently, dropout/support adjustment for those event rows must
not be interpreted as measured event coverage.

`counts_mat` and `data_mat` must have identical row names, column names, and
ordering. `GroupLabel` must follow the cell-column order.

## Transforming SF expression

Before entering the five-step workflow, transform all cells participating in
the comparison together. Do not transform each condition separately.

```r
sf_data_log <- sweep(sf_counts, 2, Matrix::colSums(sf_counts), "/")
sf_data_log <- log1p(sf_data_log * 10000)

z_cdf01_row <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  if (!is.finite(mu) || !is.finite(sigma) || sigma == 0) {
    return(rep(0.5, length(x)))
  }
  y <- pnorm((x - mu) / sigma)
  y[!is.finite(y)] <- 0.5
  y
}

sf_data <- t(apply(sf_data_log, 1, z_cdf01_row))
```

The package uses the same single internal implementation,
`transform_sf_expression_01()`. It transforms each SF independently across all
compared cells, maps zero-variance SFs to 0.5, and preserves within-SF cell
ordering. PSI is not z-scored or CDF-transformed.

## Five-step sf_event workflow

```r
library(spDNS)

# sf_data: SF x cell, values in 0-1
# event_psi: event x cell, values in 0-1
data_mat <- rbind(sf_data, event_psi)

# event_support should contain measured support when available.
counts_mat <- rbind(sf_counts, event_support)

obj <- CreatScDNSobject(
  counts = counts_mat,
  data = data_mat,
  Network = sf_event_network,
  GroupLabel = group_vector,
  network.type = "sf_event"
)

obj <- scDNS_1_CalDivs(scDNSobject = obj)

obj <- scDNS_2_creatNEAModel_v2(
  scDNSobject = obj,
  n.randNet = 5000
)

obj <- scDNS_3_GeneZscore_v2(scDNSobject = obj)

obj <- scDNS_4_scContribution(
  scDNSobject = obj,
  sigGene = unique(sf_event_network$source),
  q.th = 1
)
```

The inherited joint-density, conditional-density, JSD/KLD, coarse-graining,
GeneVariability, dropout/support adjustment, degree normalization, random
network fitting, label shuffling, node Z-score, combined Z-score, P-value/FDR,
point-wise divergence, and single-cell contribution calculations are retained.

## Outputs and interpretation

- `obj@Network`: edge-level SF-event divergence results.
- `obj@Zscore`: source and target node-level perturbation scores, P values and
  ranks. `node_type` identifies `SF` and `splicing_event`; primary
  interpretation is restricted to source-side SFs.
- `obj@scZscore`: cell-specific SF perturbation contribution matrix. In
  `sf_event` mode it is SF by cell.

The SF-level score represents condition-associated perturbation of the
relationship between an SF and its connected splicing-event module. It is not
an SF expression fold change, absolute SF activity, activation/repression
score, or mutation probability. Event-level node scores are retained as
auxiliary output.

## Original gene_gene mode

The original gene-gene workflow remains the default:

```r
obj <- CreatScDNSobject(
  counts = gene_counts,
  data = gene_data,
  Network = gene_network,
  GroupLabel = group_vector,
  network.type = "gene_gene"
)
```

Omitting `network.type` is equivalent to `network.type = "gene_gene"`.

## Citation and license

spDNS is derived from the original
[scDNS repository](https://github.com/HChaoLab/scDNS). Please cite the
original work described by its authors as *scDNS: Characterizing Gene
Perturbations in Single Cells via Network Divergence Analysis* (2025), as well
as this software adaptation. See `NOTICE` and `inst/CITATION`.

The original scDNS repository permits academic and research use and restricts
commercial use without prior written permission. spDNS is distributed under
the same restrictions; see `LICENSE`.
