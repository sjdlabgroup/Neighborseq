
<!-- README.md is generated from README.Rmd. Please edit that file -->
# Neighbor-seq

<!-- badges: start -->
<!-- badges: end -->
The spatial context of cells in a tissue and their resulting cell-cell communications influence numerous processes, including cellular differentiation, organ development and homeostasis, and immune interactions in disease. Neighbor-seq is an R package designed to infer the architecture of direct cell-cell interactions from standard single-cell RNA sequencing (scRNA-seq) data. Cell aggregates (multiplets) naturally arise in scRNA-seq experiments when two or more cells are captured in the same reaction droplet, and they typically represent at least several percent of all capture events. Neighbor-seq reconstructs physical cell-cell interactions by identifying, annotating, and analyzing cell multiplets from the undissociated cell fractions in scRNA-seq data.

![Figure 1A. A schematic representation of the Neighbor-seq workflow](Figure%201A.png)

## Installation

You can install the released version of Neighbor-seq from Github with:

``` r
devtools::install_github('sjdlabgroup/Neighborseq')
# this might take long if you don't already have Seurat installed.
```

## Neighbor-seq workflow

Neighbor-seq is a method to infer physical cell-cell communications by identifying, annotating, and analyzing cell multiplets from the undissociated cell fractions in scRNA-seq data using machine learning approaches. The Neighbor-seq algorithm consists of the following components, each further described below: (i) barcode clustering and marker gene identification, (ii) Random Forest classifier training to identify multiplets and their cell type compositions, (iii) calculating enrichment scores for cell-cell interactions, and (iv) construction of cell-cell interactome network and analysis of cell-neighbor transcriptomes. Each componenet can be programmed individually (described first) or via a few wrapper functions (described later).

This demonstration uses a subsample of small intestine scRNA-seq data from Haber et al., Nature 2017.

### Step 1: preparing the input data

If cell-type labels are not provided, Neighbor-seq utilizes a wrapper function to run Seurat functions that TP10K normalize and scale the scRNA-seq data, perform dimensionality reduction, identify cell type clusters, and find cell-type marker genes. If cell-type labels are known a priori, only normalization and marker finding functions are run. The primary outputs are (1) a gene by cell counts matrix keeping only the union of the top 50 marker genes for each cell-type, and (2) a vector of cell type identities. These steps are all implemented in a single R function.

``` r
library(Neighborseq)
library(RColorBrewer)

data = read.table('si.data.txt')
meta.data = read.table('metadata.txt', header = T)
ns.data = prep_cell_mat(ge = data, celltypes = meta.data$Cluster, logfc.threshold = 0.5)
#> Calculating cluster Goblet
#> Calculating cluster Enterocyte
#> Calculating cluster Stem
#> Calculating cluster TA
#> Calculating cluster Paneth
#> Calculating cluster EP
#> Calculating cluster Tuft
#> Calculating cluster Enteroendocrine
```

### Step 2: creating artificial multiplets

Next, all homotypic and heterotypic combinations for a default of 2 cell-types are enumerated (e.g. AA, AB, etc.). For each neighbor-type, artificial multiplets are created by randomly sampling cells from the constituent cell types and the prepared input gene by cell matrix and summing their gene counts. A default of 100 artificial multiplets is created for each neighbor-type, and these are combined with a default of 100 singlets from each cell-type.

``` r
set.seed(0) # for reproducibility
mt = multiplet_types(meta.data$Cluster)
am = artificial_multiplets(cell.mat = ns.data$cell.mat, 
                           celltypes = ns.data$celltypes, 
                           multiplet_classes = mt)
```

### Step 3: training a random forest classifier

A default of 80% of the artificial multiplet dataset set is used to train a random forest that takes as input the gene counts and predicts the barcode singlet or multiplet composition. Classification performance can be assessed using the hold-out data and plotting multivariate receiever-operator curves.

``` r
rf = multiplet_rf(am)
#> [1]  train-mlogloss:3.545750
mroc = mroc_format(rf$test, rf$pred)
mroc_plot(mroc)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="35%" height="35%" style="display: block; margin: auto;" />

### Step 4: predicting barcode type, assessing multiplet enrichment, and plotting the network

The trained random forest is then used to predict the barcode composition of all barcodes in the original dataset. Neighbor-seq calculates an enrichment score for each cell-type interaction and compares it to the distribution of enrichment scores expected by chance. The interaction enrichment score reflects the proportion of counts of a neighbor type relative to the product of the total number of edges detected from each constituent cell and all other cell types. For neighbor type C1Cn composed of cell types C1 ??? Cn, the enrichment score is specifically defined as:

$EnrichmentScore\_{C\_1C\_n} = \\frac{counts(C\_1C\_n)}{\\prod\_{i=1}^{i=n}totalEdges\_{C\_i}}$

The null hypothesis assumes that multiplet formation is random and thus the distribution of neighbor-types follows the underlying singlet population counts. As such, for each sample in a dataset, given n predicted singlets and m predicted multiplets consisting of x constituent cells, Neighbor-seq simulates the synthetic creation of m multiplets drawing without replacement from n+x cells. The resulting neighbor-types are tallied, and their enrichment scores are computed. This simulation is repeated for a default of 100 times, and for each neighbor type, lower tailed Wilcoxon testing quantifies the probability that the simulated enrichment scores have a central tendency greater than the observed enrichment score. All probabilities are adjusted using the Holm correction.

``` r
pred = xgpred(rf, ns.data$cell.mat)
result = multiplet_pop_sig(pred = pred, sample = meta.data$Mouse)
plot_interactions(result, 
                  legend.position = 'right', 
                  width_range = c(0.5,1))
#> Loading required package: ggraph
#> Loading required package: ggplot2
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="75%" height="75%" style="display: block; margin: auto;" />

## Neighbor-seq wrapper functions

The entire Neighbor-seq workflow can be run using two wrapper functions.

``` r
set.seed(0) # for reproducibility
# ns.data = prep_cell_mat(ge = data, celltypes = meta.data$Cluster, logfc.threshold = 0.5)
haber.ns = neighborseq(ns.data$cell.mat, ns.data$celltypes, meta.data$Mouse, do.mroc = F)

# plot with colored nodes
color = colorRampPalette((brewer.pal(9,"Greens")))(11)[6:12]; 
names(color) = c('Paneth','Stem','TA','EP','Enterocyte'); color=c(color, Goblet='gold3')

plot_interactions(result, 
                  color = color,
                  legend.position = 'right', 
                  width_range = c(0.5,1))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="75%" height="75%" style="display: block; margin: auto;" />

### Ensemble prediction

To increase reproducibility of results and minimize effects of randomness, Neighbor-seq can be run in ensemble mode; artificial multiplet construction and/or random forest training can be run for multiple iterations and an ensemble result can be computed. For each neighboring cell-type pair, the mean observed counts and Fisher???s combined adjusted p-value are reported. Multiplet type counts and p-value thresholds can be adjusted and interpreted in light of data quality and known biology.

``` r
set.seed(0) # for reproducibility
# ns.data = prep_cell_mat(ge = data, celltypes = meta.data$Cluster, logfc.threshold = 0.5)
ns = neighborseq(ns.data$cell.mat, ns.data$celltypes, meta.data$Mouse, iter = 5, do.mroc = F)

# plot with colored nodes
color = colorRampPalette((brewer.pal(9,"Greens")))(11)[6:12]; 
names(color) = c('Paneth','Stem','TA','EP','Enterocyte'); color=c(color, Goblet='gold3')

plot_interactions(ns$combined_result, 
                  combined = T,
                  legend.position = 'right', 
                  color = color, 
                  width_range = c(0.5, 1))
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="50%" height="50%" style="display: block; margin: auto;" />
