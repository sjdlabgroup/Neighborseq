#' Prepare Neighbor-seq inputs: gene expression matrix, clusters, and cluster marker genes
#'
#' @param ge Gene by barcode single-cell gene expression matrix
#' @param celltypes Predetermined cell cluster annotations - optional
#' @param logfc.threshold Minimum average log fold change of a gene in a cell-type compared to all other cell-types to be tested for differential expression. Default is 1; decreaseing the threshold will find weaker signals but will increase computation time. See documentation for Seurat FindMarkers function for more details.
#' @param min.pct Minimum percent of cells from a cluster to express a gene for that gene to be tested for differential expression. Default is 0.25. See documentation for Seurat FindMarkers function for more details.
#' @param max.cells Maximum number of cells to randomly select from each cell cluster during differential gene expression testing. Default is 200.
#' @param topn Number of marker genes to select for each cell-type. The union of these - there may be overlap-  will be the genes returned in the output gene expression matrix. Top arker genes are selected based on their fold-change rank.
#' @param res Resolution parameter for Seurat clustering. Larger values will find a larger number of clusters. See Seurat FindClusters function for more details.
#' @param ... Arguments passed to other methods
#'
#' @return A list containing the cell matrix, cell-types, and marker gene output
#' @export
#'
#' @import Seurat
#' @import tidyverse
#' @import xgboost
#' @import multiROC
#' @import parallel
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import ggraph
#' @import metap
#' @import RColorBrewer
#' @import stringr
#' @import ggpubr
#' @import ggplot2
#' @importFrom tidyr "separate"
#' @importFrom dplyr "sample_frac"
#' @importFrom dplyr "group_by"
#' @importFrom dplyr "top_n"
#' @importFrom dplyr "anti_join"
#' @importFrom dplyr "bind_rows"
#' @importFrom dplyr "summarize"
#' @importFrom gtools "combinations"
#' @importFrom tibble "rownames_to_column"
#' @importFrom tibble "column_to_rownames"
#'
prep_cell_mat = function(ge, celltypes = NULL, logfc.threshold = 1, min.pct = 0.25, max.cells = 200, topn = 50, res = 0.8, ...){

  ge = ge %>%
    CreateSeuratObject() %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 5000)

  if(is.null(celltypes)){
    ge = ScaleData(ge) %>% RunPCA() %>% FindNeighbors() %>% FindClusters()
    while(nlevels(ge$seurat_clusters) > 12){
      res = 0.5*res; ge = FindClusters(ge, resolution = res)
    }
    celltypes = ge$seurat_clusters
  }

  Idents(ge) = celltypes
  markers = FindAllMarkers(ge, logfc.threshold = logfc.threshold,
                           only.pos = T,
                           min.pct = min.pct,
                           max.cells.per.ident = max.cells, ...)

  wt = ifelse('avg_log2FC' %in% colnames(markers), 'avg_log2FC', 'avg_logFC')
  m = markers %>% group_by(cluster) %>% top_n(n = topn, wt = wt)
  m = unique(m$gene)

  list(cell.mat = ge@assays$RNA@counts[m, ], celltypes = celltypes, markers = markers)
}

#' Enumerate all possible cell-type combinations of n cells
#'
#' @param celltypes A vector of cell-type annotations. There must be one annotation for each barcode.
#' @param n Multiplet degree, i.e. the maximum number of cells in each multiplet. Combinations will be included for integers from 2 to n.
#' @param exclude Optional cell-type(s) to exclude when enumerating multiplet-types.
#'
#' @return A dataframe of n columns with each row containing a multiplet type. Empty columsn are filled with NA.
#' @export
#'

multiplet_types = function(celltypes, n = 2, exclude = NULL){
  celltypes = celltypes %>% as.character() %>% unique()

  if(!is.null(exclude)){celltypes = setdiff(celltypes, exclude)}
  d = list()
  for(i in 2:n){
    d[[i]] = combinations(length(celltypes), i, celltypes, repeats.allowed = T) %>% data.frame()
  }
  rbindlist(d, fill = T) %>% data.frame()
}

#' Construct a dataset of artificial multiplets
#'
#' @param cell.mat Marker gene by barcode matrix
#' @param celltypes A vector of cell-type annotations
#' @param multiplet_classes A dataframe of multiplet types
#' @param n Number of artificial multiplets to construct for each multiplet type. Will also include n singlets from each cell-type.
#'
#' @return A dataframe where each row is a singlet or multiplet and the columns contain genes counts. The last column called 'type' contains the barcode type.
#' @export
#'
artificial_multiplets = function(cell.mat, celltypes, multiplet_classes, n = 100){
  celltypes = as.character(celltypes)

  i = which(celltypes %in% (c(multiplet_classes) %>% unlist() %>% as.character() %>% unique()))
  celltypes = celltypes[i]
  cell.mat = cell.mat[, i]

  multiplet_type = apply(multiplet_classes, 1, function(x) paste(x[!is.na(x)], collapse = '_'))
  idx = sapply(unique(celltypes) %>% sort(), function(x) which(celltypes == x), simplify = F)

  am = list()
  counter = 0
  for(i in 1:nrow(multiplet_classes)){
    ii = which(!is.na(multiplet_classes[i,]))
    for(j in 1:n){
      # print(paste(i,j))
      iii = sapply(ii, function(x) idx[[multiplet_classes[i,x]]] %>% sample(1))
      xx = data.frame(cell.mat[, iii] %>% rowSums() %>% t(), type = multiplet_type[i], check.names = F)
      counter = counter + 1
      am[[counter]] = xx
    }
  }
  am = rbindlist(am)

  if(any(table(celltypes) < n)){
    print('warning: artificial multiplets constructed by sampling with replacement'); replace = T
  } else {replace = F}

  xx = lapply(idx, function(x) cell.mat[, sample(x, n, replace = replace)] %>% t() %>% data.frame(check.names=F)) %>% rbindlist()
  xx = cbind(xx, type = rep(names(idx), each = n))

  colnames(am) = str_replace_all(colnames(am), '-', '.')
  colnames(xx) = str_replace_all(colnames(xx), '-', '.')
  am = rbind(am, xx)
  am$type = factor(am$type)
  am
}

#' Train and test a random forest classifier to predict barcode composition
#'
#' @param artificial_multiplets_mat Dataframe containing artificial multiplets. Each row is a singlet or multiplet and the columns are gene counts. The last column contains the barcode type.
#' @param f Fraction of multiplets to use for training the random forest
#' @param ... Arguments passed to other methods
#'
#' @return A list containing the trained random forest, barcode-type prediction probabilities for each barcode in the test dataset, the training dataset, and the test dataset
#' @export
#'
multiplet_rf = function(artificial_multiplets_mat, f = 0.8, ...){
  artificial_multiplets_mat$type = factor(artificial_multiplets_mat$type)
  artificial_multiplets_mat = artificial_multiplets_mat %>% rownames_to_column('n')
  train = artificial_multiplets_mat %>% group_by(type) %>% sample_frac(f, replace = F)
  test = anti_join(artificial_multiplets_mat, train, by = 'n')
  artificial_multiplets_mat = column_to_rownames(artificial_multiplets_mat, 'n')
  train = column_to_rownames(train, 'n')
  test = column_to_rownames(test, 'n')

  train.type = train[, 'type']
  train.typex = train[, 'type'] %>% as.numeric(); train.typex = train.typex - 1
  train = train[,-ncol(train)] %>% as.matrix()

  ncore = detectCores()

  rf <- xgboost(data = train, label = train.typex, objective = "multi:softprob",
                eval_metric = "mlogloss",
                num_class = nlevels(artificial_multiplets_mat$type),
                nthread = ncore, nround = 1, max_depth = 20,
                num_parallel_tree = 200, subsample = 0.632,
                colsample_bytree  = sqrt(ncol(train))/ncol(train),
                colsample_bynode  = sqrt(ncol(train))/ncol(train))

  test.type = test[, 'type']
  test = test[, -ncol(test)] %>% as.matrix()

  pred = predict(rf, test)
  pred <- matrix(pred,
                 nrow = nlevels(artificial_multiplets_mat$type),
                 ncol = length(pred)/nlevels(artificial_multiplets_mat$type)) %>%
    t() %>%
    data.frame(check.names = F)

  colnames(pred) = levels(artificial_multiplets_mat$type)

  list(rf = rf, pred = pred, train = data.frame(train, type = train.type, check.names = F), test = data.frame(test, type = test.type, check.names = F))
}

#' Calculate multivariate receiver operator sensitivities and specificities
#'
#' @param test Test dataset used to evalute the multiplet random forest. The last column must be labeled 'type' and contain the barcode type.
#' @param pred Prediction probabilities for each barcode.
#'
#' @return A dataframe containing the sensitivity, specificity, and AUC for each group. Also contains the micro and macro avareages of all groups. See multiROC package documentaiton for more details.
#' @export
#'
mroc_format = function(test, pred){
  test$type = factor(test$type)
  classes = levels(test$type)
  mat = matrix(0, nrow = nrow(pred), ncol = length(classes))
  colnames(mat) = classes
  for(i in 1:nrow(mat)){
    mat[i, test$type[i]] = 1
  }
  colnames(mat) = colnames(mat) %>% str_replace_all(' ', '_') %>% paste0('_true')
  pred.mat = pred
  colnames(pred.mat) = colnames(pred.mat) %>% str_replace_all(' ', '_') %>% paste0('_pred_RF')

  df = data.frame(mat, pred.mat, check.names = F)
  mroc = suppressWarnings(multi_roc(df, force_diag = T))

  # format data.frame for plotting
  n_method <- length(unique(mroc$Methods))
  n_group <- length(unique(mroc$Groups))
  mroc_df <- data.frame(Specificity= numeric(0), Sensitivity= numeric(0),
                        Group = character(0), AUC = numeric(0),
                        Method = character(0))
  for (i in 1:n_method) {
    for (j in 1:n_group) {
      temp_data_1 <- data.frame(Specificity=mroc$Specificity[[i]][j],
                                Sensitivity=mroc$Sensitivity[[i]][j],
                                Group=unique(mroc$Groups)[j],
                                AUC=mroc$AUC[[i]][j],
                                Method = unique(mroc$Methods)[i])
      colnames(temp_data_1) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
      mroc_df <- rbind(mroc_df, temp_data_1)

    }
    temp_data_2 <- data.frame(Specificity=mroc$Specificity[[i]][n_group+1],
                              Sensitivity=mroc$Sensitivity[[i]][n_group+1],
                              Group= "Macro",
                              AUC=mroc$AUC[[i]][n_group+1],
                              Method = unique(mroc$Methods)[i])
    temp_data_3 <- data.frame(Specificity=mroc$Specificity[[i]][n_group+2],
                              Sensitivity=mroc$Sensitivity[[i]][n_group+2],
                              Group= "Micro",
                              AUC=mroc$AUC[[i]][n_group+2],
                              Method = unique(mroc$Methods)[i])
    colnames(temp_data_2) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
    colnames(temp_data_3) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
    mroc_df <- rbind(mroc_df, temp_data_2)
    mroc_df <- rbind(mroc_df, temp_data_3)
  }

  mroc_df
}

#' Plot multivariate receiver operator curves
#'
#' @param mroc_df multiROC result dataframe obtained from mroc_format function
#' @param c_size linewidth for individual barcode type curves
#' @param a_size linewidth for average ROC curve
#'
#' @return A ggplot object graphing multiROC curves
#' @export
#'
mroc_plot = function(mroc_df, c_size = 0.3, a_size = 1){
  macro.auc = mroc_df$AUC[mroc_df$Group == 'Macro'][1] %>% round(3)
  ggplot(mroc_df, aes(x = 1-Specificity, y=Sensitivity)) +
    geom_line(data = subset(mroc_df, Group != 'Macro'), aes(color = Group), size = c_size, show.legend = F) +
    scale_color_manual(values = colorRampPalette(brewer.pal(9, 'Greys'))(length(unique(mroc_df$Group)))) +
    geom_line(data = subset(mroc_df, Group == 'Macro'), color = 'Firebrick', size = a_size, show.legend = T) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    theme_classic() +
    scale_x_continuous(limits = c(-0.05, 1.05), expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0,0)) +
    ggtitle(paste('Artificial droplet classification\nAvg AUC =', macro.auc)) +
    theme(axis.text = element_text(colour = 'black', size = 10),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          title = element_text(size=9),
          legend.position = 'none')
}


#' Format dataframe to create a clustered dot-plot
#'
#' @param df A dataframe in long format containing data for the dotplot
#' @param rows Name of dataframe column containing row labels for dotplot
#' @param columns Name of dataframe column containing column labels for dotplot
#' @param values Name of dataframe column containing values for dotplot
#'
#' @return A long format dataframe with factored and leveled row and column labels
#' @export
#'
cluster_dotplot = function(df, rows, columns, values){
  xx = df %>% pivot_wider(id_cols = rows,
                          names_from = columns,
                          values_from = values) %>%
    column_to_rownames(rows)
  xx[is.na(xx)] = 0
  row.order = hclust(dist(xx))$order
  column.order = hclust(dist(t(xx)))$order
  df = df %>% mutate(!!columns := factor(!!as.name(columns), levels = colnames(xx)[column.order]))
  df = df %>% mutate(!!rows := factor(!!as.name(rows), levels = rownames(xx)[row.order]))
  return(df)
}

#' Barcode-type prediction using a trained random forest
#'
#' @param rf Random forest trained on artificial multiplets. Output of multiplet_rf
#' @param cell.mat Gene by cell matrix
#'
#' @return Dataframe of barcode-type prediction probabilities
#' @export
#'
xgpred = function(rf, cell.mat){
  # format gene expresion matrix and gene names for the random forest
  pred.mat = t(cell.mat)
  colnames(pred.mat) = str_replace_all(colnames(pred.mat), '-', '.')
  # colnames(pred.mat)[str_which(colnames(pred.mat), '^[0-9]')] =
  #   paste0('X', colnames(pred.mat)[str_which(colnames(pred.mat), '^[0-9]')])

  pred = predict(rf$rf, pred.mat, type = 'prob')
  pred <- matrix(pred,
                 nrow = nlevels(rf$train$type),
                 ncol = length(pred)/nlevels(rf$train$type)) %>%
    t() %>%
    data.frame()

  colnames(pred) = levels(rf$train$type)
  pred
}

#' Compute multiplet type enrichment scores and statistical significance
#'
#' @param pred Multiplet type prediction probabilities for each barcode
#' @param sample Optional sample label annotation for each barcode. Statistical significance is assessed relative to cell-type proportions from an individual sample only.
#' @param n Number of simulations to run. In ach simulation multiplets are randomly created from the underlying cell-type distributins.
#' @param homotypic True or False; whether to include homotypic multiplets when calculating enrichment and statistical significance
#'
#' @return A dataframe containing the multiplet significance results. Each row contains data for a multiplet type in a sample. The sample, multiplet type, enrichment score, p-value, and adjusted p-value are reported
#' @export
#'
multiplet_pop_sig = function(pred, sample = NULL, n = 100, homotypic = T){

  if(is.null(sample)){sample = rep('Sample.1', nrow(pred))}

  out = list()
  for(i in unique(sample)){

    pred.sample = pred[sample == i, ] %>% apply(1, function(x) names(x)[which.max(x)])
    celltypes = pred.sample %>% strsplit('_') %>% unlist()
    pred.sample = pred.sample[str_which(pred.sample, '_')]

    if(length(pred.sample) == 0){next}

    df = table(pred.sample) %>% data.frame() %>% separate('pred.sample', into=c('Cell_1','Cell_2'))
    colnames(df)[3] = 'Counts'

    if(homotypic == F){
      j = apply(df, 1, function(x) length(unique(x[1:(length(x)-1)])))
      df = df[j>1, ]
      rownames(df) = NULL
      pred.sample = pred.sample[j>1]
    }

    y = sapply(1:n, function(x){
      sample(celltypes, max(4, 2*length(pred.sample)), replace = F) %>%
        matrix(ncol = 2) %>%
        apply(1, sort) %>%
        t() %>%
        apply(1, paste0, collapse = '_') %>%
        table() %>%
        as.matrix() %>%
        t() %>%
        data.frame(check.names = F)
    }) %>% bind_rows(); y[is.na(y)] = 0


    types = paste0(df$Cell_1, '_', df$Cell_2)
    total_edges = c()
    for(j in unique(celltypes)){
      total_edges[j] = sum(df$Counts[str_which(types, j)])
    }

    e2 = c()
    for(j in 1:nrow(df)){
      e2[j] = df$Counts[j]/(total_edges[df$Cell_1[j]] * total_edges[df$Cell_2[j]])
    }

    names(e2) = types; e2=e2*10^4

    df$EnrichmentScore = e2

    # enrichment score for randomized data
    cn = colnames(y) %>% data.frame() %>% separate(col='.', into = c('cell1', 'cell2'))
    ct = unique(c(cn$cell1, cn$cell2))
    yedges = matrix(0, nrow = nrow(y), ncol = length(ct))
    colnames(yedges) = ct
    for(j in 1:nrow(y)){
      for(k in ct) {
        yedges[j,k] = y[j, str_which(colnames(y), k)] %>% sum()
      }
    }

    for(j in 1:nrow(y)){
      for(k in 1:ncol(y)){
        y[j,k] = y[j,k]/(yedges[j, cn$cell1[k]] * yedges[j, cn$cell2[k]])
      }
    }
    y = y*10^4


    sdiff = setdiff(names(e2), colnames(y))
    if(length(sdiff) > 0){y = cbind(y, matrix(0, nrow = nrow(y), ncol = length(sdiff),
                                              dimnames = list(rownames=NULL, colnames=sdiff)))}

    pval = sapply(names(e2), function(x) wilcox.test(y[, x], mu = e2[x], alternative = 'less')$p.value)
    padj = p.adjust(pval)


    df = data.frame(sample = i, df, pval = pval, padj); rownames(df) = NULL
    out[[i]] = df
  }
  out = rbindlist(out) %>% data.frame()

  out
}

#' Function to run the entire Neighbor-seq pipeline
#'
#' @param cell.mat Gene by cell matrix. Performance is improved if only marker genes are kept.
#' @param celltypes Vector containing cell-type annotations for each barcode
#' @param sample Optional vector containing sample annotations for each barcode
#' @param iter Optional number of iterations to train the random forest
#' @param exclude Optional cell-type to exclude from the artificial multiplet dataset. Can be uesful when known multiplets are present.
#' @param multiplet.degree Maximum number of cells in a multiplet
#' @param n.am Number of artificial multiplets to create for each multiplet type
#' @param f Fraction of artificial multiplets data to use for training the random forest
#' @param nsim Number of simulations to run to calculate multiplet enrichment significance
#' @param do.mroc True or False; whether to calculate multivariate receiver operator curves
#' @param homotypic True or False; whether to include homotypic multiplets when assessing enrichment and significance
#'
#' @return A list containing: 1. A combined result dataframe with multiplet enrichment scores and Fisher combined p-values if iter>1, 2. a list or dataframe with multiplet enrichment scores and p-values from each iteration, 3. a list or dataframe containing barcode type prediction probabilities for all barcodes in the dataset, 4. a list or dataframe containing the multiROC analysis
#' @export
#'
neighborseq = function(cell.mat, celltypes, sample = NULL, iter = 1, exclude = NULL, multiplet.degree = 2, n.am = 100, f = 0.8, nsim = 100, do.mroc = T, homotypic = T){
  start = Sys.time()
  mt = multiplet_types(celltypes, n = multiplet.degree, exclude = exclude)

  print('Creating artificial multiplets')
  am = artificial_multiplets(cell.mat, celltypes, mt, n = n.am)

  rf = list()
  mroc = list()
  pred = list()
  result = list()
  for(i in 1:iter){
    print(paste('Training random forest', i, '/', iter))
    rf[[i]] = multiplet_rf(am, f = f)

    if(do.mroc == T){
      print('Computing multi-ROC')
      mroc[[i]] = mroc_format(rf[[i]]$test, rf[[i]]$pred)
    }
    else{mroc = NULL}

    pred[[i]] = xgpred(rf[[i]], cell.mat)

    print('Computing population significance')
    result[[i]] = multiplet_pop_sig(pred[[i]], sample = sample, n = nsim, homotypic = homotypic)

  }

  if(iter > 1){
    combined_result = rbindlist(result) %>%
      group_by(sample, Cell_1, Cell_2) %>%
      summarize(MeanCounts = mean(Counts),
                CombinedP = suppressWarnings(sumlog(padj)$p))
  } else {
    result = result[[1]]
    pred = pred[[1]]
    rf = rf[[1]]
    combined_result = NULL
    if(!is.null(mroc)){mroc = mroc[[1]]}
  }

  print(round(Sys.time() - start), 1)

  list(combined_result = combined_result, result = result, pred = pred, rf = rf, mroc = mroc)
}

#' Plot the cell-cell interaction network
#'
#' @param df Dataframe from containing the multiplet enrichment and significance results
#' @param minCounts_sample Minimum number of multiplet type counts in a sample to include the plot
#' @param minCounts_total Minimum number of multiplet type counts in all samples to include the plot
#' @param maxP Maximum multiplet type p-value to include in the plot
#' @param padding Line padding around cell-type labels
#' @param combined True or False; whether the dataframe contains an ensemble Neighbor-seq result or a single run
#' @param layout Network graph layout; see ggraph documentation for more information
#' @param color A named vector of color lables for each cell-type
#' @param color_breaks Edge color breaks. See ggplot documentation for more details
#' @param width_breaks Edge width breaks. See ggplot documentation for more details
#' @param width_range Edge width range. See ggplot documentation for more details
#' @param min_color Edge color for smallest edge / multiplet type with fewest counts
#' @param max_color Edge color for largest edge / multiplet type with highest counts
#' @param legend.position Legend position; either 'top', 'bottom', 'left', 'right', 'none'
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object containing the cell-cell interaction network
#' @export
#'
plot_interactions = function(df, minCounts_sample = 10, minCounts_total = 10, maxP = 0.05, padding = 1, combined = F,
                             layout = 'stress', color = NULL, color_breaks = NULL, width_breaks = NULL,
                             width_range = c(0.1,1), min_color = 'grey50', max_color = 'black',
                             legend.position = 'none', ...){

  require(ggraph)

  df$type = paste0(df$Cell_1, '_', df$Cell_2)

  n = padding
  if(combined == T){
    df = subset(df, MeanCounts > minCounts_sample & CombinedP < maxP)
    df = df %>%
      group_by(type, Cell_1, Cell_2) %>%
      summarize(MeanCounts = sum(MeanCounts, na.rm=T),
                CombinedP = min(CombinedP, na.rm=T),
                .groups = 'keep')
    edge_width = -1*log10(df$CombinedP)
    edge_color = df$MeanCounts
  } else {
    df = subset(df, Counts > minCounts_sample & padj < maxP)
    df = df %>%
      group_by(type, Cell_1, Cell_2) %>%
      summarize(Counts = sum(Counts, na.rm=T),
                padj = min(padj),
                .groups = 'keep')
    edge_width = -1*log10(df$padj)
    edge_color = df$Counts
  }
  edge_width = edge_width[df$Cell_1 != df$Cell_2]
  edge_color = edge_color[df$Cell_1 != df$Cell_2]
  df = df[df$Cell_1 != df$Cell_2, ]

  graph = graph_from_data_frame(df[,2:4])

  if(!is.null(color)){
    color = color[V(graph) %>% as.factor() %>% names()]
  } else {
    color = rep('black', length(V(graph)))
  }

  if(is.null(color_breaks)){color_breaks = waiver()}
  if(is.null(width_breaks)){width_breaks = waiver()}

  ggraph(graph, layout = layout) +
    geom_edge_link(aes(start_cap = label_rect(node1.name, padding = ggplot2::margin(n,n,n,n,'mm')),
                       end_cap = label_rect(node2.name, padding = ggplot2::margin(n,n,n,n,'mm')),
                       width = edge_width,
                       color = edge_color)) +
    scale_edge_width(range = width_range, name = '-log10(P)', breaks = width_breaks) +
    scale_edge_color_gradient(low = min_color, high = max_color, name = 'Counts', breaks = color_breaks) +
    geom_node_text(aes(label = name), color = color, size = 2.5, repel = T, point.padding = NA, box.padding = 0, force = 0.001)+
    theme(panel.background = element_rect(fill = 'white'),
          legend.position = legend.position,
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size=7),
          legend.key=element_blank())
}


