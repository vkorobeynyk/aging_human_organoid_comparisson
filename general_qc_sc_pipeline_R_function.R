#' Compute percent of counts coming from a gene set (e.g. mitochondrial, ribosomal) per cell
#' @param SO Seurat object
#' @param gene_pattern regex matched against gene symbols (e.g. "^MT-")
#' @return named numeric vector of percentages per cell (NAs, e.g. from zero-count cells, set to 0)
compute_percent_geneset = function(SO, gene_pattern) {
  genes = grep(pattern = gene_pattern, x = rownames(SO@assays$RNA$counts), value = TRUE)
  pct = (Matrix::colSums(SO@assays$RNA$counts[genes, , drop = FALSE]) / Matrix::colSums(SO@assays$RNA$counts)) * 100
  pct[is.na(pct)] = 0
  pct
}

#' Histogram of a QC metric, with optional vertical threshold lines
#' @param dt data.frame (typically colData(sce) as.data.frame())
#' @param column column name to plot
#' @param xlab x-axis label
#' @param thresholds numeric vector of x-intercepts to draw as vertical lines (NA values skipped)
plot_qc_histogram = function(dt, column, xlab, thresholds = NULL) {
  p = ggplot(dt, aes(x = .data[[column]])) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.5) +
    xlab(xlab) +
    ylab("Number of cells")
  
  for (t in thresholds) {
    if (!is.na(t)) p = p + geom_vline(xintercept = t, col = "blue")
  }
  
  p
}

#' Run the standard single-cell QC pipeline: compute mito/ribo content, filter out very
#' low-count/low-feature cells, compute per-cell QC metrics, flag outliers via adaptive
#' MADs (library size, detected features, mitochondrial %), plot before/after diagnostics,
#' and optionally filter to QC-passing cells only.
#'
#' @param SO Seurat object
#' @param to_filter if TRUE, return only cells passing QC; if FALSE, return all cells
#'   (with qc_pass metadata still computed/attached)
#' @param nmads_lower,nmads_higher MAD multiplier used for lower/higher adaptive outlier
#'   thresholds on library size and detected features (mito % only uses nmads_higher)
#' @param min_features,min_counts hard pre-filter applied before QC metric calculation -
#'   cells below these are dropped outright (doublet-detection tools generally need a
#'   minimum number of detected features to work reliably, and cells this sparse tend to
#'   be low quality regardless of where the adaptive MADs thresholds land)
#' @return filtered (or unfiltered, if to_filter = FALSE) Seurat object with QC metadata attached
run_general_qc_sc_pipeline_R = function(SO, to_filter = TRUE, nmads_lower = 2.5, nmads_higher = 2.5,
                                        min_features = 1000, min_counts = 2000)
{
  ############
  #### QC ####
  ############
  percent.mito = compute_percent_geneset(SO, "^MT-")
  SO = AddMetaData(object = SO, metadata = percent.mito, col.name = "percent.mito")
  
  percent.ribo = compute_percent_geneset(SO, "^RP[S,L]")
  SO = AddMetaData(object = SO, metadata = percent.ribo, col.name = "percent.ribo")
  
  print(VlnPlot(SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4))
  
  ##########################
  #### load data as sce ####
  ##########################
  
  # remove cells with very few features/counts: doublet-detection tools generally need a
  # minimum number of detected features to work reliably, and at very low counts/features
  # the rest of the QC (adaptive outlier detection) becomes unreliable too
  SO = SO[, SO$nFeature_RNA > min_features]
  SO = SO[, SO$nCount_RNA > min_counts]
  
  sce = as.SingleCellExperiment(SO)
  
  #########################
  #### QC b/ filtering ####
  #########################
  message("############# QC b/ filtering #############")
  #### Per cell and per feature QC
  
  # recalculate MT genes on the (now count/feature-filtered) object, since which genes
  # are present can change after the cell filtering step above
  mito.genes = grep(pattern = "^MT-", x = rownames(counts(sce)), value = TRUE)
  
  sce = addPerCellQC(sce, subsets = list(mito = mito.genes))
  
  # plot counts vs features detected
  print(plotColData(sce, x = "sum", y = "detected"))
  
  rowData(sce)$genes = rownames(sce)
  
  #### Plot explanatory Variables
  sce = logNormCounts(sce)
  
  #### Filtering parameters: detect outliers via adaptive MADs
  libsize.drop_lower = isOutlier(sce$total, nmads = nmads_lower, type = "lower", log = TRUE)
  libsize.drop_upper = isOutlier(sce$total, nmads = nmads_higher, type = "higher", log = TRUE)
  
  feature.drop_lower = isOutlier(sce$detected, nmads = nmads_lower, type = "lower", log = TRUE)
  feature.drop_upper = isOutlier(sce$detected, nmads = nmads_higher, type = "higher", log = TRUE)
  
  mito.drop = isOutlier(sce$subsets_mito_percent, nmads = nmads_higher, type = "higher")
  
  ### Quality metrics and outlier detection
  pass = !(libsize.drop_lower | libsize.drop_upper | feature.drop_lower | feature.drop_upper | mito.drop)
  sce$qc_pass = factor(ifelse(test = pass, yes = "QC_pass", no = "QC_fail"))
  
  #### Overall quality metrics before QC
  dt = colData(sce) %>% as.data.frame()
  
  print(plot_qc_histogram(
    dt, "total", "Library size",
    thresholds = c(attr(libsize.drop_lower, "thresholds")["lower"], attr(libsize.drop_upper, "thresholds")["higher"])
  ))
  
  print(plot_qc_histogram(
    dt, "detected", "Number of expressed genes",
    thresholds = c(attr(feature.drop_lower, "thresholds")["lower"], attr(feature.drop_upper, "thresholds")["higher"])
  ))
  
  print(plot_qc_histogram(
    dt, "subsets_mito_percent", "Mitochondrial proportion (%)",
    thresholds = attr(mito.drop, "thresholds")["higher"]
  ))
  
  #### QC-pass and QC-fail
  print(plotColData(sce, x = "sum", y = "detected", colour_by = "qc_pass") +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position = "right", aspect.ratio = 1))
  
  #### QC-pass only
  print(plotColData(sce[, colData(sce)$qc_pass == "QC_pass"],
                    x = "sum", y = "detected", colour_by = "qc_pass") +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position = "right", aspect.ratio = 1))
  
  #### QC-fail only
  print(plotColData(sce[, colData(sce)$qc_pass == "QC_fail"],
                    x = "sum", y = "detected", colour_by = "qc_pass") +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position = "right", aspect.ratio = 1))
  
  # Small pause: without this, in some rendering contexts the plot above can appear to
  # render after the "QC a/ filtering" message below rather than before it.
  Sys.sleep(2)
  
  #########################
  #### QC a/ filtering ####
  #########################
  message("############# QC a/ filtering #############")
  
  sce_filtered = if (to_filter) sce[, pass] else sce
  
  #### Overall quality metrics after QC
  dt = colData(sce_filtered) %>% as.data.frame()
  
  print(plot_qc_histogram(dt, "total", "Library size"))
  print(plot_qc_histogram(dt, "detected", "Number of expressed genes"))
  print(plot_qc_histogram(dt, "subsets_mito_percent", "Mitochondrial proportion (%)"))
  
  # Small pause: same rendering-order reason as above.
  Sys.sleep(2)
  
  ###############################
  #### Additional QC metrics ####
  ###############################
  message("############# Additional QC metrics #############")
  
  plot(density(sce_filtered$total / 1e3),
       xlab = "library size (thousands)",
       main = "")
  rug(sce_filtered$total / 1e3)
  
  plot(y = sce_filtered$total / 1e3,
       x = sce_filtered$subsets_mito_percent,
       pch = 20,
       ylab = "library size (thousands)",
       xlab = "mitochondrial proportion (%)")
  
  plot(y = sce_filtered$total / 1e3,
       x = sce_filtered$detected,
       pch = 20,
       ylab = "library size (thousands)",
       xlab = "number of genes")
  
  return(as.Seurat(sce_filtered))
}