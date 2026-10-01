# -----------------------------------------------------------------------------
# adapted from DoubletFinder - Plotting disabled.
# -----------------------------------------------------------------------------
find.pK_noplot = function (sweep.stats) 
{
  "%ni%" = Negate("%in%")
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    bc.mvn = as.data.frame(matrix(0L, nrow = length(unique(sweep.stats$pK)), 
                                  ncol = 5))
    colnames(bc.mvn) = c("ParamID", "pK", "MeanBC", "VarBC", 
                         "BCmetric")
    bc.mvn$pK = unique(sweep.stats$pK)
    bc.mvn$ParamID = 1:nrow(bc.mvn)
    x = 0
    for (i in unique(bc.mvn$pK)) {
      x = x + 1
      ind = which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] = mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] = sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] = mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, 
                                                                            "BCreal"])^2)
    }
    return(bc.mvn)
  }
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    bc.mvn = as.data.frame(matrix(0L, nrow = length(unique(sweep.stats$pK)), 
                                  ncol = 6))
    colnames(bc.mvn) = c("ParamID", "pK", "MeanAUC", "MeanBC", 
                         "VarBC", "BCmetric")
    bc.mvn$pK = unique(sweep.stats$pK)
    bc.mvn$ParamID = 1:nrow(bc.mvn)
    x = 0
    for (i in unique(bc.mvn$pK)) {
      x = x + 1
      ind = which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] = mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] = mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] = sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] = mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, 
                                                                            "BCreal"])^2)
    }
    
    return(bc.mvn)
  }
}

#' Run GSEA on a ranked gene list (ranked by avg_log2FC) and optionally plot it
#' @param plot_type one of "classic", "volcano", "ridge", "bar", or anything else for no plot
#' @return the full GSEA result list, or 0 if nothing was significantly enriched
#'   (both this file's callers and the analysis scripts check `length(result) == 1`
#'   to detect this "nothing enriched" case)
gsea_analysis = function(markers,cluster ,gs  = NULL, n_pathways = NULL,plot_type)
{
  if(is.null(gs)) {stop("gsea_analysis(): please specify a geneset (gs)")}
  if(is.null(n_pathways)) {n_pathways = 100}
  # order markers based on column avg_log2FC
  markers %<>% dplyr::arrange(,desc(avg_log2FC))
  
  rank = markers$avg_log2FC 
  names(rank) = markers$gene
  
  rank = rank[!duplicated(names(rank))] # remove few duplicated gene IDs
  
  # run GSEA
  safe_gsea = purrr::possibly(genGSEA_copy, otherwise = NULL)
  gse = safe_gsea(genelist = rank,geneset = gs, set_seed = 1)
  if(is.null(gse)) {return(0)} # nothing significantly enriched
  gse_original = gse
  pathways = gse$gsea_df %>% dplyr::arrange(,desc(abs(NES))) %>% dplyr::select(1) %>% unlist %>% magrittr::extract(1:n_pathways)
  pathways = pathways[!is.na(pathways)]
  if(plot_type == "classic") {plotGSEA(gse, plot_type = "classic", show_pathway = pathways, label_by = "description") # , show_gene = genes
  } else if(plot_type == "volcano") {plotGSEA(gse, plot_type = "volcano", show_pathway = pathways, label_by = "description")
  } else if(plot_type == "ridge") {plotGSEA(gse, plot_type = "ridge", show_pathway = pathways, label_by = "description")
  } else if(plot_type == "bar") {
    # while plot == "bar" it always displays all pathways, which become very crowded if we have many entries
    gse$gsea_df = gse$gsea_df %>% dplyr::arrange(,desc(abs(NES))) %>% dplyr::slice(1:n_pathways)
    plotGSEA(gse, plot_type = "bar", colour = c("turquoise3", "orange"), label_by = "description")}
  
  return(gse_result = gse_original)
}

#' Fit a single cv.glmnet model for the given family, with the lasso-like settings
#' @param train_data data.frame/matrix of predictors (no label column)
#' @param train_label response vector
#' @param model_family "gaussian" or "binomial"
fit_glmnet = function(train_data, train_label, model_family) {
  if (model_family == "gaussian") {
    cv.glmnet(as.matrix(train_data), train_label,
              family = "gaussian", type.measure = "mse", keep = TRUE)
  } else if (model_family == "binomial") {
    cv.glmnet(as.matrix(train_data), factor(train_label, levels = c("ICE", "wt")), # ICE as baseline
              family = "binomial", type.measure = "auc", keep = TRUE)
  } else {
    stop("fit_glmnet(): model_family must be 'gaussian' or 'binomial', got: ", model_family)
  }
}

#' Build the per-seed result (model params, coefficients, predictions) for one
#' trained glmnet model, in the family-specific shape train_glmnet_model()'s callers expect
build_seed_result = function(cvfit, test_label, preds, seed, model_family, test_metadata = NULL) {
  index_lambda = cvfit$index["min", ]
  result = list(model_params = list(lambda.min = cvfit$lambda.min, seed = seed), coeffs = list())
  
  if (model_family == "gaussian") {
    result$model_params$MSE_cv = cvfit$cvm[index_lambda]
    result$model_params$MSE_cv_sd = cvfit$cvsd[index_lambda]
    result$model_params$MSE_test = ModelMetrics::mse(test_label, preds[, 1])
    result$coeffs$gene_imp = coef(cvfit) %>% as.matrix() %>% as.data.frame() %>% round(3)
    result$predictions = data.frame(
      pred = preds[, 1],
      chronological_age = test_label,
      seed = seed,
      sex = test_metadata$Sex,
      diagnosis = test_metadata$Diagnosis,
      sample_source = test_metadata$Sample_Source
    )
  } else if (model_family == "binomial") {
    result$model_params$AUC_cv = cvfit$cvm[index_lambda]
    result$model_params$AUC_cv_sd = cvfit$cvsd[index_lambda]
    result$coeffs$gene_imp = coef(cvfit) %>% as.data.frame()
    result$predictions = data.frame(pred = preds[, 1], phenotype = test_label, seed = seed)
  }
  
  result
}

#' Train a glmnet model predicting `df$label` from the remaining columns of `df`
#' @param df data.frame with predictor columns plus a 'label' column
#' @param model_family "gaussian" (continuous label, e.g. Age) or "binomial" (e.g. phenotype)
#' @param n_times number of random 80/20 train/test splits to repeat when test = TRUE
#' @param metadata full metadata (same row order as df) - Sex/Diagnosis/Sample_Source are
#'   attached to the gaussian predictions for downstream covariate plots
#' @param test if TRUE, repeat n_times train/test splits and return per-seed results in
#'   $lst; if FALSE, fit a single model on all of df (no held-out test set) and return an
#'   empty $lst
#' @return list(lst = per-seed results (test=TRUE only), cvfit = final/last fitted cv.glmnet)
train_glmnet_model = function(df, model_family, n_times, metadata, test = TRUE)
{
  if (test) {
    lst = list()
    for (seed in seq_len(n_times)) {
      message(seed)
      set.seed(seed)
      
      train_index = createDataPartition(df$label, p = 0.8, list = FALSE)
      train_data = df[train_index, ]
      test_data  = df[-train_index, ]
      test_metadata = metadata[-train_index, ]
      
      train_label = train_data$label
      train_data$label = NULL
      cvfit = fit_glmnet(train_data, train_label, model_family)
      
      ################
      ### predict  ###
      ################
      test_label = test_data$label
      test_data$label = NULL
      preds = predict(cvfit, newx = as.matrix(test_data), s = "lambda.min", type = "response")
      
      lst[[paste0("seed_", seed)]] = build_seed_result(cvfit, test_label, preds, seed, model_family, test_metadata)
    }
  } else {
    lst = list()
    train_label = df$label
    df$label = NULL
    cvfit = fit_glmnet(df, train_label, model_family)
  }
  
  return(list(lst = lst, cvfit = cvfit))
}

#' Metacell-level DESeq2 DEA + PCA/correlation-heatmap diagnostics, used throughout
#' the human and organoid metacell analysis scripts.
#' @param organism "human" or "organoid" - controls which metadata columns are used
#'   for the results contrast and the heatmap annotation
metacell_DESeq2_DEA = function(sce, phenotype_column, dataset_column, design_formula, title_corHeatmap = FALSE, contrast, organism)
{
  dds = DESeqDataSetFromMatrix(counts(sce), 
                               colData = colData(sce), 
                               design = design_formula)
  dds = estimateSizeFactors(dds)
  
  
  vst = DESeq2::vst(dds, blind=TRUE)
  
  # Plot PCA
  plt_phenotype = DESeq2::plotPCA(vst, intgroup = phenotype_column)
  plt_dataset = DESeq2::plotPCA(vst, intgroup = dataset_column)
  
  # Run DESeq2 differential expression analysis
  dds = DESeq(dds)
  
  # Extract the vst matrix from the object and compute pairwise correlation values
  vst_mat = assay(vst)
  vst_cor = cor(vst_mat)
  
  md = colData(sce) %>% as.data.frame
  if(organism == "human")
  {
    res = results(dds ,name = contrast) %>% 
      as.data.frame %>%
      dplyr::rename("avg_log2FC" = "log2FoldChange") %>%
      mutate(gene = rownames(.)) 
    
    # for heatmap
    df = data.frame(phenotype = md %>% pull(phenotype_column), 
                    diagnosis     = md %>% pull(Diagnosis),
                    Sample_source = md %>% pull(Sample_Source),
                    row.names     = colnames(vst_mat))
    colnames(df)[1] = phenotype_column # change to "Age"
    
    diagnosis_levels = sort(unique(as.character(df$diagnosis)))
    sample_source_levels = sort(unique(as.character(df$Sample_source)))
    
    annotation_colors = list(
      diagnosis = setNames(RColorBrewer::brewer.pal(max(3, length(diagnosis_levels)), "Set1")[seq_along(diagnosis_levels)], diagnosis_levels),
      Sample_source = setNames(scales::hue_pal()(length(sample_source_levels)), sample_source_levels)
    )
    annotation_colors[[phenotype_column]] = c("white", "#08519C") # continuous covariate: white (low) -> dark blue (high)
    
  } else if(organism == "organoid") {
    res = results(dds ,contrast = contrast) %>% 
      as.data.frame %>%
      dplyr::rename("avg_log2FC" = "log2FoldChange") %>%
      mutate(gene = rownames(.)) 
    
    # for heatmap
    df = data.frame(phenotype = md %>% pull(phenotype_column) , 
                    maturation_stage = md %>% pull(maturation_stage), 
                    dataset = md %>% pull(dataset), 
                    row.names = colnames(vst_mat))
    
    phenotype_levels = sort(unique(as.character(df$phenotype)))
    maturation_levels = sort(unique(as.character(df$maturation_stage)))
    dataset_levels = sort(unique(as.character(df$dataset)))
    
    annotation_colors = list(
      phenotype = setNames(c("#D85A30", "#1D9E75")[seq_along(phenotype_levels)], phenotype_levels),
      maturation_stage = setNames(scales::hue_pal()(length(maturation_levels)), maturation_levels),
      dataset = setNames(scales::hue_pal()(length(dataset_levels)), dataset_levels)
    )
  }
  
  # generate heatmap
  heatmap = pheatmap(vst_cor, annotation_col = df, 
                     annotation_colors = annotation_colors,
                     show_colnames = FALSE, 
                     show_rownames = FALSE,
                     main = if (isFALSE(title_corHeatmap)) NA else title_corHeatmap,
                     silent = TRUE)
  
  return(list(res = results(dds),
              markers = res,
              vst_transformed_deseq2 = vst,
              dds = dds,
              plt_phenotype = plt_phenotype,
              plt_dataset = plt_dataset,
              heatmap = heatmap))
}


reduceGOredundancy = function(gsea_output, threshold = 0.7)
{
  simMatrix = calculateSimMatrix(gsea_output$hsapiens_BP_ID,
                                 orgdb="org.Hs.eg.db",
                                 ont="BP",
                                 method="Rel")
  
  # grouping of terms based on similarity
  scores = setNames(-log10(gsea_output$p.adjust), gsea_output$hsapiens_BP_ID)
  reducedTerms = reduceSimMatrix(simMatrix,
                                 scores,
                                 threshold=threshold,
                                 orgdb="org.Hs.eg.db")
  
  return(reducedTerms)
}

#' Build a treemap-ready reduced-GO-terms table for one NES direction of a GSEA
#' result, and plot the treemap as a side effect (except in the single-term case,
#' where the original clustering can't run on just one term).
#' @return data.frame ready for treemap indexing, or NULL if nothing enriched in this direction
summarize_and_plot_GOterms_direction = function(gsea_df_subset, title_label, direction_label, ont)
{
  if (nrow(gsea_df_subset) == 0) return(NULL)
  
  if (nrow(gsea_df_subset) == 1) {
    # Can't cluster a single term - just reshape it into the same structure
    # reduceGOredundancy() would normally produce, so callers get a consistent shape.
    return(data.frame(
      go = gsea_df_subset[, 1],
      cluster = 1,
      parent = gsea_df_subset[, 1],
      score = 1,
      size = 1,
      term = gsea_df_subset$Description,
      parentTerm = gsea_df_subset$Description,
      termUniqueness = 1,
      termUniquenessWithinCluster = 1,
      termDispensability = 1
    ))
  }
  
  if (ont != "BP") return(NULL)
  
  reduced = gsea_df_subset %>% dplyr::select(c("hsapiens_BP_ID", "p.adjust")) %>% reduceGOredundancy()
  
  reduced %>%
    treemap::treemap(
      index = c("parent", "term"),
      vSize = "score", # size ~ significance
      title = paste0(title_label, "| ", direction_label),
      fontsize.labels = c(10, 12),
      align.labels = list(c("center", "top"), c("center", "center"))
    )
  
  reduced
}

plot_ReducedGOterms = function(gsea_df, title_label, ont = "BP")
{
  if(!ont %in% "BP") {return("Please specify correct ontolotgy - BP")}
  
  neg = summarize_and_plot_GOterms_direction(gsea_df %>% filter(NES < 0), title_label, "NES < 0", ont)
  pos = summarize_and_plot_GOterms_direction(gsea_df %>% filter(NES > 0), title_label, "NES > 0", ont)
  
  return(list(neg = neg, pos = pos))
}

# n -> how many metacells to compute
# size -> how many real cells to aggregate expression of to create a metacell
# df must contain cells as rows and genes as columns
# Code adapted from Anne Brunet recent transcriptomic clock paper
bootstrap.metacells = function(df, size=15, n=100, replace="dynamic") {
  set.seed(1)
  metacells = c()
  # If dynamic then only sample with replacement if required due to shortage of cells.
  if (replace == "dynamic") {
    if (nrow(df) <= size) {replace = TRUE} else {replace = FALSE}
  }
  for (i in c(1:n)) {
    batch = df[sample(1:nrow(df), size = size, replace = replace), ]
    metacells = rbind(metacells, colSums(batch))
  }
  colnames(metacells) = colnames(df)
  return(metacells)
}

#' Aggregate single cells into "metacells" (bootstrapped pseudobulk expression profiles),
#' separately per dataset.
#'
#' A metacell is built by summing raw counts across a random sample of
#' `amountOfCells_toAggregate` real cells from the same dataset (see
#' `bootstrap.metacells()`). This reduces the sparsity/dropout noise of single-cell data
#' while still preserving per-dataset variability, which makes downstream pseudobulk-style
#' DESeq2 DEA and glmnet modelling more robust than working on raw single cells directly.
#'
#' @param SO Seurat object, already subset to a single Celltype
#' @param amountOfCells_toAggregate number of real cells summed together per metacell
#' @param N_metacells_generate number of bootstrapped metacells generated per dataset
#' @param min_gene_count minimum total count (summed across all resulting metacells) for
#'   a gene to be kept - removes genes that are essentially unexpressed even after
#'   aggregation, before they add noise to downstream analyses
#' @return a SingleCellExperiment of metacells (log-normalized, low-count genes removed)
create_metacells = function(SO, amountOfCells_toAggregate, N_metacells_generate, min_gene_count = 10) {
  sce = as.SingleCellExperiment(SO)
  raw_counts = counts(sce)
  cell_metadata = colData(sce) %>% as.data.frame()
  
  required_cols = c("dataset", "Celltype", "phenotype", "maturation_stage", "timepoint", "Protocol")
  missing_cols = setdiff(required_cols, colnames(cell_metadata))
  if (length(missing_cols) > 0) {
    stop("create_metacells(): SO is missing required colData column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  celltype = unique(cell_metadata$Celltype)
  # This function assumes SO was already subset to a single Celltype before being passed
  if (length(celltype) != 1) {
    stop("create_metacells(): expected SO to contain exactly one Celltype, found: ", paste(celltype, collapse = ", "))
  }
  
  # One summary row per dataset: numeric columns averaged, character columns take the
  # first value (assumed constant within a dataset, e.g. phenotype/timepoint/Protocol)
  dataset_summary = cell_metadata %>%
    group_by(dataset) %>%
    dplyr::summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      across(where(is.character), ~ dplyr::first(.x))
    )
  
  # Bootstrap-aggregate cells into metacells, independently per dataset, then stack
  # all datasets' metacells together into one counts matrix (metacells x genes)
  pseudobulk_counts = lapply(dataset_summary$dataset, function(dataset_name) {
    cells_in_dataset = rownames(cell_metadata)[cell_metadata$dataset == dataset_name]
    counts_cells_by_genes = t(raw_counts[, cells_in_dataset, drop = FALSE])
    bootstrap.metacells(df = counts_cells_by_genes, size = amountOfCells_toAggregate, n = N_metacells_generate)
  }) %>% do.call(rbind.data.frame, .)
  
  # Expand the per-dataset summary metadata to one row per generated metacell
  n_metacells_total = nrow(dataset_summary) * N_metacells_generate
  dataset_per_metacell = rep(dataset_summary$dataset, each = N_metacells_generate)
  
  metacell_metadata = data.frame(
    dataset = dataset_per_metacell,
    Celltype = celltype,
    phenotype = rep(dataset_summary$phenotype, each = N_metacells_generate),
    maturation_stage = rep(dataset_summary$maturation_stage, each = N_metacells_generate),
    timepoint = rep(dataset_summary$timepoint, each = N_metacells_generate),
    Protocol = rep(dataset_summary$Protocol, each = N_metacells_generate),
    row.names = paste0(dataset_per_metacell, "_Cell_", seq_len(n_metacells_total), "_", celltype)
  )
  
  rownames(pseudobulk_counts) = rownames(metacell_metadata)
  
  sce = SingleCellExperiment(assays = list(counts = t(pseudobulk_counts)), colData = metacell_metadata)
  sce = logNormCounts(sce)
  
  # remove genes that are essentially unexpressed even after aggregating into metacells
  keep = rowSums(assay(sce, "counts")) > min_gene_count
  sce = sce[keep, ]
  
  return(sce)
}

# -----------------------------------------------------------------------------
# adapted from the `geneset` package - the original genGSEA() throws an
# error in our environment, so this is kept as an exact copy to work around that. 
# -----------------------------------------------------------------------------
genGSEA_copy = function (genelist, geneset, padj_method = "BH", p_cutoff = 0.05, 
                         q_cutoff = 0.05, min_gset_size = 10, max_gset_size = 500, 
                         set_seed = FALSE) 
{
  id <- as.character(names(genelist))
  if (missing(geneset)) 
    stop("Please provide gene set...\nWe recommend to use package `geneset` to select available gene set or make new one.")
  genesetType <- geneset$type
  transToSym <- ifelse(genesetType %in% c("enrichrdb", "bp", 
                                          "mf", "cc", "covid19"), TRUE, FALSE)
  org <- geneset$organism
  rareOrg <- TRUE
  tryCatch({
    ens_org <- tryCatch(mapEnsOrg(org), error = function(e) NULL)
    if (is.null(ens_org)) 
      stop("Please make sure the gene types of genelist and geneset are the same.")
    keyType <- gentype(id = id, org = ens_org)
    if (transToSym) {
      id_dat <- suppressMessages(transId(id, "symbol", 
                                         ens_org, unique = T))
      genelist <- genelist[names(genelist) %in% id_dat$input_id]
      names(genelist) <- id_dat$symbol
    }
    else if (keyType != "ENTREZID") {
      id_dat <- suppressMessages(transId(id, "entrezid", 
                                         ens_org, unique = T))
      genelist <- genelist[names(genelist) %in% id_dat$input_id]
      names(genelist) <- id_dat$entrezid
    }
    else if (keyType == "ENTREZID") {
      id_dat <- suppressMessages(transId(id, "symbol", 
                                         ens_org, unique = T)) %>% dplyr::relocate(input_id, 
                                                                                   .after = symbol)
      genelist <- genelist[names(genelist) %in% id_dat$input_id]
    }
  }, error = function(e) {
    ens_org <- org
    transToSym <- FALSE
  })
  fcs <- suppressWarnings(clusterProfiler::GSEA(geneList = genelist, pvalueCutoff = p_cutoff, 
                                                pAdjustMethod = padj_method, minGSSize = min_gset_size, 
                                                maxGSSize = max_gset_size, TERM2GENE = geneset$geneset, 
                                                TERM2NAME = geneset$geneset_name, exponent = 1, eps = 0, 
                                                verbose = FALSE, seed = set_seed, by = "fgsea"))
  if (nrow(as.data.frame(fcs)) == 0) {
    stop("No terms enriched ...")
  } else {
    exponent <- fcs@params[["exponent"]]
    fcs <- fcs %>% as.data.frame() %>% as.enrichdat() %>% 
      dplyr::select(-GeneRatio) %>% dplyr::filter(qvalue < 
                                                    q_cutoff)
  }
  if (!transToSym && !rareOrg) {
    if (keyType != "SYMBOL") {
      if (keyType == "ENTREZID") {
        new_geneID <- get_symbol(fcs$geneID, ens_org)
        new_fcs <- fcs %>% dplyr::mutate(geneID_symbol = new_geneID) %>% 
          dplyr::relocate(geneID_symbol, .after = geneID)
      }
      else {
        old_geneID <- replace_id(id_dat, fcs$geneID)
        new_geneID <- get_symbol(fcs$geneID, ens_org)
        new_fcs <- fcs %>% dplyr::mutate(geneID_symbol = new_geneID) %>% 
          dplyr::mutate(geneID = old_geneID) %>% dplyr::relocate(geneID_symbol, 
                                                                 .after = geneID)
      }
    }
    else {
      old_geneID <- replace_id(id_dat, fcs$geneID)
      new_fcs <- fcs %>% dplyr::mutate(geneID = old_geneID)
    }
  } else if (!rareOrg) {
    if (keyType != "SYMBOL") {
      old_geneID <- replace_id(id_dat, fcs$geneID)
      new_fcs <- fcs %>% dplyr::mutate(geneID_symbol = geneID) %>% 
        dplyr::mutate(geneID = old_geneID) %>% dplyr::relocate(geneID_symbol, 
                                                               .after = geneID)
    } else {
      new_fcs <- fcs
    }
  } else {
    new_fcs <- fcs
  }
  if (!rareOrg) {
    bioc_org <- ensOrg_name %>% dplyr::filter(tolower(latin_short_name) %in% 
                                                geneset$organism) %>% dplyr::pull(bioc_name) %>% 
      stringr::str_to_sentence()
  }
  if (genesetType %in% c("bp", "cc", "mf") && !rareOrg) {
    colnames(new_fcs)[1] = paste0(bioc_org, "_", toupper(genesetType), 
                                  "_ID")
  } else if (genesetType %in% c("bp", "cc", "mf") && rareOrg) {
    colnames(new_fcs)[1] = paste0(org, "_", toupper(genesetType), 
                                  "_ID")
  }
  genelist_df = data.frame(ID = names(genelist), logfc = genelist)
  exponent = data.frame(exponent = exponent)
  org = data.frame(org = org)
  new_geneset <- geneset$geneset
  res <- list(gsea_df = new_fcs, genelist = genelist_df, geneset = new_geneset, 
              exponent = exponent, org = org)
  return(res)
}