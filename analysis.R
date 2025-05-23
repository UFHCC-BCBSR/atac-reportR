# Annotate peaks
row.names(peaks_anno) <- peaks_anno$interval
dge$genes <- peaks_anno[rownames(dge), ]

# Run differential analysis per contrast
# Extract contrasts from strings like "1. Heat vs Control"
contrast_list <- lapply(contrast_strings, function(x) {
  # Remove leading "1.", "2.", etc.
  x_clean <- sub("^\\d+\\.\\s*", "", x)
  # Split on "vs" and trim whitespace
  parts <- trimws(unlist(strsplit(x_clean, "\\s+vs\\s+", perl = TRUE)))
  make.names(parts)  # Ensure safe for model.matrix
})
rownames(dge$samples) <- colnames(dge) <- make.names(colnames(dge))

results_list <- list()
contrast_list <- lapply(contrast_list, function(x) gsub("-", "_", x))

for (i in seq_along(contrast_list)) {
  group1 <- make.names(contrast_list[[i]][1])
  group2 <- make.names(contrast_list[[i]][2])
  contrast_name <- paste0(group1, "_vs_", group2)
  
    # Identify which column (variable) in dge$samples contains this contrast
      var_match <- NULL
    for (var in colnames(dge$samples)) {
        vals <- make.names(unique(as.character(dge$samples[[var]])))
        if (all(c(group1, group2) %in% vals)) {
            var_match <- var
            break
          }
      }
    if (is.null(var_match)) stop(paste("Could not match contrast:", contrast_name))
  
    # Subset samples for this contrast
      dge$samples$Condition <- make.names(as.character(dge$samples[[var_match]]))
    relevant_samples <- which(dge$samples$Condition %in% c(group1, group2))
    dge_contrast <- dge[, relevant_samples]
  
  # Redefine design for this contrast
  dge_contrast$samples$Condition <- droplevels(factor(dge_contrast$samples$Condition))
  design_contrast <- model.matrix(~0 + Condition, data = dge_contrast$samples)
  colnames(design_contrast) <- levels(dge_contrast$samples$Condition)
  
  # Filter peaks
  keep <- filterByExpr(dge_contrast, design_contrast, min.count = report_params$min_count_for_filtering, min.prop = report_params$min_prop_for_filtering)
  dge_contrast <- dge_contrast[keep, , keep.lib.sizes = TRUE]
  
  # Model fitting
  dge_contrast <- estimateDisp(dge_contrast, design_contrast)
  fit_contrast <- glmQLFit(dge_contrast, design_contrast)
  
  # Define contrast
  contrast_vec <- numeric(ncol(design_contrast))
  contrast_vec[which(colnames(design_contrast) == group1)] <- 1
  contrast_vec[which(colnames(design_contrast) == group2)] <- -1
  
  # Run test and store results
  result <- glmQLFTest(fit_contrast, contrast = contrast_vec)
  results_list[[contrast_name]] <- topTags(result, n = Inf)
}


