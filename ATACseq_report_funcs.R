library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(DT)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(plotly)
library(heatmaply)
library(org.Mm.eg.db)

download_button_png <- function(plot_object, output_name = "plot", width = 12, height = 6, dpi = 300) {
  # Step 1: Save the ggplot object as a temporary PNG file
  temp_file <- tempfile(fileext = ".png")
  ggsave(temp_file, plot = plot_object, width = width, height = height, dpi = dpi)
  
  # Step 2: Convert the image to Base64 encoding
  encoded_img <- base64enc::dataURI(file = temp_file, mime = "image/png")
  
  # Step 3: Create a styled HTML download button
  download_button <- HTML(paste0(
    '<a href="', encoded_img, '" download="', output_name, '.png" ',
    'class="btn btn-primary" style="padding:10px; font-size:16px; text-decoration:none; ',
    'color:white; background-color:#007BFF; border-radius:5px;">',
    '📥 Download ', output_name, '</a>'
  ))
  
  # Return the HTML download button
  return(download_button)
}

summarize_atac_sample_qc <- function(dds, render_table = TRUE) {
  samples <- colnames(dds)
  
  # Identify extra metadata columns to include
  exclude_cols <- c("sample", "sizeFactor", "fastq_1", "fastq_2", "replicate")
  extra_cols <- setdiff(colnames(dds@colData), exclude_cols)
  
  # === 1. OPEN BASES ===
  open_bases_df <- purrr::map_dfr(samples, function(sample) {
    peak_file <- file.path(resolved_seq_dir, paste0(sample, ".mLb.clN_peaks.broadPeak"))
    
    
    
    if (file.exists(peak_file)) {
      peaks <- readr::read_tsv(peak_file, col_names = FALSE, col_types = readr::cols_only(
        X1 = readr::col_character(), X2 = readr::col_double(), X3 = readr::col_double()
      ))
      total_bases <- sum(peaks$X3 - peaks$X2)
      tibble::tibble(Sample = sample, Open_bases = total_bases)
    } else {
      warning(paste("Missing peak file for", sample))
      tibble::tibble(Sample = sample, Open_bases = NA_real_)
    }
  })
  
  # === 2. ALIGNED READS, PEAK COUNT, FACTOR ===
  summary_table <- purrr::map_dfr(samples, function(sample) {
    flagstat_path <- file.path(resolved_seq_dir, paste0(sample, ".mLb.clN.sorted.bam.flagstat"))
    peaks_path <- file.path(resolved_seq_dir, paste0(sample, ".mLb.clN_peaks.broadPeak"))
    
    aligned_reads <- readLines(flagstat_path) %>%
      stringr::str_subset("mapped \\(") %>%
      stringr::str_extract("^\\d+") %>%
      as.numeric()
    
    peak_count <- length(readLines(peaks_path))
    factor <- round(colData(dds)[sample, "sizeFactor", drop = TRUE], 3)
    
    tibble::tibble(
      Sample = sample,
      Aligned_reads = aligned_reads,
      Factor = factor,
      Peaks = peak_count
    )
  }) %>% dplyr::distinct()
  
  # === 3. FRIP SCORES ===
  # MultiQC FRiP raw
  frip_raw <- readr::read_tsv(
    file.path(resolved_seq_dir, "multiqc_mlib_frip_score-plot.txt"),
    col_names = TRUE
  )
  frip_df <- frip_raw %>%
    dplyr::rowwise() %>%
    dplyr::mutate(FRIP = as.numeric(dplyr::cur_data()[[Sample]])) %>%
    dplyr::ungroup() %>%
    dplyr::select(Sample, FRIP) %>%
    dplyr::mutate(FRIP = round(FRIP * 100, 1))
  
  # === 4. JOIN ALL ===
  full_qc_df <- summary_table %>%
    dplyr::left_join(open_bases_df, by = "Sample") %>%
    dplyr::left_join(frip_df, by = "Sample") %>%
    dplyr::select(Sample, Aligned_reads, Open_bases, Factor, Peaks, FRIP) %>%
    dplyr::left_join(as_tibble(colData(dds), rownames = "Sample") %>%
                       dplyr::select(Sample, tidyselect::all_of(extra_cols)),
                     by = "Sample")
  
  # === 5. Optional Render ===
  if (render_table) {
    if (exists("spikeIn_dds")) {
      qc_caption <- "Number of aligned reads and open bases in each sample. ‘Factor’ is the spike-in normalization factor..."
    } else {
      qc_caption <- "Number of aligned reads and open bases in each sample. ‘Factor’ is the TMM normalization factor..."
    }
    
    return(DT::datatable(
      full_qc_df,
      caption = qc_caption,
      options = list(pageLength = 20, dom = 't'),
      rownames = FALSE
    ))
  }
  
  return(full_qc_df)
}


generate_consensus_peak_summary <- function(consensus_bed_path, contrast_list) {
  library(readr)
  library(dplyr)
  library(DT)
  
  # Read consensus BED file
  consensus_peaks <- read_tsv(
    consensus_bed_path,
    col_names = c("chr", "start", "end", "name", "score", "strand"),
    col_types = cols(.default = col_character())
  ) %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    mutate(width = end - start)
  
  # Number of consensus peaks and their average size
  n_peaks <- nrow(consensus_peaks)
  avg_size <- mean(consensus_peaks$width)
  
  # Create the summary for each contrast (same values repeated)
  summary_df <- purrr::map_dfr(contrast_list, function(pair) {
    tibble(
      Test = pair[[1]],
      Ctrl = pair[[2]],
      `Number of peaks` = format(n_peaks, big.mark = ","),
      `Avg size in Test` = paste0(round(avg_size, 1), " bp"),
      `Avg size in Ctrl` = paste0(round(avg_size, 1), " bp")
    )
  })
  
  # Render datatable
  datatable(summary_df,
            caption = "",  # You can add your own caption here
            rownames = FALSE,
            options = list(pageLength = 10, dom = 't'))
}

parse_contrasts <- function(x) {
  # Determine whether `x` is a file path or a string of contrasts
  if (file.exists(x)) {
    contrast_lines <- readLines(x, warn = FALSE)
  } else {
    # Assume comma-separated list
    contrast_lines <- unlist(strsplit(x, ","))
  }
  
  # Clean and split each contrast
  contrast_lines <- trimws(contrast_lines)
  contrast_lines <- contrast_lines[contrast_lines != ""]
  
  contrast_list <- lapply(contrast_lines, function(line) {
    parts <- strsplit(line, "_vs_")[[1]]
    if (length(parts) != 2) stop(paste("Invalid contrast:", line))
    parts
  })
  
  return(contrast_list)
}


# Usage
#contrast_file <- report_params$contrasts
#contrast_list <- parse_contrasts(contrast_file)

convert_all_chr_to_ucsc <- function(results_list) {
  convert_to_ucsc_chr <- function(chr_vec) {
    chr_vec <- gsub("^chr", "", chr_vec)  # strip existing 'chr' if present
    chr_ucsc <- ifelse(chr_vec %in% c("MT", "M"), "chrM", paste0("chr", chr_vec))
    return(chr_ucsc)
  }
  
  for (contrast in names(results_list)) {
    if (!is.null(results_list[[contrast]]$table$Chr)) {
      results_list[[contrast]]$table$Chr <- convert_to_ucsc_chr(results_list[[contrast]]$table$Chr)
    }
  }
  
  return(results_list)
}

# Path to the bash script
build_hub_script <- "build_ucsc_trackhub.sh"

generate_enrichment_plot_atac <- function(gene_lists, de_results_df, universe_entrez, ont_category, significance_threshold = 0.05, top_n = 10, annotation_db) {
  
  # Load the annotation object (same as before)
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db)) 
  names(gene_lists) <- gsub("efit_|_results_df","",names(gene_lists))
  
  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)
  
  # Map Ensembl IDs in de_results_df to ENTREZID
  de_results_df <- de_results_df %>%
    mutate(ENTREZID = mapIds(
      annotation_obj, 
      keys = Entrez.ID,  # Use Ensembl IDs from 'Entrez.ID' column
      column = "ENTREZID",  # Map to Entrez IDs
      keytype = "ENSEMBL",  # Correct keytype is ENSEMBL for Ensembl IDs
      multiVals = "first"
    ))
  # Prepare data for compareCluster
  data <- bind_rows(lapply(contrast_order, function(contrast) {
    genes <- gene_lists[[contrast]]
    if (length(genes) == 0 || all(is.na(genes))) {
      return(NULL)  # Skip empty contrasts
    }
    data.frame(
      Entrez = genes,
      Contrast = contrast,
      stringsAsFactors = FALSE
    )
  }))
  
  
  # Ensure no duplicates and valid formatting
  data <- data %>%
    distinct(Entrez, Contrast, .keep_all = TRUE)  # Remove duplicate Entrez entries within each contrast
  data$Entrez <- as.character(data$Entrez)  # Ensure Entrez IDs are characters
  
  # Run GO enrichment using compareCluster() for ATAC-seq contrasts
  formula_res <- compareCluster(
    Entrez ~ Contrast,
    data = data,
    fun = "enrichGO",
    universe = na.omit(universe_entrez),
    OrgDb = annotation_obj,
    keyType = "ENTREZID",
    ont = ont_category,
    pvalueCutoff = significance_threshold
  )
  # Handle no results case
  if (is.null(formula_res) || nrow(formula_res@compareClusterResult) == 0) {
    formula_res <- new("compareClusterResult",
                       compareClusterResult = data.frame(
                         Cluster = factor(),
                         ID = character(),
                         Description = character(),
                         GeneRatio = character(),
                         BgRatio = character(),
                         pvalue = numeric(),
                         p.adjust = numeric(),
                         qvalue = numeric(),
                         geneID = character(),
                         Count = integer(),
                         stringsAsFactors = FALSE
                       ))
  }
  # Ensure clusters are ordered correctly
  formula_res@compareClusterResult$Cluster <- factor(
    formula_res@compareClusterResult$Cluster,
    levels = contrast_order
  )
  
  # Filter by significance threshold
  filtered_results <- subset(
    formula_res@compareClusterResult,
    p.adjust <= significance_threshold
  )
  # Handle the special case: no enrichment found
  if (nrow(filtered_results) == 0) {
    message_plot <- ggplot() +
      annotate("text", x = 1, y = 1, label = paste0("No significant GO enrichment found\n(", ont_category, ")"), size = 6, hjust = 0.5) +
      theme_void() +
      ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = ""))
    
    interactive_plot <- ggplotly(message_plot)
    static_plot <- message_plot
    
    return(list(
      interactive_plot = interactive_plot,
      static_plot = static_plot,
      go_results = NULL
    ))
  }
  
  # Step 1: Gather all Entrez IDs from filtered_results
  all_entrez_ids <- unique(unlist(strsplit(filtered_results$geneID, "/")))
  
  # Step 2: Map all Entrez IDs to gene symbols in one go
  gene_symbols_batch <- mapIds(annotation_obj,
                               keys = all_entrez_ids,  # Pass all Entrez IDs at once
                               column = "SYMBOL",
                               keytype = "ENTREZID",
                               multiVals = "first") %>%
    na.omit()  # Remove any NAs from the mapping
  
  # Step 3: Create a named vector of gene symbols (Entrez IDs as keys)
  gene_symbols_vector <- gene_symbols_batch[all_entrez_ids]  # Map each Entrez ID to its symbol
  
  # Step 4: Assign gene symbols based on top peaks per contrast (row-wise)
  # Step 1: Pre-filter de_results_df once for relevant contrasts
  # Step 1: Pre-filter de_results_df once for relevant contrasts
  de_results_df_filtered <- de_results_df %>%
    filter(!is.na(ENTREZID) & ENTREZID != "")  # Remove rows with NA or empty ENTREZID
  
  # Step 1: Load the parallel package
  library(parallel)
  
  # Step 2: Create a function that assigns gene symbols to a list of Entrez IDs (same as before)
  assign_gene_symbols <- function(peak_list, contrast_full) {
    # Extract base contrast and direction
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)
    
    # Split the peak list into Entrez IDs
    entrez_ids <- unlist(strsplit(peak_list, "/"))
    
    # Filter de_results_df based on contrast and direction
    de_sub <- de_results_df_filtered %>%
      filter(grepl(contrast_base, Contrast),  # Match the base name
             ENTREZID %in% entrez_ids,         # Ensure ENTREZID is valid
             case_when(
               direction == "up" ~ logFC > 0,
               direction == "down" ~ logFC < 0
             ))
    
    # Sort by FDR and get the top `n` peaks
    top_peaks <- de_sub %>%
      arrange(FDR) %>%
      slice_head(n = 20) %>%
      pull(ENTREZID)
    top_peaks <- unlist(top_peaks)
    
    # Map Entrez IDs to gene symbols
    gene_symbols <- gene_symbols_vector[top_peaks]
    
    # Return gene symbols as a string (with line breaks)
    if (length(gene_symbols) > 0) {
      paste(unique(gene_symbols), collapse = "<br>")
      
    } else {
      NA_character_
    }
  }
  
  # Step 3: Parallelize the function across all rows in filtered_results using mclapply
  # Determine number of available cores
  num_cores <- 5  # Leave one core free for the system
  
  # Parallelize using mclapply (unlike lapply, mclapply runs in parallel)
  filtered_results$GeneSymbols <- mclapply(seq_len(nrow(filtered_results)), function(i) {
    peak_list <- filtered_results$geneID[i]
    contrast_full <- as.character(filtered_results$Cluster[i])
    assign_gene_symbols(peak_list, contrast_full)
  }, mc.cores = num_cores)  # Specify the number of cores to use
  
  # Save filtered GO results for download
  download_go_results <- filtered_results %>%
    select(Cluster, Description, p.adjust, GeneSymbols, everything())
  
  # Identify top `n` GO terms for plotting
  top_GO_terms <- filtered_results %>%
    group_by(Cluster) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  # Filter for plotting (only top GO terms)
  formula_res@compareClusterResult <- filtered_results %>%
    filter(Description %in% top_GO_terms)
  
  # Convert GeneRatio to numeric
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  # Reorder GO terms using hierarchical clustering (same as before)
  reorder_GO_terms <- function(df) {
    term_matrix <- table(df$Description, df$Cluster)  
    
    if (nrow(term_matrix) > 1 && length(unique(df$Description)) > 1) {
      term_dist <- dist(term_matrix, method = "binary")  
      term_hclust <- hclust(term_dist, method = "ward.D2")  
      term_order <- rownames(term_matrix)[term_hclust$order]
    } else {
      term_order <- unique(df$Description)  # Ensure a valid factor level
    }
    
    df$Description <- factor(df$Description, levels = term_order)  
    return(df)
  }
  
  formula_res@compareClusterResult <- reorder_GO_terms(formula_res@compareClusterResult)
  
  # Define GeneRatio bins and size mapping
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(formula_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)  
  bin_labels <- c("≤0.01", "0.01 - 0.05", "0.05 - 0.10", "≥0.10")
  
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))
  
  size_mapping <- c("≤0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, "≥0.10" = 8)
  
  # Correct the mutate statement
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      plot_label = ifelse(
        sapply(strsplit(as.character(Description), " "), length) > 6,
        sapply(strsplit(as.character(Description), " "), function(words) {
          paste(c(words[1:3], "...", tail(words, 3)), collapse = " ")
        }),
        as.character(Description)
      )
    ) %>%
    mutate(
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",  
        "GO Term: ", as.character(Description), "<br>",  
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)  
      )
    )
  
  # Create interactive plot
  p <- ggplot(formula_res@compareClusterResult, aes(
    x = Cluster, 
    y = plot_label, 
    size = GeneRatioCategory,  
    color = p.adjust,
    text = tooltip_text  
  )) +
    geom_point(alpha = 0.8) +  
    scale_size_manual(name = "Gene Ratio", values = size_mapping) +  
    scale_color_gradient(low = "red", high = "blue", 
                         limits = c(min(formula_res@compareClusterResult$p.adjust, na.rm = TRUE), 
                                    max(formula_res@compareClusterResult$p.adjust, na.rm = TRUE)),
                         name = "p.adjust") +  
    guides(color = guide_colorbar(title = "p.adjust")) +  
    ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    xlab("DE Gene list") +
    ylab("GO Term") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  # **Convert to interactive plot (ensure tooltips work)**
  interactive_plot <- ggplotly(p, tooltip = "text") %>%  
    layout(legend = list(title = list(text = "Gene Ratio")))
  
  # **Static High-Resolution Plot (for manuscript)**
  static_plot <- dotplot(formula_res, showCategory = top_n) +
    ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Return interactive plot, static plot, and GO results
  return(list(
    interactive_plot = interactive_plot, 
    static_plot = static_plot, 
    go_results = download_go_results
  ))
}

generate_kegg_enrichment_plot_atac <- function(gene_lists, de_results_df, universe_entrez, significance_threshold = 0.05, top_n = 10, annotation_db) {
  
  # Load the annotation object (same as before)
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db)) 
  names(gene_lists) <- gsub("efit_|_results_df","",names(gene_lists))
  
  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)
  
  # Map Ensembl IDs in de_results_df to ENTREZID
  de_results_df <- de_results_df %>%
    mutate(ENTREZID = mapIds(
      annotation_obj, 
      keys = Entrez.ID,  # Use Ensembl IDs from 'Entrez.ID' column
      column = "ENTREZID",  # Map to Entrez IDs
      keytype = "ENSEMBL",  # Correct keytype is ENSEMBL for Ensembl IDs
      multiVals = "first"
    ))
  
  # Prepare data for compareCluster
  data <- bind_rows(lapply(contrast_order, function(contrast) {
    data.frame(
      Entrez = gene_lists[[contrast]],
      Contrast = contrast
    )
  }))
  
  # Ensure no duplicates and valid formatting
  data <- data %>%
    distinct(Entrez, Contrast, .keep_all = TRUE)  # Remove duplicate Entrez entries within each contrast
  data$Entrez <- as.character(data$Entrez)  # Ensure Entrez IDs are characters
  
  # Run KEGG enrichment using compareCluster() for ATAC-seq contrasts
  kegg_res <- compareCluster(
    Entrez ~ Contrast,
    data = data,
    fun = "enrichKEGG",
    universe = na.omit(universe_entrez),
    organism = "hsa",  # Assuming human (Homo sapiens) for KEGG enrichment
    keyType = "kegg",
    pvalueCutoff = significance_threshold
  )
  
  # Ensure clusters are ordered correctly
  kegg_res@compareClusterResult$Cluster <- factor(
    kegg_res@compareClusterResult$Cluster,
    levels = contrast_order
  )
  
  # **Filter only by significance threshold FIRST (to preserve all significant results)**
  filtered_results <- subset(
    kegg_res@compareClusterResult,
    p.adjust <= significance_threshold
  )
  
  # Step 1: Gather all Entrez IDs from filtered_results
  all_entrez_ids <- unique(unlist(strsplit(filtered_results$geneID, "/")))
  
  # Step 2: Map all Entrez IDs to gene symbols in one go
  gene_symbols_batch <- mapIds(annotation_obj,
                               keys = all_entrez_ids,  # Pass all Entrez IDs at once
                               column = "SYMBOL",
                               keytype = "ENTREZID",
                               multiVals = "first") %>%
    na.omit()  # Remove any NAs from the mapping
  
  # Step 3: Create a named vector of gene symbols (Entrez IDs as keys)
  gene_symbols_vector <- gene_symbols_batch[all_entrez_ids]  # Map each Entrez ID to its symbol
  
  # Step 4: Assign gene symbols based on top peaks per contrast (parallelized)
  
  # Step 1: Pre-filter de_results_df once for relevant contrasts
  de_results_df_filtered <- de_results_df %>%
    filter(!is.na(ENTREZID) & ENTREZID != "")  # Remove rows with NA or empty ENTREZID
  
  # Step 2: Parallelize the function across all rows in filtered_results using mclapply
  library(parallel)
  num_cores <- 5  # Set the number of cores (adjust as needed)
  
  # Create a function to assign gene symbols to a list of Entrez IDs
  assign_gene_symbols <- function(peak_list, contrast_full) {
    # Extract base contrast and direction
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)
    
    # Split the peak list into Entrez IDs
    entrez_ids <- unlist(strsplit(peak_list, "/"))
    
    # Filter de_results_df based on contrast and direction
    de_sub <- de_results_df_filtered %>%
      filter(grepl(contrast_base, Contrast),  # Match the base name
             ENTREZID %in% entrez_ids,         # Ensure ENTREZID is valid
             case_when(
               direction == "up" ~ logFC > 0,
               direction == "down" ~ logFC < 0
             ))
    
    # Sort by FDR and get the top `n` peaks
    top_peaks <- de_sub %>%
      arrange(FDR) %>%
      slice_head(n = 20) %>%
      pull(ENTREZID)
    top_peaks <- unlist(top_peaks)
    
    # Map Entrez IDs to gene symbols
    gene_symbols <- gene_symbols_vector[top_peaks]
    
    # Return gene symbols as a string (with line breaks)
    if (length(gene_symbols) > 0) {
      paste(unique(gene_symbols), collapse = "<br>")
    } else {
      NA_character_
    }
  }
  
  # Parallelize using mclapply (unlike lapply, mclapply runs in parallel)
  filtered_results$GeneSymbols <- mclapply(seq_len(nrow(filtered_results)), function(i) {
    peak_list <- filtered_results$geneID[i]
    contrast_full <- as.character(filtered_results$Cluster[i])
    assign_gene_symbols(peak_list, contrast_full)
  }, mc.cores = num_cores)  # Specify the number of cores to use
  
  # Save filtered KEGG results for download
  download_kegg_results <- filtered_results %>%
    select(Cluster, Description, p.adjust, GeneSymbols, everything())  # Reorder columns
  
  # **Identify top `n` KEGG pathways for plotting**
  top_KEGG_terms <- filtered_results %>%
    group_by(Cluster) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  # **Filter for plotting (only keeping top KEGG pathways)**
  kegg_res@compareClusterResult <- filtered_results %>%
    filter(Description %in% top_KEGG_terms)
  
  # Convert GeneRatio to numeric
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  # **Restore KEGG Pathway Ordering Using Hierarchical Clustering**
  reorder_KEGG_terms <- function(df) {
    term_matrix <- table(df$Description, df$Cluster)  
    
    if (nrow(term_matrix) > 1) {
      term_dist <- dist(term_matrix, method = "binary")  
      term_hclust <- hclust(term_dist, method = "ward.D2")  
      term_order <- rownames(term_matrix)[term_hclust$order]
    } else {
      term_order <- df$Description  
    }
    
    df$Description <- factor(df$Description, levels = term_order)  
    return(df)
  }
  
  kegg_res@compareClusterResult <- reorder_KEGG_terms(kegg_res@compareClusterResult)
  
  # **Define GeneRatio bins and size mapping**
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(kegg_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)  
  bin_labels <- c("≤0.01", "0.01 - 0.05", "0.05 - 0.10", "≥0.10")
  
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))
  
  # **Define reference dot sizes**
  size_mapping <- c("≤0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, "≥0.10" = 8)
  
  # **Ensure all necessary variables are characters/numeric**
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",  
        "KEGG Pathway: ", as.character(Description), "<br>",  
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)  
      )
    )
  
  # **Create ggplot object**
  p <- ggplot(kegg_res@compareClusterResult, aes(
    x = Cluster, 
    y = Description, 
    size = GeneRatioCategory, 
    color = p.adjust,
    text = tooltip_text
  )) +
    geom_point(alpha = 0.8) +  
    scale_size_manual(name = "Gene Ratio", values = size_mapping) +  
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +  
    guides(color = guide_colorbar(title = "p.adjust")) +  
    ggtitle("KEGG Pathway Enrichment") +
    xlab("DE Gene list") +
    ylab("KEGG Pathway") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) 
  
  interactive_plot <- ggplotly(p, tooltip = "text") %>%  
    layout(legend = list(
      title = list(text = "Gene Ratio"),  # Ensure Gene Ratio is labeled
      colorbar = list(title = "p.adjust") # Ensure p.adjust is labeled for color
    ))
  
  
  # **Static KEGG Plot**
  static_plot <- dotplot(kegg_res, showCategory = top_n) +
    ggtitle("KEGG Pathway Enrichment") +
    theme_minimal()
  
  return(list(interactive_plot = interactive_plot, static_plot = static_plot, kegg_results = download_kegg_results))
}

write_bigbed_de_peaks <- function(results_list, report_params, chrom_sizes_file, bedToBigBed_path) {
  public_web_dir <- "/orange/cancercenter-dept/web/public/BCB-SR/trackhubs"
  
  for (contrast in names(results_list)) {
    contrast_clean <- gsub("\\.", "-", contrast)
    contrast_clean <- gsub("_", "-", contrast_clean)
    contrast_clean <- gsub("-vs-", "_vs_", contrast_clean)
    
    trackhub_dir <- file.path(public_web_dir, report_params$seqID, contrast_clean)
    if (!dir.exists(trackhub_dir)) dir.create(trackhub_dir, recursive = TRUE)
    
    result_table <- results_list[[contrast]]$table
    de_table <- subset(result_table, FDR < 0.05 & abs(logFC) > 1 & logCPM > 2 & Gene.Name != "")
    
    # Read valid UCSC chromosome names
    valid_chroms <- read.table(chrom_sizes_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
    
    # Filter to valid chromosomes
    de_table <- de_table[de_table$Chr %in% valid_chroms, ]
    
    if (nrow(de_table) == 0) next
    
    # Safe fallback for missing names
    peak_names <- ifelse(
      is.na(de_table$Gene.Name) | de_table$Gene.Name == "",
      if ("interval" %in% colnames(de_table)) {
        de_table$interval
      } else {
        paste0(de_table$Chr, "_", de_table$Start, "_", de_table$End)
      },
      de_table$Gene.Name
    )
    
    # Construct BED dataframe
    bed_df <- data.frame(
      chrom      = de_table$Chr,
      chromStart = de_table$Start - 1,  # BED is 0-based
      chromEnd   = de_table$End,
      name       = peak_names,
      score      = if ("Peak.Score" %in% colnames(de_table)) as.integer(de_table$Peak.Score) else 0L,
      strand     = de_table$Strand,
      stringsAsFactors = FALSE
    )
    
    bed_file <- file.path(trackhub_dir, paste0(contrast_clean, "_DA_peaks.bed"))
    write.table(bed_df, bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Sort BED
    sort_cmd <- sprintf("LC_COLLATE=C sort -k1,1 -k2,2n %s -o %s", bed_file, bed_file)
    system(sort_cmd)
    
    # Convert to BigBed
    bb_file <- file.path(trackhub_dir, paste0(contrast_clean, "_DA_peaks.bb"))
    cmd <- sprintf("%s %s %s %s", bedToBigBed_path, bed_file, chrom_sizes_file, bb_file)
    system(cmd, intern = TRUE)
    
    file.remove(bed_file)
  }
}

find_variable_for_contrast <- function(group1, group2, metadata) {
  for (var in colnames(metadata)) {
    vals <- make.names(unique(as.character(metadata[[var]])))
    if (all(make.names(c(group1, group2)) %in% vals)) return(var)
  }
  stop("No matching variable found for contrast.")
}

# Read the report parameters from the text file
parse_params <- function(filepath) {
  lines <- readLines(filepath)
  param_list <- list()
  
  for (line in lines) {
    if (grepl("^--", line)) {
      parts <- strsplit(line, " ", fixed = TRUE)[[1]]
      key <- gsub("^--", "", parts[1])  # Remove leading --
      value <- paste(parts[-1], collapse = " ")
      value <- gsub('^"|"$', '', value)  # Remove surrounding quotes
      
      # Convert to numeric if it looks like a number
      if (!is.na(suppressWarnings(as.numeric(value)))) {
        value <- as.numeric(value)
      }
      
      param_list[[key]] <- value
    }
  }
  
  return(param_list)
}

get_log_matrix <- function(dge) {
  if (!is.null(dge$E_corrected)) {
    return(dge$E_corrected)
  } else if (!is.null(dge$counts)) {
    return(cpm(dge, log = TRUE))
  } else {
    stop("Input DGE object must contain either 'counts' or 'E_corrected'.")
  }
}

library(plotly)
library(RColorBrewer)
library(car)

plot_pca <- function(dge, title, grp_var=report_params$group_var, show_legend = TRUE, combine_plots = FALSE) {
  # Extract log-transformed CPM values
  PCA_DATA <- t(get_log_matrix(dge))
  rownames(PCA_DATA) <- dge$samples$sample_id
  numsonly <- as.data.frame(PCA_DATA)
  numsonly <- as.data.frame(lapply(numsonly, as.numeric))  # Ensure numeric values
  
  # Perform PCA
  pca_res <- prcomp(numsonly, center = TRUE, scale. = FALSE, rank. = 2)
  scores <- pca_res$x
  
  # Extract group information
  group_labels <- dge$samples[[grp_var]]  
  colors <- RColorBrewer::brewer.pal(length(unique(group_labels)), "Set1")
  
  # Initialize plotly object
  fig <- plot_ly()
  
  # Function to plot ellipses
  plot_ellipse_function <- function(target) {
    target_data <- as.data.frame(scores[group_labels == target, 1:2])
    target_color <- colors[which(unique(group_labels) == target)]
    
    ellipse_coords <- car::dataEllipse(target_data$PC1, target_data$PC2, 
                                       levels = 0.68, plot.points = FALSE, add = TRUE, draw = FALSE)
    
    fig <<- fig %>%
      add_polygons(x = ellipse_coords[,1], y = ellipse_coords[,2],
                   line = list(color = target_color, dash = "dot"),
                   fillcolor = target_color, opacity = 0.3,
                   showlegend = FALSE, hoverinfo = "skip")
  }
  
  # Function to plot points
  plot_points_function <- function(target) {
    target_data <- as.data.frame(scores[group_labels == target, 1:2])
    rownames(target_data) <- rownames(PCA_DATA)[group_labels == target]
    target_color <- colors[which(unique(group_labels) == target)]
    
    text_data <- as.character(rownames(target_data))
    PC1 <- as.numeric(target_data$PC1)
    PC2 <- as.numeric(target_data$PC2)
    
    fig <<- fig %>%
      add_trace(data = target_data, x = ~PC1, y = ~PC2, 
                type = "scatter", mode = "markers", name = target,
                marker = list(color = target_color),
                text = ~text_data, hoverinfo = "text", showlegend = show_legend)
  }
  
  # Add ellipses and points
  invisible(lapply(unique(group_labels), plot_ellipse_function))
  invisible(lapply(unique(group_labels), plot_points_function))
  
  # Add title
  # Get % variance explained
  percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  # Add axis titles with % variance
  fig <- fig %>% layout(
    title = title,
    xaxis = list(title = paste0("PC1 (", percentVar[1], "%)")),
    yaxis = list(title = paste0("PC2 (", percentVar[2], "%)"))
  )
  
  return(fig)
}