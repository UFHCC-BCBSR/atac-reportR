server <- function(input, output, session) {
  app_base_dir <- getwd()
  params_path_global <- reactiveVal(NULL)
  sample_sheet_path <- reactiveVal(NULL)
  sample_sheet_temp_file <- reactiveVal(NULL)
  current_params <- reactiveVal(list())
  report_ready <- reactiveVal(FALSE)
  params_generated <- reactiveVal(FALSE)
  contrasts_path <- reactiveVal(NULL)
  host_appdir <- Sys.getenv("HOST_APPDIR")
  if (host_appdir == "") host_appdir <- getwd()  # fallback if unset
  data_dir <- reactive({
    req(input$seqID)
    Sys.glob(file.path('/blue/cancercenter-dept/privapps/data/atac', '*', input$seqID))
  })  
  output$sample_sheet_source <- renderText({
    req(input$seqID)
    if (!is.null(sample_sheet_path())) {
      return(paste("✅ Using uploaded sample sheet:", sample_sheet_path()))
    } else if (!is.null(input$sample_sheet) && file.exists(input$sample_sheet)) {
      return(paste("✅ Using sample sheet from text field:", input$sample_sheet))
    } else {
      fallback <- Sys.glob(file.path('/blue/cancercenter-dept/privapps/data/atac', '*', input$seqID, 'samplesheet.valid.csv'))[1]
      if (!is.na(fallback) && file.exists(fallback)) {
        return(paste("✅ Using fallback sample sheet from derived path:", fallback))
      }
    }
    return("❌ No valid sample sheet found.")
  })
  
  
  dq <- function(x) paste0('"', x, '"')
  
  output$download_sample_sheet <- downloadHandler(
    filename = function() paste0(input$seqID, "_sample_sheet.csv"),
    content = function(file) {
      df <- tryCatch(sample_sheet_data(), error = function(e) NULL)
      if (is.null(df)) {
        stop("No sample sheet available to download.")
      }
      
      write.csv(df, file, row.names = FALSE)
    }
  )
  observeEvent(input$contrasts, {
    req(input$contrasts)
    dest <- file.path(data_dir(), paste0(input$seqID, "_contrasts.csv"))
    file.copy(input$contrasts$datapath, dest, overwrite = TRUE)
    contrasts_path(dest)
  })
  
  sample_sheet_data <- reactive({
    # Priority 1: user uploaded a new file
    if (!is.null(sample_sheet_path())) {
      return(read.csv(sample_sheet_path(), stringsAsFactors = FALSE))
    }
    
    # Priority 2: user manually typed a path
    if (!is.null(input$sample_sheet) && input$sample_sheet != "") {
      if (file.exists(input$sample_sheet)) {
        return(read.csv(input$sample_sheet, stringsAsFactors = FALSE))
      }
    }
    
    # Priority 3: fallback from nf-core output dir
    fallback <- Sys.glob(file.path('/blue/cancercenter-dept/privapps/data/atac', '*', input$seqID, 'samplesheet.valid.csv'))[1]
    if (!is.na(fallback) && file.exists(fallback)) {
      return(read.csv(fallback, stringsAsFactors = FALSE))
    }
    
    return(NULL)
  })
  
  
  validate_sample_sheet <- function() {
    df <- tryCatch(sample_sheet_data(), error = function(e) NULL)
    if (is.null(df)) {
      output$sample_validation_status <- renderText("❌ No sample sheet available from upload, params file, or nf-core output.")
      report_ready(FALSE)
      return()
    }
    messages <- c()
    if (!"sample_id" %in% colnames(df)) {
      messages <- c(messages, "⚠️ No 'sample_id' column found. The 'sample' column will be used instead.")
    }
    contrast_lines <- NULL
    if (!is.null(input$contrasts)) {
      # input$contrasts is a fileInput; read from its $datapath
      contrast_lines <- trimws(readLines(input$contrasts$datapath, warn = FALSE))
      contrast_lines <- contrast_lines[contrast_lines != ""]
    }
    if (!is.null(contrast_lines) && length(contrast_lines) > 0) {
      for (contrast in contrast_lines) {
        parts <- unlist(strsplit(contrast, "_vs_"))
        found <- any(sapply(df, function(col) all(parts %in% unique(col))))
        if (!found) {
          messages <- c(messages, paste0("❌ No column found with both contrast values: ", paste(parts, collapse = " vs ")))
        }
      }
    }
    if (length(messages) == 0) {
      output$sample_validation_status <- renderText("✅ Sample sheet loaded and validated.")
      report_ready(TRUE)
    } else {
      output$sample_validation_status <- renderText(paste(messages, collapse = "\n"))
      report_ready(FALSE)
    }
  }
  observeEvent(input$validate_sheet, {
    validate_sample_sheet()
  })
  
  
  observeEvent(input$upload_sample_sheet, {
    tryCatch({
      req(input$upload_sample_sheet)
      df <- read.csv(input$upload_sample_sheet$datapath, stringsAsFactors = FALSE)
      path <- file.path(data_dir(), paste0(input$seqID, "edited_sample_sheet.csv"))
      write.csv(df, path, row.names = FALSE)
      sample_sheet_path(path)
      report_ready(FALSE)
      validate_sample_sheet() # << automatically validate the new upload!
    }, error = function(e) {
      output$sample_validation_status <- renderText(paste("❌ Error uploading sample sheet:", e$message))
    })
  })
  observe({
    # This will trigger when sample_sheet_path() and input$seqID change
    if (is.null(sample_sheet_path())) {
      fallback <- Sys.glob(file.path('/blue/cancercenter-dept/privapps/data/atac', '*', input$seqID, 'samplesheet.valid.csv'))[1]
      if (!is.na(fallback) && file.exists(fallback)) {
        validate_sample_sheet()
      }
    }
  })
  
  generate_param_lines <- function(input) {
    params <- list(
      paste("--report_title", dq(input$report_title)),
      paste("--group_var", dq(input$group_var)),
      paste("--organism", dq(input$organism)),
      paste("--annotation_db", dq(input$annotation_db)),
      paste("--PI", dq(input$PI)),
      paste("--Institution", dq(input$Institution)),
      paste("--Department", dq(input$Department)),
      paste("--Study_Contact", dq(input$Study_Contact)),
      paste("--Project_Title", dq(input$Project_Title)),
      paste("--Study_Summary", dq(input$Study_Summary)),
      paste("--Sample_Types", dq(input$Sample_Types)),
      paste("--Organism", dq(input$Organism_Text)),
      paste("--Analysis_Goals", dq(input$Analysis_Goals)),
      paste("--Report_Prepared_By", dq(input$Report_Prepared_By)),
      paste("--Report_Reviewed_By", dq(input$Report_Reviewed_By)),
      paste("--raw_seq_URL", dq(input$raw_seq_URL)),
      paste("--multiqc_url", dq(input$multiqc_url)),
      paste("--seqID", dq(input$seqID)),
      paste("--min_count_for_filtering", dq(input$min_count_for_filtering)),
      paste("--min_prop_for_filtering", dq(input$min_prop_for_filtering)),
      paste("--test_mode", dq(ifelse(input$test_mode, "TRUE", "FALSE")))
    )
    if (!is.null(contrasts_path())) {
      params <- append(params, paste("--contrasts", dq(contrasts_path())))
    }
    if (!is.null(input$nfcore_spikein_dir)) {
      params <- append(params, paste("--nfcore_spikein_dir", dq(input$nfcore_spikein_dir)))
    }
    if (!is.null(sample_sheet_path())) {
      params <- append(params, paste("--sample_sheet", dq(sample_sheet_path())))
    } else if (input$sample_sheet == "" || is.null(input$sample_sheet)) {
      fallback_path <- Sys.glob(file.path('/blue/cancercenter-dept/privapps/data/atac', '*', input$seqID, 'samplesheet.valid.csv'))[1]
      if (!is.na(fallback_path) && file.exists(fallback_path)) {
        params <- append(params, paste("--sample_sheet", dq(fallback_path)))
      }
    }
    return(params)
  }
  
  observeEvent(input$generate_params, {
    tryCatch({
      params <- generate_param_lines(input)
      # --- Check that group_var is a column in the sample sheet ---
      sample_df <- sample_sheet_data()
      group_var <- input$group_var
      
      if (is.null(sample_df) || !(group_var %in% colnames(sample_df))) {
        output$status <- renderText(paste0("❌ Error: Group variable '", group_var, "' is not a column in the sample sheet."))
        return()
      }
      
      
      param_data <- list()
      for (line in params) {
        parts <- strsplit(gsub("^--", "", line), " ", fixed = TRUE)[[1]]
        key <- parts[1]
        value <- paste(parts[-1], collapse = " ")
        value <- gsub('^"|"$', '', value)
        param_data[[key]] <- value
      }
      current_params(param_data)
      
      params_path <- file.path(data_dir(), paste0(input$seqID, "_params.txt"))
      params_path_global(params_path)
      writeLines(unlist(params), params_path)
      
      output$status <- renderText({
        paste("Params file saved to:", params_path, "\n\nContents:\n", paste(readLines(params_path), collapse = "\n"))
      })
      
      output$download_params <- downloadHandler(
        filename = function() paste0(input$seqID, "_params.txt"),
        content = function(file) {
          req(params_path_global())
          file.copy(params_path_global(), file)
        }
      )
      params_generated(TRUE)
      
    }, error = function(e) {
      output$status <- renderText(paste("Error generating params:", e$message))
    })
  })
  
  
  observeEvent(input$load_params, {
    req(input$load_params)
    tryCatch({
      sample_sheet_path(NULL)
      param_lines <- readLines(input$load_params$datapath)
      param_list <- list()
      
      for (line in param_lines) {
        if (grepl("^--", line)) {
          parts <- strsplit(line, " ", fixed = TRUE)[[1]]
          key <- gsub("^--", "", parts[1])
          value <- paste(parts[-1], collapse = " ")
          value <- gsub('^"|"$', '', value)
          param_list[[key]] <- value
        }
      }
      
      current_params(param_list)
      
      # UI updates (unchanged)
      updateTextInput(session, "report_title", value = param_list$report_title)
      updateTextInput(session, "group_var", value = param_list$group_var)
      updateTextInput(session, "contrasts", value = param_list$contrasts)
      updateSelectInput(session, "organism", selected = param_list$organism)
      updateSelectInput(session, "annotation_db", selected = param_list$annotation_db)
      updateTextInput(session, "seqID", value = param_list$seqID)
      updateTextInput(session, "sample_sheet", value = param_list$sample_sheet)
      updateTextInput(session, "PI", value = param_list$PI)
      updateTextInput(session, "Institution", value = param_list$Institution)
      updateTextInput(session, "Department", value = param_list$Department)
      updateTextInput(session, "Study_Contact", value = param_list$Study_Contact)
      updateTextInput(session, "Project_Title", value = param_list$Project_Title)
      updateTextAreaInput(session, "Study_Summary", value = param_list$Study_Summary)
      updateTextInput(session, "Sample_Types", value = param_list$Sample_Types)
      updateTextInput(session, "Organism_Text", value = param_list$Organism)
      updateTextInput(session, "Analysis_Goals", value = param_list$Analysis_Goals)
      updateTextInput(session, "Report_Prepared_By", value = param_list$Report_Prepared_By)
      updateTextInput(session, "Report_Reviewed_By", value = param_list$Report_Reviewed_By)
      updateTextInput(session, "raw_seq_URL", value = param_list$raw_seq_URL)
      updateTextInput(session, "multiqc_url", value = param_list$multiqc_url)
      updateTextInput(session, "nfcore_spikein_dir", value = param_list$nfcore_spikein_dir)
      updateNumericInput(session, "min_count_for_filtering", value = as.numeric(param_list$min_count_for_filtering))
      updateNumericInput(session, "min_prop_for_filtering", value = as.numeric(param_list$min_prop_for_filtering))
      updateCheckboxInput(session, "test_mode", value = tolower(param_list$test_mode) == "true")
      
      output$status <- renderText("Params file loaded successfully.")
      
    }, error = function(e) {
      output$status <- renderText(paste("Error loading params file:", e$message))
    })
  })
  observe({
    if (report_ready() && params_generated()) {
      shinyjs::enable("run_report")
    } else {
      shinyjs::disable("run_report")
    }
  })
  observe({
    if (report_ready()) {
      shinyjs::enable("generate_params")
    } else {
      shinyjs::disable("generate_params")
    }
  })
  
  
  
  observeEvent(input$run_report, {
    output$status <- renderText("Running report...")
    withProgress(message = "Rendering report...", value = 0.1, {
      tryCatch({
        params_path <- params_path_global()
        if (is.null(params_path) || !file.exists(params_path)) {
          stop("Params file not found. Please click 'Generate Params File' first.")
        }
        report_path <- file.path(data_dir(), paste0(input$seqID, ".Report.html"))
        
        incProgress(0.6, detail = "Rendering R Markdown...")
        rmarkdown::render(
          input = "atac_edgeR_report.Rmd",
          output_file = report_path,
          params = list(
            params_file = params_path,
            report_title = input$report_title,
            test_mode = input$test_mode
          ),
          envir = new.env(parent = globalenv()),
          quiet = FALSE
        )
        
        output$download_report <- downloadHandler(
          filename = function() paste0(input$seqID, ".Report.html"),
          content = function(file) {
            file.copy(report_path, file, overwrite = TRUE)
          }
        )
        
        incProgress(0.1)
        output$status <- renderText({
          paste("Rendered report:", report_path)
        })
        
      }, error = function(e) {
        output$status <- renderText(paste("Error rendering report:", e$message))
      })
    })
  })
  observeEvent(input$clear_params, {
    # Reset all inputs
    updateTextInput(session, "report_title", value = "")
    updateTextInput(session, "nfcore_spikein_dir", value = "")
    updateTextInput(session, "group_var", value = "")
    updateTextInput(session, "contrasts", value = "")
    updateTextInput(session, "seqID", value = "")
    updateTextInput(session, "sample_sheet", value = "")
    
    updateSelectInput(session, "organism", selected = "mmu")
    updateSelectInput(session, "annotation_db", selected = "org.Mm.eg.db")
    
    updateNumericInput(session, "min_count_for_filtering", value = 10)
    updateNumericInput(session, "min_prop_for_filtering", value = 0.5)
    updateCheckboxInput(session, "test_mode", value = TRUE)
    
    # Optional metadata fields
    updateTextInput(session, "PI", value = "")
    updateTextInput(session, "Institution", value = "")
    updateTextInput(session, "Department", value = "")
    updateTextInput(session, "Study_Contact", value = "")
    updateTextInput(session, "Project_Title", value = "")
    updateTextAreaInput(session, "Study_Summary", value = "")
    updateTextInput(session, "Sample_Types", value = "")
    updateTextInput(session, "Organism_Text", value = "")
    updateTextInput(session, "Analysis_Goals", value = "")
    updateTextInput(session, "Report_Prepared_By", value = "")
    updateTextInput(session, "Report_Reviewed_By", value = "")
    updateTextInput(session, "raw_seq_URL", value = "")
    updateTextInput(session, "multiqc_url", value = "")
    
    # Reset reactive flags
    report_ready(FALSE)
    params_generated(FALSE)
    
    # Disable buttons and clear messages
    shinyjs::disable("run_report")
    output$sample_validation_status <- renderText("")
    output$status <- renderText("")
    sample_sheet_path(NULL)
    sample_sheet_temp_file(NULL)
    params_path_global(NULL)
  })
  
}
