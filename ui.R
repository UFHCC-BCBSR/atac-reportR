library(bs4Dash)
library(shinyjs)


ui <- dashboardPage(
  
  header = dashboardHeader(
    title = dashboardBrand(
      title = strong("UFHCC BCB-SR"),
      color = "lightblue",
      href = "https://cancer.ufl.edu/research/shared-resources/biostatistics-computational-biology-shared-resource/",
      image = "https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/UFH_CancerCenter.jpg"
    ),
    skin = "light",
    status = "white"
  ),
  
  sidebar = dashboardSidebar(
    skin = "light",
    status = "primary",
    width = 300,
    sidebarMenu(
      id = "sidebarMenu",
      menuItem("Report Generator", tabName = "main", icon = icon("cogs")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  )
  ,
  
  body = dashboardBody(
    useShinyjs(),
    tabItems(
      tabItem(tabName = "main",
              fluidRow(
                box(
                  title = "ATAC-seq Report Generator",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  
                  fileInput("load_params", "Upload Existing Params File", accept = ".txt"),
                  actionButton("clear_params", "Reset Params", class = "btn btn-outline-danger btn-block"),
                  tags$h4("Analysis Parameters"),
                  helpText("These fields are required to generate the report. Please ensure contrast groups and project ID are accurate."),
                  textInput("seqID", "Project ID", value = "NS4167_Male_GoodSamples"),
                  helpText("Project ID must exactly match the directory in /blue where your atac-seq results are",
                           "See the 'About' page for more details"),
                  textInput("report_title", "Report Title", value = "Differential Accessibility of Heat vs. Control in Mouse (QC-passed samples only)"),
                  textInput("group_var", "Group Variable (must not be 'group')", value = "Condition"),
                  # In your UI layout
                  fileInput("contrasts_file", "Upload Contrasts File (optional)", accept = ".csv"),
                  textInput("contrasts_text", "Or Enter Contrasts (comma-separated)", 
                            placeholder = "e.g., Group1_vs_Group2, Var1_vs_Var2"),
                  helpText("You may either upload a plain text `.csv` file (no header, one contrast per line), ",
                           "OR enter a comma-separated list of contrasts in the text box above. ",
                           "The app will automatically detect the variable(s) used in each contrast."),
                  textInput("nfcore_spikein_dir", "Spike-in Dir (Optional. If not specified, TMM normalization will be performed.)", value = ""),
                  textOutput("sample_sheet_source"),
                  downloadButton("download_sample_sheet", "Download and Edit Sample Information (Optional)", class = "btn-secondary btn-sm"),
                  fileInput("upload_sample_sheet", "Upload and Save Edited Sample Sheet (Optional)", accept = ".csv"),
                  actionButton("validate_sheet", "Validate Sample Sheet"),
                  verbatimTextOutput("sample_validation_status"),
                  
                  selectInput("organism", "Organism", choices = c("mmu", "hsa"), selected = "mmu"),
                  selectInput("annotation_db", "Annotation DB", choices = c("org.Mm.eg.db", "org.Hs.eg.db"), selected = "org.Mm.eg.db"),
                  numericInput("min_count_for_filtering", "Min Count for Filtering", value = 10, min = 0),
                  numericInput("min_prop_for_filtering", "Min Prop for Filtering", value = 0.5, min = 0, max = 1, step = 0.05),
                  
                  tags$h4("Optional Parameters"),
                  helpText("Optional parameters to add information to the report. The report can be run without these."),
                  
                  textInput("PI", "Principal Investigator", value = "[Dr. Thomas Clanton](mailto:tclanton@hhp.ufl.edu)"),
                  textInput("Institution", "Institution", value = "University of Florida"),
                  textInput("Department", "Department", value = "UF Health Cancer Center"),
                  textInput("Study_Contact", "Study Contact Email", value = "[Jason Brant](mailto:jobrant@ufl.edu)"),
                  textInput("Project_Title", "Project Title", value = "Accessible Epigenetic biomarkers for climate Readiness from Mouse to Human"),
                  textAreaInput("Study_Summary", "Study Summary",
                                value = "This project will lay the biological groundwork for developing a rapid test from human blood cells that accurately reports whether a service member is sufficiently acclimatized to high heat and humidity...",
                                height = "100px"),
                  textInput("Sample_Types", "Sample Types", value = "Mouse (male) lymphocytes"),
                  textInput("Organism_Text", "Organism (For display)", value = "Mus musculus"),
                  textInput("Analysis_Goals", "Analysis Goals", value = ""),
                  textInput("Report_Prepared_By", "Report Prepared By", value = "[Dr. Heather Kates](mailto:hkates@ufl.edu), Bioinformatics Analyst"),
                  textInput("Report_Reviewed_By", "Report Reviewed By", value = "[Dr. Jason Brant](mailto:jobrant@ufl.edu), Assistant Professor"),
                  textInput("raw_seq_URL", "Raw Sequencing URL", value = ""),
                  textInput("multiqc_url", "MultiQC URL", value = ""),
                  checkboxInput("test_mode", "Test Mode (Only render setup chunks)", value = TRUE),
                  
                  
                  tags$hr(),
                  actionButton("generate_params", "Generate Params File", class = "btn-secondary btn-block",disabled=TRUE),
                  actionButton("run_report", "Run Report", class = "btn btn-primary btn-block",disabled=TRUE),
                  downloadButton("download_report", "Download Report"),
                  downloadButton("download_params", "Download Params File"),
                  
                  br(), verbatimTextOutput("status"),
                  uiOutput("download_link")
                )
              )
      )
      ,
      tabItem(tabName = "about",
              box(
                title = "About This Tool",
                status = "info",
                width = 12,
                solidHeader = TRUE,
                h3("What is ATAC-seq Differential Accessibility Analysis?"),
                p("ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) identifies regions of open chromatin in the genome. 
      Differential accessibility (DA) analysis compares chromatin openness across conditions, such as treatment vs. control, to identify regions that may be associated with regulatory activity."),
                
                h3("About This Report Generator"),
                p("This application allows users to create a differential accessibility analysis report from nf-core/atacseq output. 
      It supports optional spike-in normalization, generates contrast-specific results, and creates a downloadable HTML report."),
                
                h3("How to Use"),
                tags$ol(
                  tags$li(
                    tagList(
                      "Navigate to ",
                      tags$code("/blue/cancercenter-dept/privapps/data/atac"),
                      " and run:",
                      tags$br(),
                      tags$code("bash retrieve-atac-results.bash --output <your_nf-core/atacseq outdir> --destination /blue/cancercenter-dept/privapps/data/atac/<your_PI>/<Project ID>"),
                      tags$br(),
                      "to copy your nf-core/atacseq output to the directory accessible to this app.",
                      tags$br(),
                      "Project ID must exactly match the Project ID value entered in the app (PI can be anything, but it must exist; Project ID must be a subdirectory of PI)"
                    )
                  ),
                  tags$li("Upload a pre-configured parameters file or manually fill out the form."),
                  tags$li("Optionally download, edit, and re-upload a sample sheet to customize sample names or replicates."),
                  tags$li("Validate your sample sheet to ensure contrasts can be tested."),
                  tags$li("Click 'Generate Params File' and then 'Run Report'."),
                  tags$li("Download the rendered HTML report when complete.")
                ),
                
                h3("Need Help?"),
                p("Contact the UFHCC BCB-SR Bioinformatics team at ",
                  tags$a(href = "mailto:hkates@ufl.edu", "hkates@ufl.edu"))
              )
      )
      
    )
  ),
  
  footer = dashboardFooter(
    left = "Contact: hkates@ufl.edu",
    right = "UFHCC BCB-SR",
    fixed = FALSE
  )
)
