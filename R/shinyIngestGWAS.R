library(shiny)
library(data.table)
library(readxl)

#' @export
shinyIngestGWAS <- function() {
 shiny::runApp(system.file("shiny/ingest_app", package = "gact"))
}

canon_fields <- c("marker","chr","pos","ea","nea","eaf","b","seb","p","n","ncase","ncontrol","info")

ui <- fluidPage(
 titlePanel("GWAS Ingestion App"),

 sidebarLayout(
  sidebarPanel(

   # =======================
   # DATABASE
   # =======================
   h4("Database"),
   fileInput("db_file", "Select GAlist (.rds)", accept = ".rds"),
   textOutput("db_status"),

   hr(),

   # =======================
   # GWAS FILES
   # =======================
   fileInput("files", "Upload GWAS files", multiple = TRUE),

   hr(),

   # =======================
   # MANIFEST (NEW)
   # =======================
   h4("Manifest (batch ingest)"),
   fileInput("manifest", "Upload Excel manifest", accept = c(".xlsx", ".xls")),
   actionButton("load_manifest", "Load manifest"),
   uiOutput("manifest_selector"),
   actionButton("ingest_manifest", "Ingest selected traits"),

   hr(),

   # =======================
   # TEMPLATES
   # =======================
   h4("Metadata Templates"),
   selectInput("template", "Template",
               choices = c("None", "UKB", "FinnGen", "SAIGE"),
               selected = "None"),
   checkboxInput("auto_template", "Auto-detect template", TRUE),
   textOutput("template_status"),

   hr(),

   # =======================
   # METADATA
   # =======================
   h4("Metadata"),
   textInput("source", "Source"),
   textInput("trait", "Trait"),
   selectInput("type", "Type", c("quantitative", "binary")),
   textInput("ancestry", "Ancestry"),
   textInput("build", "Build", "GRCh37"),
   numericInput("n", "Sample size", value = NA),
   numericInput("ncase", "Cases", value = 0),
   numericInput("ncontrol", "Controls", value = 0),
   textInput("reference", "Reference"),
   textInput("comments", "Comments"),

   hr(),
   actionButton("detect", "Detect schema"),
   actionButton("ingest_one", "Ingest current file"),
   actionButton("ingest_all", "Ingest ALL files")
  ),

  mainPanel(
   selectInput("file_select", "Select file", choices = NULL),

   tabsetPanel(
    tabPanel("Preview", tableOutput("preview")),
    tabPanel("Schema", verbatimTextOutput("schema")),
    tabPanel("Mapping Editor", uiOutput("mapping_ui")),
    tabPanel("Validation", verbatimTextOutput("validation")),
    tabPanel("Transformations", verbatimTextOutput("transformations")),
    tabPanel("Manifest", tableOutput("manifest_preview"))
   )
  )
 )
)

server <- function(input, output, session) {

 options(shiny.maxRequestSize = 5000 * 1024^2)

 # =======================
 # STORAGE
 # =======================
 GAlist <- reactiveVal(NULL)
 db_path <- reactiveVal(NULL)
 datasets <- reactiveVal(list())
 schemas <- reactiveVal(list())
 manifest_data <- reactiveVal(NULL)

 # =======================
 # LOAD DB
 # =======================
 observeEvent(input$db_file, {
  req(input$db_file)

  tmp <- input$db_file$datapath
  stable <- file.path(tempdir(), input$db_file$name)
  file.copy(tmp, stable, overwrite = TRUE)

  GAlist(readRDS(stable))
  db_path(stable)

  output$db_status <- renderText({
   paste("Loaded DB:", input$db_file$name)
  })
 })

 # =======================
 # TEMPLATE LOGIC
 # =======================
 detect_template_from_data <- function(df, fname = "") {
  cols <- tolower(colnames(df))
  if (any(grepl("mlogp|neglog10", cols))) return("FinnGen")
  if ("#chrom" %in% cols || "markername" %in% cols) return("UKB")
  if ("markerid" %in% cols || "p.value" %in% cols) return("SAIGE")
  "None"
 }

 apply_template <- function(template) {
  if (template == "UKB") {
   updateTextInput(session, "source", value = "UK Biobank")
   updateTextInput(session, "ancestry", value = "EUR")
   updateTextInput(session, "build", value = "GRCh37")
  }
  if (template == "FinnGen") {
   updateTextInput(session, "source", value = "FinnGen")
   updateTextInput(session, "ancestry", value = "FIN")
   updateTextInput(session, "build", value = "GRCh38")
   updateSelectInput(session, "type", selected = "binary")
  }
  if (template == "SAIGE") {
   updateTextInput(session, "source", value = "SAIGE")
   updateTextInput(session, "build", value = "GRCh38")
  }
 }

 observeEvent(input$template, {
  if (input$template != "None") apply_template(input$template)
 })

 output$template_status <- renderText({
  paste("Template:", input$template)
 })

 # =======================
 # LOAD GWAS FILES
 # =======================
 observeEvent(input$files, {
  data_list <- lapply(input$files$datapath, fread, data.table = FALSE)
  names(data_list) <- input$files$name
  datasets(data_list)

  updateSelectInput(session, "file_select",
                    choices = names(data_list),
                    selected = names(data_list)[1])

  if (isTRUE(input$auto_template)) {
   tmpl <- detect_template_from_data(data_list[[1]])
   updateSelectInput(session, "template", selected = tmpl)
  }
 })

 current_data <- reactive({
  req(input$file_select)
  datasets()[[input$file_select]]
 })

 output$preview <- renderTable(head(current_data()))

 # =======================
 # MANIFEST LOAD
 # =======================
 observeEvent(input$load_manifest, {
  req(input$manifest)

  df <- read_excel(input$manifest$datapath)
  manifest_data(df)

  output$manifest_preview <- renderTable(head(df))

  showNotification("Manifest loaded")
 })

 output$manifest_selector <- renderUI({
  req(manifest_data())

  selectInput("manifest_rows",
              "Select traits",
              choices = manifest_data()$phenocode,
              multiple = TRUE)
 })

 # =======================
 # AUTO-FILL FROM MANIFEST
 # =======================
 observeEvent(input$manifest_rows, {
  req(manifest_data())

  row <- manifest_data()[manifest_data()$phenocode == input$manifest_rows[1], ]

  if (nrow(row) == 0) return()

  updateTextInput(session, "trait", value = row$Description)
  updateTextInput(session, "source", value = "UK Biobank")
  updateTextInput(session, "ancestry", value = "EUR")
  updateNumericInput(session, "n", value = row$n_cases_full_cohort_both_sexes)
 })

 # =======================
 # SCHEMA
 # =======================
 observeEvent(input$detect, {
  sch <- detectStatSchema(current_data())
  s <- schemas()
  s[[input$file_select]] <- sch
  schemas(s)
 })

 current_schema <- reactive({
  schemas()[[input$file_select]]
 })

 output$schema <- renderPrint({
  req(current_schema())
  current_schema()[c("mapping","confidence","warnings")]
 })

 # =======================
 # MAPPING
 # =======================
 output$mapping_ui <- renderUI({
  cols <- colnames(current_data())
  sch <- current_schema()

  lapply(canon_fields, function(f) {
   selected <- if (!is.null(sch) && !is.na(sch$mapping[[f]]))
    sch$mapping[[f]] else ""

   selectInput(paste0("map_", f), f, c("", cols), selected = selected)
  })
 })

 custom_mapping <- reactive({
  vals <- sapply(canon_fields, function(f) input[[paste0("map_", f)]], simplify = FALSE)
  vals <- lapply(vals, function(x) if (is.null(x) || x == "") NA else x)
  setNames(unlist(vals), canon_fields)
 })

 normalized <- reactive({
  normalizeStatSchema(current_data(), custom_mapping())
 })

 validation_res <- reactive({
  validateStatSchema(normalized())
 })

 output$validation <- renderPrint(validation_res())

 output$transformations <- renderPrint({
  attr(normalized(), "transformations")
 })

 # =======================
 # INGEST MANIFEST (🔥)
 # =======================
 observeEvent(input$ingest_manifest, {

  req(GAlist(), db_path(), manifest_data(), input$manifest_rows)

  df <- manifest_data()
  G <- GAlist()

  selected <- df[df$phenocode %in% input$manifest_rows, ]

  for (i in seq_len(nrow(selected))) {

   row <- selected[i, ]
   tmp <- tempfile(fileext = ".tsv.gz")

   tryCatch({

    download.file(row$aws_link, tmp, mode = "wb")

    dat <- fread(tmp)

    sch <- detectStatSchema(dat)
    dat_norm <- normalizeStatSchema(dat, sch)

    if (!validateStatSchema(dat_norm)$ok) next

    G <- updateStatDB(
     GAlist = G,
     stat = dat_norm,
     source = "UKB",
     trait = row$Description,
     type = "quantitative",
     gender = "both",
     ancestry = "EUR",
     build = "GRCh37",
     reference = row$trait_efo_terms,
     n = row$n_cases_full_cohort_both_sexes,
     comments = paste("phenocode:", row$phenocode),
     writeStatDB = TRUE
    )

   }, error = function(e) {
    showNotification(paste("Failed:", row$phenocode), type = "error")
   })
  }

  GAlist(G)
  saveRDS(G, db_path(), compress = FALSE)

  showNotification("Manifest ingestion complete!")
 })

}

shinyApp(ui, server)


