library(shiny)
library(data.table)
library(DT)
library(gact)
library(jsonlite)

canon_fields <- c(
 "marker","chr","pos","ea","nea","eaf",
 "b","seb","p","n","ncase","ncontrol","info"
)

shinyApp(

 ui = fluidPage(
  titlePanel("GWAS Ingestion Pipeline"),

  sidebarLayout(
   sidebarPanel(

    fileInput("db_file", "Select GAlist (.rds)", accept = ".rds"),
    textOutput("db_status"),

    hr(),

    fileInput(
     "files", "Upload GWAS files",
     multiple = TRUE,
     accept = c(".txt", ".tsv", ".csv", ".gz")
    ),

    selectInput("file_select", "Select uploaded file", choices = NULL),

    hr(),

    actionButton("detect", "Auto-detect schema"),
    actionButton("process", "Preview normalization"),
    actionButton("save_recipe", "Save recipe"),
    actionButton("ingest", "Ingest into DB"),

    hr(),

    h4("Metadata"),

    textInput("trait", "Trait"),
    selectInput("type", "Type", choices = c("quantitative", "binary")),
    textInput("category", "Category"),
    textInput("efo", "EFO"),

    numericInput("n", "N", NA),
    numericInput("ncase", "N cases", NA),
    numericInput("ncontrol", "N controls", NA),

    textInput("ancestry", "Ancestry"),
    selectInput("build", "Genome build", choices = c("GRCh37", "GRCh38")),
    selectInput("gender", "Gender", choices = c("both", "males", "females"), selected = "both"),

    textInput("source", "Source"),
    textInput("reference", "Reference"),
    textInput("comments", "Comments")
   ),

   mainPanel(
    tabsetPanel(
     tabPanel("Preview", DTOutput("preview")),
     tabPanel("Schema", verbatimTextOutput("schema")),
     tabPanel("Mapping Editor", uiOutput("mapping_ui")),
     tabPanel("Mapping Summary", DTOutput("mapping_summary")),  # ⭐ NEW
     tabPanel("Normalized", DTOutput("normalized")),
     tabPanel("Validation", verbatimTextOutput("validation")),
     tabPanel("Transformations", verbatimTextOutput("transformations"))
    )
   )
  )
 ),

 server = function(input, output, session) {

  options(shiny.maxRequestSize = 2000 * 1024^2)

  GAlist_rv <- reactiveVal(NULL)
  db_path_rv <- reactiveVal(NULL)
  schema_rv <- reactiveVal(NULL)

  # =========================
  # Load DB
  # =========================
  observeEvent(input$db_file, {
   req(input$db_file)

   tmp <- input$db_file$datapath
   stable <- file.path(tempdir(), input$db_file$name)
   file.copy(tmp, stable, overwrite = TRUE)

   GAlist_rv(readRDS(stable))
   db_path_rv(stable)

   output$db_status <- renderText({
    paste("Loaded DB:", input$db_file$name)
   })
  })

  # =========================
  # FAST LOAD (only 50 rows)
  # =========================
  data_list <- reactive({
   req(input$files)

   out <- lapply(input$files$datapath, function(f) {
    fread(f, data.table = FALSE, nrows = 50)  # ⭐ FAST
   })
   names(out) <- input$files$name
   out
  })

  observeEvent(input$files, {
   req(data_list())

   updateSelectInput(
    session,
    "file_select",
    choices = names(data_list()),
    selected = names(data_list())[1]
   )

   schema_rv(NULL)
  })

  current_data <- reactive({
   req(data_list(), input$file_select)
   data_list()[[input$file_select]]
  })

  output$preview <- renderDT({
   req(current_data())
   datatable(current_data(), options = list(scrollX = TRUE))
  })

  # =========================
  # Schema detection
  # =========================
  observeEvent(input$detect, {
   req(current_data())
   schema_rv(detectStatSchema(current_data()))
  })

  output$schema <- renderPrint({
   req(schema_rv())
   list(
    mapping = schema_rv()$mapping,
    confidence = schema_rv()$confidence,
    warnings = schema_rv()$warnings
   )
  })

  # =========================
  # Mapping UI
  # =========================
  output$mapping_ui <- renderUI({
   req(current_data())

   cols <- colnames(current_data())
   sch <- schema_rv()

   lapply(canon_fields, function(f) {

    selected <- if (!is.null(sch) && !is.na(sch$mapping[[f]])) {
     sch$mapping[[f]]
    } else {
     "IGNORE"
    }

    selectInput(
     inputId = paste0("map_", f),
     label = f,
     choices = c("IGNORE", cols),
     selected = selected
    )
   })
  })

  mapping <- reactive({
   m <- sapply(canon_fields, function(f) {
    val <- input[[paste0("map_", f)]]
    if (is.null(val) || identical(val, "IGNORE")) NA_character_ else val
   }, USE.NAMES = TRUE)

   setNames(m, canon_fields)
  })

  # =========================
  # ⭐ Mapping summary (NEW)
  # =========================
  output$mapping_summary <- renderDT({
   req(mapping())

   df <- data.frame(
    canonical = names(mapping()),
    mapped_to = as.character(mapping()),
    stringsAsFactors = FALSE
   )

   df$mapped_to[is.na(df$mapped_to)] <- "NOT MAPPED"

   datatable(
    df,
    rownames = FALSE,
    options = list(scrollX = TRUE),
    callback = JS(
     "function(row, data) {",
     " if (data[1] === 'NOT MAPPED') {",
     "   $('td:eq(1)', row).css('color', 'red');",
     " }",
     "}"
    )
   )
  })

  # =========================
  # Normalize
  # =========================
  norm_data <- reactive({
   req(current_data(), mapping())

   normalizeStatSchema(
    stat = current_data(),
    schema = mapping()
   )
  })

  output$normalized <- renderDT({
   req(norm_data())
   datatable(head(norm_data(), 50), options = list(scrollX = TRUE))
  })

  # =========================
  # Validation
  # =========================
  validation_rv <- reactive({
   req(norm_data())
   validateStatSchema(norm_data())
  })

  output$validation <- renderPrint({
   req(validation_rv())
   validation_rv()
  })

  output$transformations <- renderPrint({
   req(norm_data())
   attr(norm_data(), "transformations")
  })

  # =========================
  # Save recipe
  # =========================
  observeEvent(input$save_recipe, {

   req(current_data(), mapping(), GAlist_rv())

   dirs <- GAlist_rv()$dirs

   if (!"recipes" %in% names(dirs)) {
    db_root <- dirname(dirs[1])
    recipes_dir <- file.path(db_root, "recipes")
    dir.create(recipes_dir, recursive = TRUE, showWarnings = FALSE)

    G <- GAlist_rv()
    G$dirs["recipes"] <- recipes_dir
    GAlist_rv(G)

   } else {
    recipes_dir <- dirs["recipes"]
   }

   clean_field <- function(x, default = NA) {
    if (is.null(x) || x == "") return(default)
    x
   }

   recipe <- list(
    file_name = input$file_select,

    trait = input$trait,
    type = input$type,
    category = clean_field(input$category, "unknown"),
    efo = clean_field(input$efo, NA),

    gender = input$gender,
    ancestry = clean_field(input$ancestry, "unknown"),
    build = input$build,

    n = input$n,
    ncase = input$ncase,
    ncontrol = input$ncontrol,

    source = clean_field(input$source, basename(input$file_select)),
    reference = clean_field(input$reference, NA),
    comments = clean_field(input$comments, ""),

    mapping = as.list(mapping()),
    transformations = attr(norm_data(), "transformations"),
    validation = validation_rv(),

    created_at = as.character(Sys.time())
   )

   base_name <- tools::file_path_sans_ext(basename(input$file_select))

   recipe_file <- file.path(
    recipes_dir,
    paste0(base_name, "_recipe.json")
   )

   write_json(recipe, recipe_file, pretty = TRUE, auto_unbox = TRUE)

   showNotification(paste("Recipe saved:", recipe_file), duration = 6)
  })

  # =========================
  # Ingest
  # =========================
  observeEvent(input$ingest, {

   req(norm_data(), validation_rv(), GAlist_rv(), db_path_rv())

   if (!isTRUE(validation_rv()$ok)) {
    showNotification("Validation failed", type = "error")
    return()
   }

   G_updated <- updateStatDB(
    GAlist = GAlist_rv(),
    stat = norm_data(),
    source = input$source,
    trait = input$trait,
    type = input$type,
    category = input$category,
    efo = input$efo,
    gender = input$gender,
    ancestry = input$ancestry,
    build = input$build,
    n = input$n,
    ncase = input$ncase,
    ncontrol = input$ncontrol,
    reference = input$reference,
    comments = input$comments,
    writeStatDB = TRUE
   )

   GAlist_rv(G_updated)
   saveRDS(G_updated, db_path_rv(), compress = FALSE)

   showNotification("Ingestion successful", type = "message")
  })
 }
)
