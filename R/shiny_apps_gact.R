#' @export
#'

shinyAppsDB <- function(GAlist=NULL, what="DrugDBTables") {

 if(what=="DrugDBTables") {
  require(shiny)
  require(data.table)
  require(DT)

  shinyApp(
   ui <- fluidPage(
    sidebarLayout(
     sidebarPanel(
      selectInput("what", "Select a file to load:",
                  choices = c("Drug Labels PharmGKB",
                              "Drug Labels Genes PharmGKB",
                              "Variants PharmGKB",
                              "Automated Annotations PharmGKB",
                              "Drugs DGI",
                              "Genes DGI",
                              "Categories DGI",
                              "Interactions DGI",
                              "Direct Associations OpenTarget",
                              "Direct Associations By Data Type OpenTarget",
                              "Direct Associations By Data Source OpenTarget"))
     ),
     mainPanel(
      DTOutput("data")
     )
    )
   )
   ,
   server <- function(input, output) {

    observeEvent(input$what, {
     if (input$what == "Drug Labels PharmGKB") {
      drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "drugLabels.tsv")
     } else if (input$what == "Drug Labels Genes PharmGKB") {
      drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "drugLabels.byGene.tsv")
     } else if (input$what == "Variants PharmGKB") {
      drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "clinicalVariants.tsv")
     } else if (input$what == "Automated Annotations PharmGKB") {
      drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "automated_annotations.tsv")
     } else if (input$what == "Drugs DGI") {
      drugdb_file <- file.path(GAlist$dirs["dgidb"], "drugs.tsv")
     } else if (input$what == "Genes DGI") {
      drugdb_file <- file.path(GAlist$dirs["dgidb"], "genes.tsv")
     } else if (input$what == "Categories DGI") {
      drugdb_file <- file.path(GAlist$dirs["dgidb"], "categories.tsv")
     } else if (input$what == "Interactions DGI") {
      drugdb_file <- file.path(GAlist$dirs["dgidb"], "interactions.tsv")
     } else if (input$what == "Direct Associations OpenTarget") {
      drugdb_file <- file.path(GAlist$dirs["opentargets"], "associationByOverallDirect.tsv")
     } else if (input$what == "Direct Associations By Data Type OpenTarget") {
      drugdb_file <- file.path(GAlist$dirs["opentargets"], "associationByDatatypeDirect.tsv")
     } else if (input$what == "Direct Associations By Data Source OpenTarget") {
      drugdb_file <- file.path(GAlist$dirs["opentargets"], "associationByDatasourceDirect.tsv")
     }
     df <- fread(drugdb_file, quote = "", data.table = FALSE)
     output$data <- renderDT(df,
                             options = list(pageLength = 50, dom = 'Bfrtip',
                                            buttons = c('copy', 'csv')),
                             extensions = "Buttons")
    })
   }

  )
 }
}
