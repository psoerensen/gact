#' @title Launch the installed GWAS ingestion application
#'
#' @description Runs the packaged Shiny application located at inst/shiny/ingest_app.
#'
#'
#' @details The app is used to select local files and ingest them into a GACT database. Call it explicitly in an interactive session.
#'
#' @return The return value of shiny::runApp(); the call launches an interactive app and blocks until it exits.
#' @export
shinyIngestGWAS <- function() {
 shiny::runApp(system.file("shiny/ingest_app", package = "gact"))
}
